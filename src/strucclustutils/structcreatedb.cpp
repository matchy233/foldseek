//
// Created by Martin Steinegger on 11/14/20.
//
#include "Command.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "DBWriter.h"
#include "FastSort.h"

#include "structureto3di.h"
#include "SubstitutionMatrix.h"
#include "GemmiWrapper.h"
#include "PulchraWrapper.h"
#include "microtar.h"
#include "PatternCompiler.h"
#include "Coordinate16.h"
#include "itoa.h"
#include "MathUtil.h"

#ifdef HAVE_PROSTT5
#include "ProstT5.h"
#if !defined(__CYGWIN__) && !defined(__EMSCRIPTEN__) && !defined(__APPLE__)
#include "ProstT5ForkRunner.h"
#define FORK_RUNNER 1
#else
#define FORK_RUNNER 0
#define prostt5Forking(...)
#endif
#endif

#include <iostream>
#include <dirent.h>
#include <sys/stat.h>

#ifdef HAVE_GCS
#include "google/cloud/storage/client.h"
#endif

#ifdef OPENMP
#include <omp.h>
#endif


#ifdef HAVE_ZLIB
#include <zlib.h>
static int structure_file_gzread(mtar_t *tar, void *data, size_t size) {
    size_t res = gzread((gzFile)tar->stream, data, size);
    return (res == size) ? MTAR_ESUCCESS : MTAR_EREADFAIL;
}

static int structure_file_gzseek(mtar_t *tar, long offset, int whence) {
    int res = gzseek((gzFile)tar->stream, offset, whence);
    return (res != -1) ? MTAR_ESUCCESS : MTAR_ESEEKFAIL;
}

static int structure_file_gzclose(mtar_t *tar) {
    gzclose((gzFile)tar->stream);
    return MTAR_ESUCCESS;
}

int structure_mtar_gzopen(mtar_t *tar, const char *filename) {
    // Init tar struct and functions
    memset(tar, 0, sizeof(*tar));
    tar->read = structure_file_gzread;
    tar->seek = structure_file_gzseek;
    tar->close = structure_file_gzclose;
    tar->isFinished = 0;
    // Open file
    tar->stream = gzopen(filename, "rb");
    if (!tar->stream) {
        return MTAR_EOPENFAIL;
    }

#if defined(ZLIB_VERNUM) && ZLIB_VERNUM >= 0x1240
    if (gzbuffer((gzFile)tar->stream, 1 * 1024 * 1024) != 0) {
        Debug(Debug::WARNING) << "Could not set gzbuffer size, performance might be bad\n";
    }
#endif

    return MTAR_ESUCCESS;
}
#endif

#define ZSTD_STATIC_LINKING_ONLY
#include <zstd.h>

struct zstd_file_t {
    FILE *fp;
    ZSTD_DStream *dctx;
    char *inBuf;
    size_t inCap;
    size_t inPos;
    size_t inSize;
    char *outBuf;
    size_t outCap;
    size_t outPos;
    size_t outSize;
};

static int file_zstdread(mtar_t *tar, void *data, size_t size) {
    zstd_file_t *s = static_cast<zstd_file_t*>(tar->stream);
    char *dst = static_cast<char*>(data);
    size_t need = size;

    while (need) {
        size_t avail = s->outSize - s->outPos;
        if (avail) {
            size_t take = (avail < need) ? avail : need;
            memcpy(dst, s->outBuf + s->outPos, take);
            s->outPos += take;
            dst += take;
            need -= take;
            continue;
        }

        if (s->inPos == s->inSize) {
            s->inSize = fread(s->inBuf, 1, s->inCap, s->fp);
            s->inPos = 0;
            if (s->inSize == 0) {
                break;
            }
        }

        ZSTD_inBuffer  in  { s->inBuf + s->inPos, s->inSize - s->inPos, 0 };
        ZSTD_outBuffer out { s->outBuf, s->outCap, 0 };

        const size_t ret = ZSTD_decompressStream(s->dctx, &out, &in);
        if (ZSTD_isError(ret)) {
            return MTAR_EREADFAIL;
        }

        s->inPos  += in.pos;
        s->outSize = out.pos;
        s->outPos  = 0;

        if (ret == 0 && s->outSize == 0) {
            break;
        }
    }
    return (need == 0) ? MTAR_ESUCCESS : MTAR_EREADFAIL;
}

static int file_zstdseek(mtar_t *tar, long offset, int whence) {
    if (whence != SEEK_CUR || offset < 0) {
        return MTAR_ESEEKFAIL;
    }

    char tmp[8192];
    long left = offset;
    while (left) {
        size_t chunk = (left > (long)sizeof(tmp)) ? sizeof(tmp) : (size_t)left;
        if (file_zstdread(tar, tmp, chunk) != MTAR_ESUCCESS) {
            return MTAR_ESEEKFAIL;
        }
        left -= (long)chunk;
    }
    return MTAR_ESUCCESS;
}

static int file_zstdclose(mtar_t *tar) {
    zstd_file_t *s = static_cast<zstd_file_t*>(tar->stream);
    if (s) {
        ZSTD_freeDStream(s->dctx);
        fclose(s->fp);
        free(s->inBuf);
        free(s->outBuf);
        delete s;
        tar->stream = nullptr;
    }
    return MTAR_ESUCCESS;
}

int mtar_zstdopen(mtar_t *tar, const char *filename) {
    std::memset(tar, 0, sizeof(*tar));
    tar->read  = file_zstdread;
    tar->seek  = file_zstdseek;
    tar->close = file_zstdclose;
    tar->isFinished = 0;

    zstd_file_t *s = new zstd_file_t{};
    tar->stream = s;

    s->fp = fopen(filename, "rb");
    if (!s->fp) {
        delete s;
        return MTAR_EOPENFAIL;
    }

    s->dctx = ZSTD_createDStream();
    if (!s->dctx) {
        std::fclose(s->fp);
        delete s;
        return MTAR_EOPENFAIL;
    }

    s->inCap  = ZSTD_DStreamInSize();
    s->outCap = ZSTD_DStreamOutSize();
    s->inBuf  = static_cast<char*>(malloc(s->inCap));
    s->outBuf = static_cast<char*>(malloc(s->outCap));
    if (!s->inBuf || !s->outBuf) {
        file_zstdclose(tar);
        return MTAR_EOPENFAIL;
    }

    s->inSize = s->inPos = s->outSize = s->outPos = 0;
    return MTAR_ESUCCESS;
}

template <typename T, typename U>
static inline bool compareByFirst(const std::pair<T, U>& a, const std::pair<T, U>& b) {
    return a.first < b.first;
}

std::string removeModel(const std::string& input) {
    size_t modelIndex = input.find("MODEL");
    if (modelIndex == std::string::npos) {
        return input;
    }
    std::string prefix = input.substr(0, modelIndex);
    size_t secondUnderscoreIndex = input.find('_', modelIndex + 6);
    if (secondUnderscoreIndex == std::string::npos) {
        return prefix;
    }
    std::string suffix = input.substr(secondUnderscoreIndex+1);
    return prefix + suffix;
}

void addMissingAtomsInStructure(GemmiWrapper &readStructure, PulchraWrapper &pulchra) {
    std::vector<std::string> chainNames;
    for (size_t ch = 0; ch < readStructure.chain.size(); ch++) {
        size_t chainStart = readStructure.chain[ch].first;
        size_t chainEnd = readStructure.chain[ch].second;
        size_t chainLen = chainEnd - chainStart;
        bool allX = true;
        for (size_t pos = 0; pos < chainLen; pos++) {
            const char aa = readStructure.ami[chainStart+pos];
            if (aa != 'X' && aa != 'x') {
                allX = false;
                break;
            }
        }
        if (allX) {
            chainNames.push_back("SKIP");
            continue;
        }
        chainNames.push_back(readStructure.chainNames[ch]);
        // Detect if structure is Ca only
        if (std::isnan(readStructure.n[chainStart + 0].x) &&
            std::isnan(readStructure.n[chainStart + 1].x) &&
            std::isnan(readStructure.n[chainStart + 2].x) &&
            std::isnan(readStructure.n[chainStart + 3].x) &&
            std::isnan(readStructure.c[chainStart + 0].x) &&
            std::isnan(readStructure.c[chainStart + 1].x) &&
            std::isnan(readStructure.c[chainStart + 2].x) &&
            std::isnan(readStructure.c[chainStart + 3].x)) {
            pulchra.rebuildBackbone(&readStructure.ca[chainStart],
                                    &readStructure.n[chainStart],
                                    &readStructure.c[chainStart],
                                    &readStructure.ami[chainStart],
                                    chainLen);
        
        }
    }
    readStructure.chainNames.clear();
    readStructure.chainNames.insert(readStructure.chainNames.begin(), chainNames.begin(), chainNames.end());
    chainNames.clear();
}

void findInterfaceResidues(GemmiWrapper &readStructure, std::pair<size_t, size_t> res1, std::pair<size_t, size_t> res2,
                           std::vector<size_t> & resIdx1, float distanceThreshold)
{
    std::vector<Vec3> coord1, coord2;
    size_t sameRes = 0;
    size_t chainLen = res1.second - res1.first;
    bool noSameRes = true;
    const float squareThreshold = distanceThreshold * distanceThreshold;
    for (size_t res1Idx = res1.first; res1Idx < res1.second; res1Idx++) {
        float x1, y1, z1;
        if (readStructure.ami[res1Idx] == 'G') {
            x1 = readStructure.ca[res1Idx].x;
            y1 = readStructure.ca[res1Idx].y;
            z1 = readStructure.ca[res1Idx].z;
        }
        else {
            x1 = readStructure.cb[res1Idx].x;
            y1 = readStructure.cb[res1Idx].y;
            z1 = readStructure.cb[res1Idx].z;
        }
        for (size_t res2Idx = res2.first; res2Idx < res2.second; res2Idx++) {
            float x2, y2, z2;
            if (readStructure.ami[res2Idx] == 'G') {
                x2 = readStructure.ca[res2Idx].x;
                y2 = readStructure.ca[res2Idx].y;
                z2 = readStructure.ca[res2Idx].z;
            }
            else {
                x2 = readStructure.cb[res2Idx].x;
                y2 = readStructure.cb[res2Idx].y;
                z2 = readStructure.cb[res2Idx].z;
            }
            float distance = MathUtil::squareDist(x1, y1, z1, x2, y2, z2);
            if (distance < 0.01) {
                noSameRes = false;
                sameRes++;
                break;
            }
            if (distance < squareThreshold) {
                resIdx1.push_back(res1Idx);
                if (noSameRes) {
                    break;
                }
            } 
        }
    }
    if (sameRes / chainLen > 0.9){
        resIdx1.clear();
    }
}

void compute3DiInterfaces(GemmiWrapper &readStructure, PulchraWrapper &pulchra, StructureTo3Di &structureTo3Di, SubstitutionMatrix & mat3Di, int chainNameMode, float distanceThreshold) {
    size_t prevInterfaceChainLen = 0;
    std::vector<char> interfaceSeq3di, interfaceAmi;
    std::vector<size_t> resIdx1, resIdx2;
    std::vector<std::string> notProteinChains;
    std::vector<int> interfacetaxIds;
    std::vector<unsigned int> interfaceModelIndices;
    std::vector<Vec3> ca, n, c, cb;
    std::vector<Vec3> interfaceCa;
    std::vector<std::string> interfaceNames, interfaceChainNames;
    std::vector<std::pair<size_t, size_t>> interfaceChain;
    std::unordered_map<unsigned int, unsigned int> modelToInterfaceNum;
    addMissingAtomsInStructure(readStructure, pulchra);
    for (size_t ch1 = 0; ch1 < readStructure.chain.size(); ch1++) {
        if (readStructure.chainNames[ch1] == "SKIP") {
            interfaceCa.push_back(Vec3(0,0,0));
            interfaceAmi.push_back('X');
            interfaceSeq3di.push_back('X');
            interfaceChain.push_back(std::make_pair(prevInterfaceChainLen, prevInterfaceChainLen+1));
            interfaceNames.push_back("ALLX");
            interfaceChainNames.push_back(readStructure.chainNames[ch1]);
            interfacetaxIds.push_back(readStructure.taxIds[ch1]);
            interfaceModelIndices.push_back(readStructure.modelIndices[ch1]);
            prevInterfaceChainLen++;
            continue;
        }
        for (size_t ch2 = ch1 + 1; ch2 < readStructure.chain.size(); ch2++) {
            if (readStructure.chainNames[ch2] == "SKIP") {
                continue;
            }
            if (readStructure.modelIndices[ch1] == readStructure.modelIndices[ch2]) {
                findInterfaceResidues(readStructure, readStructure.chain[ch1], readStructure.chain[ch2], resIdx1, distanceThreshold);
                findInterfaceResidues(readStructure, readStructure.chain[ch2], readStructure.chain[ch1], resIdx2, distanceThreshold);
                if (resIdx1.size() >= 4 && resIdx2.size() >= 4) {
                    modelToInterfaceNum[readStructure.modelIndices[ch2]]++;
                    for (size_t i = 0; i < resIdx1.size(); i++) {
                        ca.push_back(readStructure.ca[resIdx1[i]]);
                        n.push_back(readStructure.n[resIdx1[i]]);
                        c.push_back(readStructure.c[resIdx1[i]]);
                        cb.push_back(readStructure.cb[resIdx1[i]]);
                    }
                    // std::sort(resIdx2.begin(), resIdx2.end());
                    for (size_t i = 0; i < resIdx2.size(); i++) {
                        ca.push_back(readStructure.ca[resIdx2[i]]);
                        n.push_back(readStructure.n[resIdx2[i]]);
                        c.push_back(readStructure.c[resIdx2[i]]);
                        cb.push_back(readStructure.cb[resIdx2[i]]);
                    }
                    char *states = structureTo3Di.structure2states(ca.data(),
                                                                n.data(),
                                                                c.data(),
                                                                cb.data(),
                                                                resIdx1.size() + resIdx2.size());
                    for (size_t i = 0; i < resIdx1.size(); i++) {
                        interfaceSeq3di.push_back(mat3Di.num2aa[static_cast<int>(states[i])]);
                        interfaceAmi.push_back(readStructure.ami[resIdx1[i]]);
                        interfaceCa.push_back(readStructure.ca[resIdx1[i]]);
                    }
                    for (size_t i = 0; i < resIdx2.size(); i++) {
                        interfaceSeq3di.push_back(mat3Di.num2aa[static_cast<int>(states[resIdx1.size()+i])]);
                        interfaceAmi.push_back(readStructure.ami[resIdx2[i]]);
                        interfaceCa.push_back(readStructure.ca[resIdx2[i]]);
                    }
                    interfaceChain.push_back(std::make_pair(prevInterfaceChainLen, prevInterfaceChainLen + resIdx1.size()));
                    prevInterfaceChainLen += resIdx1.size();
                    interfaceChain.push_back(std::make_pair(prevInterfaceChainLen, prevInterfaceChainLen + resIdx2.size()));
                    prevInterfaceChainLen += resIdx2.size();
                    std::string interfaceName;
                    if (Util::endsWith(".gz", readStructure.names[ch1]) || Util::endsWith(".zstd", readStructure.names[ch1]) || Util::endsWith(".zst", readStructure.names[ch1])) {
                        interfaceName.append(Util::remove_extension(Util::remove_extension(readStructure.names[ch1])));
                    } else {
                        interfaceName.append(Util::remove_extension(readStructure.names[ch1]));
                    } 
                    if (readStructure.modelCount > 1) {
                        interfaceName.append("_MODEL_");
                        interfaceName.append(SSTR(readStructure.modelIndices[ch1]));
                    }
                    interfaceModelIndices.push_back(readStructure.modelIndices[ch1]);
                    interfaceModelIndices.push_back(readStructure.modelIndices[ch2]);
                    interfaceName.append("_INT_");
                    interfaceName.append(SSTR(modelToInterfaceNum[readStructure.modelIndices[ch2]]));
                    std::string interfaceNameFirst = interfaceName;
                    std::string interfaceNameSecond = interfaceName;
                    if (chainNameMode == LocalParameters::CHAIN_MODE_ADD ||
                        (chainNameMode == LocalParameters::CHAIN_MODE_AUTO && readStructure.names.size() > 1)) {
                        interfaceNameFirst.push_back('_');
                        interfaceNameFirst.append(readStructure.chainNames[ch1]);
                        interfaceNameSecond.push_back('_');
                        interfaceNameSecond.append(readStructure.chainNames[ch2]);
                    }
                    if (readStructure.title.size() > 0) {
                        interfaceNameFirst.push_back(' ');
                        interfaceNameFirst.append(readStructure.title);
                        interfaceNameSecond.push_back(' ');
                        interfaceNameSecond.append(readStructure.title);
                    }   
                    interfaceNames.push_back(interfaceNameFirst);
                    interfaceNames.push_back(interfaceNameSecond);
                    interfaceChainNames.push_back(readStructure.chainNames[ch1]);
                    interfaceChainNames.push_back(readStructure.chainNames[ch2]);
                    interfacetaxIds.push_back(readStructure.taxIds[ch1]);
                    interfacetaxIds.push_back(readStructure.taxIds[ch2]);
                }
                else {
                    interfaceCa.push_back(Vec3(0,0,0));
                    interfaceAmi.push_back('X');
                    interfaceSeq3di.push_back('X');
                    interfaceChain.push_back(std::make_pair(prevInterfaceChainLen, prevInterfaceChainLen+1));
                    interfaceNames.push_back("SHORT");
                    interfaceChainNames.push_back(readStructure.chainNames[ch1]);
                    interfacetaxIds.push_back(readStructure.taxIds[ch1]);
                    interfaceModelIndices.push_back(readStructure.modelIndices[ch1]);
                    prevInterfaceChainLen++;
                }
                resIdx1.clear();
                resIdx2.clear();
                ca.clear();
                n.clear();
                c.clear();
                cb.clear();
            }
        }
    }
    // copy interface data to readStructure
    // this overwrites the original data (clear and push_back)
    readStructure.ca.clear();
    readStructure.ca.insert(readStructure.ca.begin(), interfaceCa.begin(), interfaceCa.end());
    readStructure.ami.clear();
    readStructure.ami.insert(readStructure.ami.begin(), interfaceAmi.begin(), interfaceAmi.end());
    readStructure.seq3di.clear();
    readStructure.seq3di.insert(readStructure.seq3di.begin(), interfaceSeq3di.begin(), interfaceSeq3di.end());
    readStructure.chain.clear();
    readStructure.chain.insert(readStructure.chain.begin(), interfaceChain.begin(), interfaceChain.end());
    readStructure.names.clear();
    readStructure.names.insert(readStructure.names.begin(), interfaceNames.begin(), interfaceNames.end());
    readStructure.chainNames.clear();
    readStructure.chainNames.insert(readStructure.chainNames.begin(), interfaceChainNames.begin(), interfaceChainNames.end());
    readStructure.taxIds.clear();
    readStructure.taxIds.insert(readStructure.taxIds.begin(), interfacetaxIds.begin(), interfacetaxIds.end());
    readStructure.modelIndices.clear();
    readStructure.modelIndices.insert(readStructure.modelIndices.begin(), interfaceModelIndices.begin(), interfaceModelIndices.end());
}

size_t
writeStructureEntry(SubstitutionMatrix & mat, GemmiWrapper & readStructure, StructureTo3Di & structureTo3Di,
                    PulchraWrapper & pulchra, std::vector<char> & alphabet3di, std::vector<char> & alphabetAA,
                    std::vector<int8_t> & camol, std::string & header, 
                    DBWriter & aadbw, DBWriter & hdbw, DBWriter & torsiondbw, DBWriter & cadbw, int chainNameMode,
                    float maskBfactorThreshold, size_t & tooShort, size_t & notProtein, size_t & globalCnt, int thread_idx, int coordStoreMode,
                    size_t & fileidCnt, std::map<std::string, std::pair<size_t, unsigned int>> & entrynameToFileId,
                    std::map<std::string, size_t> & filenameToFileId,
                    std::map<size_t, std::string> & fileIdToName,
                    DBWriter* mappingWriter) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    if (par.dbExtractionMode == LocalParameters::DB_EXTRACT_MODE_INTERFACE) {
        compute3DiInterfaces(readStructure, pulchra, structureTo3Di, mat, chainNameMode, par.distanceThreshold);
    }
    size_t id = __sync_fetch_and_add(&globalCnt, readStructure.chain.size());
    size_t entriesAdded = 0;
    for (size_t ch = 0; ch < readStructure.chain.size(); ch++) {
        size_t dbKey = id + ch;
        size_t chainStart = readStructure.chain[ch].first;
        size_t chainEnd = readStructure.chain[ch].second;
        size_t chainLen = chainEnd - chainStart;
        header.clear();
        if (par.dbExtractionMode == LocalParameters::DB_EXTRACT_MODE_CHAIN) {
            if (chainLen <= 3) {
                tooShort++;
                continue;
            }
            bool allX = true;
            for (size_t pos = 0; pos < chainLen; pos++) {
                const char aa = readStructure.ami[chainStart+pos];
                if (aa != 'X' && aa != 'x') {
                    allX = false;
                    break;
                }
            }
            if (allX) {
                notProtein++;
                continue;
            }
            // Detect if structure is Ca only
            if (std::isnan(readStructure.n[chainStart + 0].x) &&
                std::isnan(readStructure.n[chainStart + 1].x) &&
                std::isnan(readStructure.n[chainStart + 2].x) &&
                std::isnan(readStructure.n[chainStart + 3].x) &&
                std::isnan(readStructure.c[chainStart + 0].x) &&
                std::isnan(readStructure.c[chainStart + 1].x) &&
                std::isnan(readStructure.c[chainStart + 2].x) &&
                std::isnan(readStructure.c[chainStart + 3].x)) {
                pulchra.rebuildBackbone(&readStructure.ca[chainStart],
                                        &readStructure.n[chainStart],
                                        &readStructure.c[chainStart],
                                        &readStructure.ami[chainStart],
                                        chainLen);

            }
            char * states = structureTo3Di.structure2states(&readStructure.ca[chainStart],
                                                            &readStructure.n[chainStart],
                                                            &readStructure.c[chainStart],
                                                            &readStructure.cb[chainStart],
                                                            chainLen);
            for (size_t pos = 0; pos < chainLen; pos++) {
                if (readStructure.ca_bfactor[pos] < maskBfactorThreshold) {
                    alphabet3di.push_back(tolower(mat.num2aa[static_cast<int>(states[pos])]));
                    alphabetAA.push_back(tolower(readStructure.ami[chainStart+pos]));
                } else {
                    alphabet3di.push_back(mat.num2aa[static_cast<int>(states[pos])]);
                    alphabetAA.push_back(readStructure.ami[chainStart+pos]);
                }
            }
            if (Util::endsWith(".gz", readStructure.names[ch]) || Util::endsWith(".zstd", readStructure.names[ch]) || Util::endsWith(".zst", readStructure.names[ch])) {
                header.append(Util::remove_extension(Util::remove_extension(readStructure.names[ch])));
            }
            else {
                header.append(Util::remove_extension(readStructure.names[ch]));
            } 
            if (readStructure.modelCount > 1 || par.modelNameMode == LocalParameters::MODEL_MODE_ADD) {
                header.append("_MODEL_");
                header.append(SSTR(readStructure.modelIndices[ch]));
            }
            if (chainNameMode == LocalParameters::CHAIN_MODE_ADD ||
                (chainNameMode == LocalParameters::CHAIN_MODE_AUTO && readStructure.names.size() > 1)) {
                header.push_back('_');
                header.append(readStructure.chainNames[ch]);
            }
            if (readStructure.title.size() > 0) {
                header.push_back(' ');
                header.append(readStructure.title);
            }
        }
        else if (par.dbExtractionMode == LocalParameters::DB_EXTRACT_MODE_INTERFACE) {
            if (readStructure.names[ch] == "ALLX") {
                notProtein++;
                continue;
            }
            if (readStructure.names[ch] == "SHORT") {
                tooShort++;
                continue;
            }
            for (size_t pos = 0; pos < chainLen; pos++) {
                alphabet3di.push_back(readStructure.seq3di[chainStart+pos]);
                alphabetAA.push_back(readStructure.ami[chainStart+pos]);
            }
            header.append(readStructure.names[ch]);
        }
        alphabet3di.push_back('\n');
        alphabetAA.push_back('\n');
        torsiondbw.writeData(alphabet3di.data(), alphabet3di.size(), dbKey, thread_idx);
        aadbw.writeData(alphabetAA.data(), alphabetAA.size(), dbKey, thread_idx);
        header.push_back('\n');
        std::string entryName = Util::parseFastaHeader(header.c_str());
#pragma omp critical
        {
            if (par.dbExtractionMode == LocalParameters::DB_EXTRACT_MODE_CHAIN) {
                std::string filenameWithExtension;
                if (Util::endsWith(".gz", readStructure.names[ch]) || Util::endsWith(".zstd", readStructure.names[ch]) || Util::endsWith(".zst", readStructure.names[ch])) {
                    filenameWithExtension = Util::remove_extension(Util::remove_extension(readStructure.names[ch]));
                } else {
                    filenameWithExtension = Util::remove_extension(readStructure.names[ch]);
                }
                std::string filenameWithoutExtension = Util::remove_extension(filenameWithExtension);
                std::map<std::string, size_t>::iterator it = filenameToFileId.find(filenameWithoutExtension);
                size_t fileid;
                if (it != filenameToFileId.end()) {
                    fileid = it->second;
                } else {
                    fileid = fileidCnt;
                    filenameToFileId[filenameWithoutExtension] = fileid;
                    fileIdToName[fileid] = filenameWithoutExtension;
                    fileidCnt++;
                }
                entrynameToFileId[entryName] = std::make_pair(fileid, readStructure.modelIndices[ch]);
            } else if (par.dbExtractionMode == LocalParameters::DB_EXTRACT_MODE_INTERFACE) {
                std::string filenameWithoutExtension;
                if (chainNameMode == LocalParameters::CHAIN_MODE_ADD || chainNameMode == LocalParameters::CHAIN_MODE_AUTO) {
                    size_t firstUnderscore = readStructure.names[ch].find_last_of('_');
                    filenameWithoutExtension = readStructure.names[ch].substr(0, firstUnderscore);
                } else {
                    filenameWithoutExtension = readStructure.names[ch];
                }
                std::map<std::string, size_t>::iterator it = filenameToFileId.find(filenameWithoutExtension);
                size_t fileid;
                if (it != filenameToFileId.end()) {
                    fileid = it->second;
                } else {
                    fileid = fileidCnt;
                    filenameToFileId[filenameWithoutExtension] = fileid;
                    fileIdToName[fileid] = filenameWithoutExtension;
                    fileidCnt++;
                }
                entrynameToFileId[entryName] = std::make_pair(fileid, readStructure.modelIndices[ch]);
            }
        }
        hdbw.writeData(header.c_str(), header.size(), dbKey, thread_idx);

        if (mappingWriter != NULL) {
            std::string taxId = SSTR(readStructure.taxIds[ch]);
            taxId.append(1, '\n');
            mappingWriter->writeData(taxId.c_str(), taxId.size(), dbKey, thread_idx, false);
        }

        float* camolf32;
        if (coordStoreMode == LocalParameters::COORD_STORE_MODE_CA_DIFF) {
            camol.resize((chainLen - 1) * 3 * sizeof(int16_t) + 3 * sizeof(float) + 1 * sizeof(uint8_t));
            int16_t* camolf16 = reinterpret_cast<int16_t*>(camol.data());
            // check if any of the coordinates is too large to be stored as int16_t
            if (Coordinate16::convertToDiff16(chainLen, (double*)(readStructure.ca.data() + chainStart) + 0, camolf16)
             || Coordinate16::convertToDiff16(chainLen, (double*)(readStructure.ca.data() + chainStart) + 1, camolf16 + 1 * (chainLen + 1))
             || Coordinate16::convertToDiff16(chainLen, (double*)(readStructure.ca.data() + chainStart) + 2, camolf16 + 2 * (chainLen + 1))) {
                // store all coordinates as float instead
                goto overflow;
            }
            cadbw.writeData((const char*)camol.data(), (chainLen - 1) * 3 * sizeof(uint16_t) + 3 * sizeof(float) + 1 * sizeof(uint8_t), dbKey, thread_idx);
            goto cleanup;
        }
overflow:
        camol.resize(chainLen * 3 * sizeof(float));
        camolf32 = reinterpret_cast<float*>(camol.data());
        for (size_t pos = 0; pos < chainLen; pos++) {
            camolf32[(0 * chainLen) + pos] = (std::isnan(readStructure.ca[chainStart+pos].x))
                        ? 0.0 : readStructure.ca[chainStart+pos].x;
        }
        for (size_t pos = 0; pos < chainLen; pos++) {
            camolf32[(1 * chainLen) + pos] = (std::isnan(readStructure.ca[chainStart+pos].y))
                        ? 0.0 : readStructure.ca[chainStart+pos].y;
        }
        for (size_t pos = 0; pos < chainLen; pos++) {
            camolf32[(2 * chainLen) + pos] = (std::isnan(readStructure.ca[chainStart+pos].z))
                        ? 0.0 : readStructure.ca[chainStart+pos].z;
        }
        cadbw.writeData((const char*)camol.data(), chainLen * 3 * sizeof(float), dbKey, thread_idx);
cleanup:
        alphabet3di.clear();
        alphabetAA.clear();
        camol.clear();
        entriesAdded++;
    }
    return entriesAdded;
}

void sortDatafileByIdOrder(DBWriter & dbw,
                           DBReader<unsigned int> &dbr,  std::vector<std::pair<std::string, unsigned int>> & order) {
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

#pragma omp for schedule(static)
        for (size_t i = 0; i < dbr.getSize(); i++) {
            size_t id = order[i].second;
            char *data = dbr.getData(id, thread_idx);
            size_t length = dbr.getEntryLen(id);
            dbw.writeData(data, (length == 0 ? 0 : length - 1), dbr.getDbKey(id), thread_idx);
        }
    }
}

extern int createdb(int argc, const char **argv, const Command& command);
int structcreatedb(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, 0, MMseqsParameter::COMMAND_COMMON);
    std::string outputName = par.filenames.back();
    if (par.prostt5Model != "") {
#ifdef HAVE_PROSTT5
        // reset set parameters
        for (size_t i = 0; i < command.params->size(); ++i) {
            command.params->at(i)->wasSet = false;
        }
        par.shuffleDatabase = true;
        par.PARAM_SHUFFLE.wasSet = true;
        int status = createdb(argc, argv, command);
        if (status != EXIT_SUCCESS) {
            return status;
        }
        fflush(stdout);

        std::vector<std::string> prefix = { "", "/model" };
        std::vector<std::string> suffix = { "", "/prostt5-f16.gguf" };
        // bool quantized = false;
        std::string modelWeights;
        for (size_t i = 0; i < prefix.size(); ++i) {
            for (size_t j = 0; j < suffix.size(); ++j) {
                std::string tensorPath = par.prostt5Model + prefix[i] + suffix[j];
                if (FileUtil::fileExists(tensorPath.c_str()) && FileUtil::directoryExists(tensorPath.c_str()) == false) {
                    modelWeights = tensorPath;
                    break;
                }
            }
        }
        if (modelWeights.empty()) {
            std::vector<std::string> prefix = { "", "/model" };
            std::vector<std::string> suffix = { "/model.safetensors" };
            std::string modelWeights;
            for (size_t i = 0; i < prefix.size(); ++i) {
                for (size_t j = 0; j < suffix.size(); ++j) {
                    std::string tensorPath = par.prostt5Model + prefix[i] + suffix[j];
                    if (FileUtil::fileExists(tensorPath.c_str())) {
                        modelWeights = par.prostt5Model + prefix[i];
                        break;
                    }
                }
            }
            if (modelWeights.empty()) {
                Debug(Debug::ERROR) << "Could not find ProstT5 model weights. Download with `foldseek databases ProstT5 prostt5_out tmp`\n";
                return EXIT_FAILURE;
            } else {
                Debug(Debug::ERROR) << "Found ProstT5 model weights for previous Foldseek release. Download new weights with `foldseek databases ProstT5 prostt5_out tmp`\n";
                return EXIT_FAILURE;
            }
        }

        LlamaInitGuard guard(par.verbosity > 3);
        std::vector<std::string> devices = ProstT5::getDevices();
        for (std::vector<std::string>::iterator it = devices.begin(); it != devices.end(); ++it) {
            Debug(Debug::INFO) << *it << "\n";
        }
        if (par.gpu == 1 && !devices.empty()) {
            for (std::vector<std::string>::iterator it = devices.begin(); it != devices.end();) {
                if (it->find("CUDA") == std::string::npos) {
                    it = devices.erase(it); // Erase returns the next iterator
                } else {
                    ++it; // Move to the next element
                }
            }
            if (devices.size() == 0) {
                Debug(Debug::ERROR) << "No GPU devices found\n";
                return EXIT_FAILURE;
            }
        } else {
            for (size_t i = 0; i < devices.size(); i++) {
                if (devices[i] == "Metal") {
                    par.gpu = 1;
                    devices.clear();
                    devices.push_back("Metal");
                    break;
                }
            }
        }

        bool useForkRunner = FORK_RUNNER && par.gpu == 0;
        DBReader<unsigned int> reader(outputName.c_str(), (outputName+".index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        reader.open(useForkRunner ? DBReader<unsigned int>::SORT_BY_LENGTH : DBReader<unsigned int>::LINEAR_ACCCESS);

        unsigned const int MIN_SPLIT_LENGTH = 2;
        std::string ssDb = outputName + "_ss";
        std::string ssIndex = ssDb + ".index"; 
        if (useForkRunner) {
            prostt5Forking(modelWeights, par.prostt5SplitLength, MIN_SPLIT_LENGTH, reader, ssDb, ssIndex, par.threads, par.compressed);
        } else {
            DBWriter writer(ssDb.c_str(), ssIndex.c_str(), par.threads, par.compressed, reader.getDbtype());
            writer.open();

            Debug::Progress progress(reader.getSize());
#ifdef OPENMP
            size_t localThreads = par.gpu == 1 ? devices.size() : 1;
#endif
#pragma omp parallel num_threads(localThreads)
            {
                int thread_idx = 0;
#ifdef OPENMP
                thread_idx = omp_get_thread_num();
#endif
                std::string device = "none";
                int localThreads = par.threads;
                if (par.gpu == 1) {
                    device = devices[thread_idx];
                    localThreads = 1;
                }
                ProstT5Model model(modelWeights.c_str(), device);
                ProstT5 context(model, localThreads);
                const char newline = '\n';
                std::string result;
#pragma omp for schedule(dynamic, 1)
                for (size_t i = 0; i < reader.getSize(); ++i) {
                    unsigned int key = reader.getDbKey(i);
                    size_t length = reader.getSeqLen(i);
                    std::string seq = std::string(reader.getData(i, thread_idx), length);
                    result.clear();

                    // splitting input sequences longer than ProstT5 attention (current cutoff 6000 AAs)
                    unsigned int split_length = par.prostt5SplitLength;
                    // split lenght of 0 will deactivate splitting
                    if (split_length > 0 && length > split_length) {
                        unsigned int n_splits, overlap_length;
                        n_splits = int(length / split_length) + 1;
                        overlap_length = length % split_length;

                        // ensure minimum overlap length; adjustment length was not computed properly with ceil/ceilf now using simple int cast
                        if (overlap_length < MIN_SPLIT_LENGTH) {
                            split_length -= int((MIN_SPLIT_LENGTH - overlap_length) / (n_splits - 1)) + 1;
                        }

                        // loop over splits and predict
                        for (unsigned int i = 0; i < n_splits; i++){
                            unsigned int split_start = i * split_length;
                            result.append(context.predict(seq.substr(split_start, split_length)));
                        }
                    } else {
                        result.append(context.predict(seq));
                    }

                    writer.writeStart(thread_idx);
                    writer.writeAdd(result.c_str(), result.length(), thread_idx);
                    writer.writeAdd(&newline, 1, thread_idx);
                    writer.writeEnd(key, thread_idx);
                    progress.updateProgress();
                }
            }
            writer.close(true);
        }
        reader.close();

        DBReader<unsigned int> resultReader(ssDb.c_str(), (ssDb+".index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        resultReader.open(DBReader<unsigned int>::NOSORT);
        resultReader.readMmapedDataInMemory();
        const std::pair<std::string, std::string> tempDb = Util::databaseNames(ssDb + "_tmp");
        DBWriter resultWriter(tempDb.first.c_str(), tempDb.second.c_str(), par.threads, par.compressed, resultReader.getDbtype());
        resultWriter.open();
        resultWriter.sortDatafileByIdOrder(resultReader);
        resultWriter.close(true);
        resultReader.close();
        DBReader<unsigned int>::removeDb(ssDb);
        DBReader<unsigned int>::moveDb(tempDb.first, ssDb);

        return EXIT_SUCCESS;
#else
        Debug(Debug::ERROR) << "Foldseek was compiled without ProstT5 support\n";
        return EXIT_FAILURE;
#endif
    } else {
        par.printParameters(command.cmd, argc, argv, *command.params);
    }
    par.filenames.pop_back();

    PatternCompiler include(par.fileInclude.c_str());
    PatternCompiler exclude(par.fileExclude.c_str());

    for (size_t i = 1; i < par.filenames.size(); ++i) {
        if (FileUtil::directoryExists(par.filenames[i].c_str()) || Util::endsWith(".tsv", par.filenames[i].c_str())) {
            Debug(Debug::ERROR) << "Only one directory or tsv file (" << par.filenames[i] << ") or a list of files can be given\n";
            EXIT(EXIT_FAILURE);
        }
    }

    if (FileUtil::directoryExists(par.filenames[0].c_str())) {
        if (par.filenames.size() > 1) {
            Debug(Debug::ERROR) << "Only one directory can be given\n";
            EXIT(EXIT_FAILURE);
        }
        std::vector<std::string> dirs;
        dirs.push_back(par.filenames.back());
        par.filenames.pop_back();
        while (dirs.size() != 0) {
            std::string dir = dirs.back();
            dirs.pop_back();
            DIR* handle = opendir(dir.c_str());
            if (handle == NULL) {
                continue;
            }
            while (dirent* entry = readdir(handle)) {
                std::string filename(entry->d_name);
                if (filename != "." && filename != "..") {
                    std::string fullpath = dir + "/" + filename;
                    struct stat info;
                    stat(fullpath.c_str(), &info);
                    if (info.st_mode & S_IFDIR) {
                        dirs.push_back(fullpath);
                    } else if (include.isMatch(filename.c_str()) == true && exclude.isMatch(filename.c_str()) == false) {
                        par.filenames.push_back(fullpath);
                    }
                }
            }
            closedir(handle);
        }
    } else if (Util::endsWith(".tsv", par.filenames[0])) {
        if (par.filenames.size() > 1) {
            Debug(Debug::ERROR) << "Only one tsv file can be given\n";
            EXIT(EXIT_FAILURE);
        }
        std::string tsv = par.filenames.back();
        par.filenames.pop_back();

        FILE* file = FileUtil::openFileOrDie(tsv.c_str(), "r", true);
        char* line = NULL;
        size_t len = 0;
        ssize_t read;
        while ((read = getline(&line, &len, file)) != -1) {
            if (line[read - 1] == '\n') {
                line[read - 1] = '\0';
                read--;
            }
            par.filenames.push_back(line);
        }
        free(line);
        fclose(file);
    }

    Debug(Debug::INFO) << "Output file: " << outputName << "\n";
    SORT_PARALLEL(par.filenames.begin(), par.filenames.end());

    DBWriter torsiondbw((outputName+"_ss").c_str(), (outputName+"_ss.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    torsiondbw.open();
    DBWriter hdbw((outputName+"_h").c_str(), (outputName+"_h.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_GENERIC_DB);
    hdbw.open();
    DBWriter cadbw((outputName+"_ca").c_str(), (outputName+"_ca.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, LocalParameters::DBTYPE_CA_ALPHA);
    cadbw.open();
    DBWriter aadbw((outputName).c_str(), (outputName+".index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    aadbw.open();
    DBWriter* mappingWriter = NULL;
    if (par.writeMapping) {
        mappingWriter = new DBWriter((outputName+"_mapping_tmp").c_str(), (outputName+"_mapping_tmp.index").c_str(), static_cast<unsigned int>(par.threads), false, LocalParameters::DBTYPE_OMIT_FILE);
        mappingWriter->open();
    }

    SubstitutionMatrix mat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);
    Debug::Progress progress(par.filenames.size());
    std::map<std::string, std::pair<size_t, unsigned int>> entrynameToFileId;
    std::map<size_t, std::string> fileIdToName;
    std::map<std::string, size_t> filenameToFileId;

    std::vector<std::string> tarFiles;
    std::vector<std::string> looseFiles;
    std::vector<std::string> gcsPaths;
    std::vector<std::string> dbs;

    for (size_t i = 0; i < par.filenames.size(); ++i) {
        if (Util::endsWith(".tar.gz", par.filenames[i])
                || Util::endsWith(".tgz", par.filenames[i])
                || Util::endsWith(".tar.zst", par.filenames[i])
                || Util::endsWith(".tar.zstd", par.filenames[i])
                || Util::endsWith(".tar", par.filenames[i])) {
            tarFiles.push_back(par.filenames[i]);
        }
#ifdef HAVE_GCS
        else if (Util::startWith("gcs://", par.filenames[i])) {
            gcsPaths.push_back(par.filenames[i]);
        }
#endif
        else if (FileUtil::fileExists((par.filenames[i] + ".dbtype").c_str())) {
            dbs.push_back(par.filenames[i]);
        } else {
            looseFiles.push_back(par.filenames[i]);
        }
    }

    size_t globalCnt = 0;
    size_t globalFileidCnt = 0;
    size_t incorrectFiles = 0;
    size_t tooShort = 0;
    size_t notProtein = 0;
    bool needsReorderingAtTheEnd = false;
    size_t needToWriteModel = 0;

    // cannot be const for compatibility with older compilers/openmp and omp shared
    int inputFormat = par.inputFormat;

    // Process tar files!
    for (size_t i = 0; i < tarFiles.size(); i++) {
        mtar_t tar;
        if (Util::endsWith(".tar.gz", tarFiles[i]) || Util::endsWith(".tgz", tarFiles[i])) {
#ifdef HAVE_ZLIB
            if (structure_mtar_gzopen(&tar, tarFiles[i].c_str()) != MTAR_ESUCCESS) {
                Debug(Debug::ERROR) << "Cannot open file " << tarFiles[i] << "\n";
                EXIT(EXIT_FAILURE);
            }
#else
            Debug(Debug::ERROR) << "Foldseek was not compiled with zlib support. Cannot read compressed input.\n";
            EXIT(EXIT_FAILURE);
#endif
        } else if (Util::endsWith(".tar.zstd", tarFiles[i]) || Util::endsWith(".tar.zst", tarFiles[i])) {
            if (mtar_zstdopen(&tar, tarFiles[i].c_str()) != MTAR_ESUCCESS) {
                Debug(Debug::ERROR) << "Cannot open file " << tarFiles[i] << "\n";
                EXIT(EXIT_FAILURE);
            }
        } else {
            if (mtar_open(&tar, tarFiles[i].c_str(), "r") != MTAR_ESUCCESS) {
                Debug(Debug::ERROR) << "Cannot open file " << tarFiles[i] << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        progress.updateProgress();
#ifdef OPENMP
        int localThreads = par.threads;
        if (localThreads > 1) {
            needsReorderingAtTheEnd = true;
        }
#endif

#pragma omp parallel default(none) shared(tar, par, torsiondbw, hdbw, cadbw, aadbw, mat, progress, globalCnt, globalFileidCnt, entrynameToFileId, filenameToFileId, fileIdToName, mappingWriter, std::cerr, std::cout, inputFormat) num_threads(localThreads) reduction(+:incorrectFiles, tooShort, notProtein, needToWriteModel)
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#ifdef HAVE_ZLIB
            const unsigned int CHUNK = 128 * 1024;
            unsigned char in[CHUNK];
            unsigned char out[CHUNK];
            z_stream strm;
            memset(&strm, 0, sizeof(z_stream));
            strm.zalloc = Z_NULL;
            strm.zfree = Z_NULL;
            strm.opaque = Z_NULL;
            strm.next_in = in;
            strm.avail_in = 0;
            int status = inflateInit2(&strm, 15 | 32);
            if (status < 0) {
                Debug(Debug::ERROR) << "Cannot initialize zlib stream\n";
                EXIT(EXIT_FAILURE);
            }
#endif
            //recon_related
            StructureTo3Di structureTo3Di;
            PulchraWrapper pulchra;
            GemmiWrapper readStructure;
            std::vector<char> alphabet3di;
            std::vector<char> alphabetAA;
            std::vector<int8_t> camol;
            std::string header;
            std::string name;
            std::string pdbFile;
            mtar_header_t tarHeader;
            size_t bufferSize = 1024 * 1024;
            char *dataBuffer = (char *) malloc(bufferSize);
            size_t inflateSize = 1024 * 1024;
            char *inflateBuffer = (char *) malloc(inflateSize);
            bool proceed = true;
            PatternCompiler includeThread(par.fileInclude.c_str());
            PatternCompiler excludeThread(par.fileExclude.c_str());

            while (proceed) {
                bool writeEntry = true;
#pragma omp critical
                {
                    if (tar.isFinished == 0 && (mtar_read_header(&tar, &tarHeader)) != MTAR_ENULLRECORD) {
                        // GNU tar has special blocks for long filenames
                        if (tarHeader.type == MTAR_TGNU_LONGNAME || tarHeader.type == MTAR_TGNU_LONGLINK) {
                            if (tarHeader.size == 0) {
                                Debug(Debug::ERROR) << "Invalid tar input with long name/link record of size 0\n";
                                EXIT(EXIT_FAILURE);
                            }
                            if (tarHeader.size > bufferSize) {
                                bufferSize = tarHeader.size * 1.5;
                                dataBuffer = (char *) realloc(dataBuffer, bufferSize);
                            }
                            if (mtar_read_data(&tar, dataBuffer, tarHeader.size) != MTAR_ESUCCESS) {
                                Debug(Debug::ERROR) << "Cannot read entry " << tarHeader.name << "\n";
                                EXIT(EXIT_FAILURE);
                            }
                            // skip null byte
                            name.assign(dataBuffer, tarHeader.size - 1);
                            // skip to next record
                            if (mtar_read_header(&tar, &tarHeader) == MTAR_ENULLRECORD) {
                                Debug(Debug::ERROR) << "Tar truncated after entry " << name << "\n";
                                EXIT(EXIT_FAILURE);
                            }
                        } else {
                            name = tarHeader.name;
                        }
                        if (tarHeader.type == MTAR_TREG || tarHeader.type == MTAR_TCONT ||
                            tarHeader.type == MTAR_TOLDREG) {
                            if (tarHeader.size > bufferSize) {
                                bufferSize = tarHeader.size * 1.5;
                                dataBuffer = (char *) realloc(dataBuffer, bufferSize);
                            }
                            if (mtar_read_data(&tar, dataBuffer, tarHeader.size) != MTAR_ESUCCESS) {
                                Debug(Debug::ERROR) << "Cannot read entry " << name << "\n";
                                EXIT(EXIT_FAILURE);
                            }
                            if (includeThread.isMatch(name.c_str()) == false || excludeThread.isMatch(name.c_str()) == true) {
                                proceed = true;
                                writeEntry = false;
                            } else {
                                proceed = true;
                                writeEntry = true;
                            }
                        } else {
                            proceed = true;
                            writeEntry = false;
                        }
                    } else {
                        tar.isFinished = 1;
                        proceed = false;
                        writeEntry = false;
                    }
                }

                if (writeEntry) {
                    pdbFile.clear();
                    if (Util::endsWith(".gz", name)) {
#ifdef HAVE_ZLIB
                        inflateReset(&strm);
                        strm.avail_in = tarHeader.size;
                        strm.next_in = (unsigned char *) dataBuffer;
                        do {
                            unsigned have;
                            strm.avail_out = CHUNK;
                            strm.next_out = out;
                            int err = inflate(&strm, Z_NO_FLUSH);
                            switch (err) {
                                case Z_OK:
                                case Z_STREAM_END:
                                case Z_BUF_ERROR:
                                    break;
                                default:
                                    inflateEnd(&strm);
                                    Debug(Debug::ERROR) << "Gzip error " << err << " entry " << name << "\n";
                                    EXIT(EXIT_FAILURE);
                            }
                            have = CHUNK - strm.avail_out;
                            pdbFile.append((char *) out, have);
                        } while (strm.avail_out == 0);
#else
                        Debug(Debug::ERROR) << "MMseqs2 was not compiled with zlib support. Cannot read compressed input.\n";
                        EXIT(EXIT_FAILURE);
#endif
                    } else {
                        pdbFile.append(dataBuffer, tarHeader.size);
                    }
                    if (readStructure.loadFromBuffer(pdbFile.c_str(), pdbFile.size(), name, (GemmiWrapper::Format)inputFormat) == false) {
                        incorrectFiles++;
                        continue;
                    }

                    __sync_add_and_fetch(&needToWriteModel, (readStructure.modelCount > 1));
                    writeStructureEntry(
                        mat, readStructure, structureTo3Di, pulchra,
                        alphabet3di, alphabetAA, camol, header, aadbw, hdbw, torsiondbw, cadbw,
                        par.chainNameMode, par.maskBfactorThreshold, tooShort, notProtein, globalCnt, thread_idx, par.coordStoreMode,
                        globalFileidCnt, entrynameToFileId, filenameToFileId, fileIdToName,
                        mappingWriter
                    );
                }
            } // end while
            free(inflateBuffer);
            free(dataBuffer);
        } // end omp open
        mtar_close(&tar);
    } // end file for


    //===================== single_process ===================//__110710__//
#pragma omp parallel default(none) shared(par, torsiondbw, hdbw, cadbw, aadbw, mat, looseFiles, progress, globalCnt, globalFileidCnt, entrynameToFileId, filenameToFileId, fileIdToName, mappingWriter, inputFormat) reduction(+:incorrectFiles, tooShort, notProtein, needToWriteModel)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        //recon_related
        StructureTo3Di structureTo3Di;
        PulchraWrapper pulchra;
        GemmiWrapper readStructure;
        std::vector<char> alphabet3di;
        std::vector<char> alphabetAA;
        std::vector<int8_t> camol;
        std::string header;
        std::string name;

#pragma omp for schedule(static)
        for (size_t i = 0; i < looseFiles.size(); i++) {
            progress.updateProgress();

            if (readStructure.load(looseFiles[i], (GemmiWrapper::Format)inputFormat) == false) {
                incorrectFiles++;
                continue;
            }
            __sync_add_and_fetch(&needToWriteModel, (readStructure.modelCount > 1));
            // clear memory
            writeStructureEntry(
                mat, readStructure, structureTo3Di,  pulchra,
                alphabet3di, alphabetAA, camol, header, aadbw, hdbw, torsiondbw, cadbw,
                par.chainNameMode, par.maskBfactorThreshold, tooShort, notProtein, globalCnt, thread_idx, par.coordStoreMode,
                globalFileidCnt, entrynameToFileId, filenameToFileId, fileIdToName,
                mappingWriter
            );
        }
    }

#ifdef HAVE_GCS
    namespace gcs = ::google::cloud::storage;
    auto options = google::cloud::Options{}
        .set<gcs::ConnectionPoolSizeOption>(par.threads)
        .set<google::cloud::storage_experimental::HttpVersionOption>("2.0");
    auto client = gcs::Client(options);
    for (size_t i = 0; i < gcsPaths.size(); i++) {
        std::vector<std::string> parts = Util::split(gcsPaths[i], "/");
        if (parts.size() == 1) {
            Debug(Debug::WARNING) << "Skipping invalid URI " << gcsPaths[i] << "\n";
            continue;
        }
        std::string bucket_name = parts[1];

        char filter = '\0';
        if (parts.size >= 3) {
            filter = parts[2][0];
        }
        progress.reset(SIZE_MAX);
#pragma omp parallel default(none) shared(par, torsiondbw, hdbw, cadbw, aadbw, mat, gcsPaths, progress, globalCnt, globalFileidCnt, entrynameToFileId, filenameToFileId, fileIdToName, client, bucket_name, filter, mappingWriter) reduction(+:incorrectFiles, tooShort, notProtein, needToWriteModel, inputFormat)
        {
            StructureTo3Di structureTo3Di;
            PulchraWrapper pulchra;
            GemmiWrapper readStructure;
            std::vector<char> alphabet3di;
            std::vector<char> alphabetAA;
            std::vector<int8_t> camol;
            std::string header;
            std::string name;

#pragma omp single
            for (auto&& object_metadata : client.ListObjects(bucket_name, gcs::Projection::NoAcl(), gcs::MaxResults(15000))) {
                std::string obj_name = object_metadata->name();
#pragma omp task firstprivate(obj_name, alphabet3di, alphabetAA, camol, header, name, filter) private(structureTo3Di, pulchra, readStructure)
                {
                    bool skipFilter = filter != '\0' && obj_name.length() >= 9 && obj_name[8] == filter;
                    bool allowedSuffix = Util::endsWith(".cif", obj_name) || Util::endsWith(".pdb", obj_name);
                    if (skipFilter && allowedSuffix) {
                        unsigned int thread_idx = 0;
#ifdef OPENMP
                        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
                        progress.updateProgress();

                        auto reader = client.ReadObject(bucket_name, obj_name);
                        if (!reader.status().ok()) {
                            Debug(Debug::ERROR) << reader.status().message() << "\n";
                        } else {
                            std::string contents{std::istreambuf_iterator<char>{reader}, {}};
                            if (readStructure.loadFromBuffer(contents.c_str(), contents.size(), obj_name, inputFormat) == false) {
                                incorrectFiles++;
                            } else {
                                __sync_add_and_fetch(&needToWriteModel, (readStructure.modelCount > 1));
                                writeStructureEntry(
                                    mat, readStructure, structureTo3Di,  pulchra,
                                    alphabet3di, alphabetAA, camol, header, aadbw, hdbw, torsiondbw, cadbw,
                                    par.chainNameMode, par.maskBfactorThreshold, tooShort, notProtein, globalCnt, thread_idx, par.coordStoreMode,
                                    globalFileidCnt, entrynameToFileId, filenameToFileId, fileIdToName,
                                    mappingWriter
                                );
                            }
                        }
                    }
                }
            }
        }
    }
#endif

    for (size_t i = 0; i < dbs.size(); ++i) {
        DBReader<unsigned int> reader(dbs[i].c_str(), (dbs[i]+".index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_LOOKUP);
        reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
        progress.reset(reader.getSize());
#pragma omp parallel default(none) shared(par, torsiondbw, hdbw, cadbw, aadbw, mat, progress, globalCnt, globalFileidCnt, entrynameToFileId, filenameToFileId, fileIdToName, reader, mappingWriter, inputFormat) reduction(+:incorrectFiles, tooShort, notProtein, needToWriteModel)
        {
            StructureTo3Di structureTo3Di;
            PulchraWrapper pulchra;
            GemmiWrapper readStructure;
            std::vector<char> alphabet3di;
            std::vector<char> alphabetAA;
            std::vector<int8_t> camol;
            std::string header;
            std::string name;

            std::string dbname = reader.getDataFileName();

            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(static)
            for (size_t i = 0; i < reader.getSize(); i++) {
                progress.updateProgress();

                char* data = reader.getData(i, thread_idx);
                size_t len = reader.getEntryLen(i);

                size_t lookupId = reader.getLookupIdByKey(reader.getDbKey(i));
                std::string name = reader.getLookupEntryName(lookupId);

                if (readStructure.loadFromBuffer(data, len, name, (GemmiWrapper::Format)inputFormat) == false) {
                    incorrectFiles++;
                } else {
                    __sync_add_and_fetch(&needToWriteModel, (readStructure.modelCount > 1));
                    writeStructureEntry(
                        mat, readStructure, structureTo3Di,  pulchra,
                        alphabet3di, alphabetAA, camol, header, aadbw, hdbw, torsiondbw, cadbw,
                        par.chainNameMode, par.maskBfactorThreshold, tooShort, notProtein, globalCnt, thread_idx, par.coordStoreMode,
                        globalFileidCnt, entrynameToFileId, filenameToFileId, fileIdToName,
                        mappingWriter
                    );
                }
            }
        }
        reader.close();
    }

    torsiondbw.close(true);
    hdbw.close(true);
    cadbw.close(true);
    aadbw.close(true);
    if (par.writeMapping) {
        mappingWriter->close(true);
        delete mappingWriter;
    }

    if (needsReorderingAtTheEnd) {
        Debug(Debug::INFO) << "Reordering by identifier\n";
        DBReader<unsigned int> header_reorder((outputName+"_h").c_str(), (outputName+"_h.index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        header_reorder.open(DBReader<unsigned int>::NOSORT);
        header_reorder.readMmapedDataInMemory();
        std::vector<std::pair<std::string, unsigned int>> mappingOrder(header_reorder.getSize());
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(static)
            for (size_t i = 0; i < header_reorder.getSize(); i++) {
                std::string entryNameRaw = Util::parseFastaHeader(header_reorder.getData(i,thread_idx));
                std::string filenameWithoutExtension;
                if (Util::endsWith(".gz", entryNameRaw) || Util::endsWith(".zstd", entryNameRaw) || Util::endsWith(".zstd", entryNameRaw)) {
                    filenameWithoutExtension = Util::remove_extension(Util::remove_extension(Util::remove_extension(entryNameRaw)));
                } else {
                    filenameWithoutExtension = Util::remove_extension(Util::remove_extension(entryNameRaw));
                }
                mappingOrder[i].first = filenameWithoutExtension;
                mappingOrder[i].second = i;
            }
        }
        SORT_PARALLEL(mappingOrder.begin(), mappingOrder.end());
        std::string lookupFile = outputName + ".lookup";
        FILE* file = FileUtil::openAndDelete(lookupFile.c_str(), "w");
        std::string buffer;
        buffer.reserve(2048);
        for (unsigned int id = 0; id < header_reorder.getSize(); id++) {
            DBReader<unsigned int>::LookupEntry entry;
            entry.id = id;
            entry.entryName = mappingOrder[id].first;
            entry.fileNumber = 0;
            header_reorder.lookupEntryToBuffer(buffer, entry);
            size_t written = fwrite(buffer.c_str(), sizeof(char), buffer.size(), file);
            if (written != buffer.size()) {
                Debug(Debug::ERROR) << "Cannot write to lookup file " << lookupFile << "\n";
                EXIT(EXIT_FAILURE);
            }
            buffer.clear();
        }
        if (fclose(file) != 0) {
            Debug(Debug::ERROR) << "Cannot close lookup file " << lookupFile << "\n";
            EXIT(EXIT_FAILURE);
        }

        DBWriter hdbw_reorder((outputName+"_h").c_str(), (outputName+"_h.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_GENERIC_DB);
        hdbw_reorder.open();
        sortDatafileByIdOrder(hdbw_reorder, header_reorder, mappingOrder);
        hdbw_reorder.close(true);
        header_reorder.close();

        DBReader<unsigned int> torsiondbr_reorder((outputName+"_ss").c_str(), (outputName+"_ss.index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        torsiondbr_reorder.open(DBReader<unsigned int>::NOSORT);
        torsiondbr_reorder.readMmapedDataInMemory();
        DBWriter torsiondbw_reorder((outputName+"_ss").c_str(), (outputName+"_ss.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
        torsiondbw_reorder.open();
        sortDatafileByIdOrder(torsiondbw_reorder, torsiondbr_reorder, mappingOrder);
        torsiondbw_reorder.close(true);
        torsiondbr_reorder.close();

        DBReader<unsigned int> cadbr_reorder((outputName+"_ca").c_str(), (outputName+"_ca.index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        cadbr_reorder.open(DBReader<unsigned int>::NOSORT);
        cadbr_reorder.readMmapedDataInMemory();
        DBWriter cadbw_reorder((outputName+"_ca").c_str(), (outputName+"_ca.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, LocalParameters::DBTYPE_CA_ALPHA);
        cadbw_reorder.open();
        sortDatafileByIdOrder(cadbw_reorder, cadbr_reorder, mappingOrder);
        cadbw_reorder.close(true);
        cadbr_reorder.close();

        DBReader<unsigned int> aadbr_reorder((outputName).c_str(), (outputName+".index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        aadbr_reorder.open(DBReader<unsigned int>::NOSORT);
        aadbr_reorder.readMmapedDataInMemory();
        DBWriter aadbw_reorder((outputName).c_str(), (outputName+".index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
        aadbw_reorder.open();
        sortDatafileByIdOrder(aadbw_reorder, aadbr_reorder, mappingOrder);
        aadbw_reorder.close(true);
        aadbr_reorder.close();

        if (par.writeMapping) {
            DBReader<unsigned int> mappingReader_reorder((outputName+"_mapping_tmp").c_str(), (outputName+"_mapping_tmp.index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
            mappingReader_reorder.open(DBReader<unsigned int>::NOSORT);
            mappingReader_reorder.readMmapedDataInMemory();
            DBWriter mappingWriter_reorder((outputName+"_mapping_tmp").c_str(), (outputName+"_mapping_tmp.index").c_str(), static_cast<unsigned int>(par.threads), par.compressed, LocalParameters::DBTYPE_CA_ALPHA);
            mappingWriter_reorder.open();
            sortDatafileByIdOrder(mappingWriter_reorder, mappingReader_reorder, mappingOrder);
            mappingWriter_reorder.close(true);
            mappingReader_reorder.close();
        }
    } else {
        DBWriter::createRenumberedDB((outputName+"_ss").c_str(), (outputName+"_ss.index").c_str(), "", "", DBReader<unsigned int>::LINEAR_ACCCESS);
        DBWriter::createRenumberedDB((outputName+"_h").c_str(), (outputName+"_h.index").c_str(), "", "", DBReader<unsigned int>::LINEAR_ACCCESS);
        DBWriter::createRenumberedDB((outputName+"_ca").c_str(), (outputName+"_ca.index").c_str(), "", "", DBReader<unsigned int>::LINEAR_ACCCESS);
        DBWriter::createRenumberedDB((outputName).c_str(), (outputName+".index").c_str(), "", "", DBReader<unsigned int>::LINEAR_ACCCESS);
        if (par.writeMapping) {
            DBWriter::createRenumberedDB((outputName+"_mapping_tmp").c_str(), (outputName+"_mapping_tmp.index").c_str(), "", "", DBReader<unsigned int>::LINEAR_ACCCESS);
        }
    }

    if (par.writeMapping) {
        std::string mappingFile = outputName + "_mapping";
        FILE* file = FileUtil::openAndDelete(mappingFile.c_str(), "w");
        char buffer[1024];
        DBReader<unsigned int> mappingReader((outputName+"_mapping_tmp").c_str(), (outputName+"_mapping_tmp.index").c_str(), 1, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        mappingReader.open(DBReader<unsigned int>::NOSORT);
        for (size_t i = 0; i < mappingReader.getSize(); ++i) {
            unsigned int key = mappingReader.getDbKey(i);
            char* data = mappingReader.getData(i, 0);
            size_t entryLength = mappingReader.getEntryLen(i) - 1;
            char* basePos = buffer;
            char* tmpBuff = Itoa::u32toa_sse2(static_cast<uint32_t>(key), buffer);
            *(tmpBuff-1) = '\t';
            strncpy(tmpBuff, data, entryLength);
            tmpBuff = tmpBuff + entryLength;
            *(tmpBuff) = '\n';
            size_t length = tmpBuff - basePos + 1;
            size_t written = fwrite(buffer, sizeof(char), length, file);
            if (written != length) {
                Debug(Debug::ERROR) << "Error writing to file " << mappingFile << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        mappingReader.close();
        if (fclose(file) != 0) {
            Debug(Debug::ERROR) << "Cannot close file " << mappingFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        DBReader<unsigned int>::removeDb(outputName + "_mapping_tmp");
    }

    if (par.writeLookup == true) {
        DBReader<unsigned int> readerHeader((outputName + "_h").c_str(), (outputName + "_h.index").c_str(),
                                            1, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        readerHeader.open(DBReader<unsigned int>::NOSORT);
        // create lookup file
        std::string lookupFile = outputName + ".lookup";
        FILE *file = FileUtil::openAndDelete(lookupFile.c_str(), "w");
        std::string buffer;
        size_t maxFileId = 0;
        buffer.reserve(2048);
        DBReader<unsigned int>::LookupEntry entry;
        std::string sourceFilename = outputName + ".source";
        FILE *sourceFile = FileUtil::openAndDelete(sourceFilename.c_str(), "w");
        std::map< std::pair<size_t, unsigned int>, size_t> modelFileIdLookup;
        size_t globalFileNumber = 0;
        bool needToWriteSource = false;
        for (unsigned int id = 0; id < readerHeader.getSize(); id++) {
            char *header = readerHeader.getData(id, 0);
            entry.id = readerHeader.getDbKey(id);
            std::string entryNameRaw = Util::parseFastaHeader(header);
            std::pair<size_t, unsigned int> fileIdModelEntry = entrynameToFileId[entryNameRaw];

            if(par.modelNameMode == LocalParameters::MODEL_MODE_ADD) {
                size_t modelIndex = entryNameRaw.find("MODEL");
                if (modelIndex == std::string::npos) {
                    entryNameRaw.append("_MODEL_");
                    entryNameRaw.append(SSTR(fileIdModelEntry.second));
                }
            }
            entry.entryName = entryNameRaw;
            size_t fileId = fileIdModelEntry.first;
            if (modelFileIdLookup.find(fileIdModelEntry) == modelFileIdLookup.end()) {
                modelFileIdLookup[fileIdModelEntry] = globalFileNumber;
                entry.fileNumber = globalFileNumber;
                globalFileNumber++;
                needToWriteSource = true;
            } else {
                entry.fileNumber = modelFileIdLookup[fileIdModelEntry];
                needToWriteSource = false;
            }
            maxFileId = std::max(fileId, maxFileId);
            readerHeader.lookupEntryToBuffer(buffer, entry);
            size_t written = fwrite(buffer.c_str(), sizeof(char), buffer.size(), file);
            if (written != buffer.size()) {
                Debug(Debug::ERROR) << "Cannot write to lookup file " << lookupFile << "\n";
                EXIT(EXIT_FAILURE);
            }
            buffer.clear();
            if (needToWriteSource) {
                std::string filename = FileUtil::baseName(fileIdToName[fileId]);
                if (needToWriteModel) {
                    filename.append("_MODEL_");
                    filename.append(SSTR(fileIdModelEntry.second));
                }
                buffer.append(SSTR(entry.fileNumber));
                buffer.push_back('\t');
                buffer.append(filename);
                buffer.push_back('\n');
                written = fwrite(buffer.c_str(), sizeof(char), buffer.size(), sourceFile);
                if (written != buffer.size()) {
                    Debug(Debug::ERROR) << "Cannot write to lookup file " << lookupFile << "\n";
                    EXIT(EXIT_FAILURE);
                }
                buffer.clear();
            }
        }
        if (fclose(file) != 0) {
            Debug(Debug::ERROR) << "Cannot close lookup file " << lookupFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (fclose(sourceFile) != 0) {
            Debug(Debug::ERROR) << "Cannot close source file " << sourceFilename << "\n";
            EXIT(EXIT_FAILURE);
        }
        readerHeader.close();
    }
    if (globalCnt == 0) {
        Debug(Debug::ERROR) << "No structures found in given input.\n";
        EXIT(EXIT_FAILURE);
    }
    
    Debug(Debug::INFO) << "Ignore " << (tooShort+incorrectFiles+notProtein) << " out of " << globalCnt << ".\n";
    Debug(Debug::INFO) << "Too short: " << tooShort << ", incorrect: " << incorrectFiles << ", not proteins: " << notProtein << ".\n";
    return EXIT_SUCCESS;
}
