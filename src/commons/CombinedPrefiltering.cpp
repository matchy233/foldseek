//
// Created by minghang on 10/7/24.
//

#include "CombinedPrefiltering.h"
#include "DBWriter.h"
#include "QueryMatcher.h"
#include "QueryMatcherTaxonomyHook.h"

bool CombinedPrefiltering::runSplit(const std::string &resultDB, const std::string &resultDBIndex, size_t split,
                                    bool merge) {
    Debug(Debug::INFO) << "Process prefiltering step " << (split + 1) << " of " << splits << "\n\n";

    size_t dbFrom = 0;
    size_t dbSize = tdbr->getSize();
    size_t queryFrom = 0;
    size_t querySize = qdbr->getSize();

    // create index table based on split parameter
    if (splitMode == Parameters::TARGET_DB_SPLIT) {
        tdbr->decomposeDomainByAminoAcid(split, splits, &dbFrom, &dbSize);
        if (dbSize == 0) {
            return false;
        }

        if (indexTable != nullptr) {
            delete indexTable;
            indexTable = nullptr;
        }

        if (sequenceLookup != nullptr) {
            delete sequenceLookup;
            sequenceLookup = nullptr;
        }

        getIndexTable(split, dbFrom, dbSize);
    } else if (splitMode == Parameters::QUERY_DB_SPLIT) {
        qdbr->decomposeDomainByAminoAcid(split, splits, &queryFrom, &querySize);
        if (querySize == 0) {
            return false;
        }
    }

    Debug(Debug::INFO) << "k-mer similarity threshold: " << kmerThr << "\n";

    double kmersPerPos = 0;
    size_t dbMatches = 0;
    size_t doubleMatches = 0;
    size_t querySeqLenSum = 0;
    size_t resSize = 0;
    size_t realResSize = 0;
    size_t diagonalOverflow = 0;
    size_t totalQueryDBSize = querySize;

    size_t localThreads = 1;
#ifdef OPENMP
    localThreads = std::max(std::min((size_t) threads, querySize), (size_t) 1);
#endif

    DBWriter tmpDbw(resultDB.c_str(), resultDBIndex.c_str(), localThreads, compressed,
                    Parameters::DBTYPE_PREFILTER_RES);
    tmpDbw.open();

    // init all thread-specific data structures
    char *notEmpty = new char[querySize];
    memset(notEmpty, 0, querySize * sizeof(char)); // init notEmpty

    auto **reslens = new std::list<int> *[localThreads];
    for (size_t i = 0; i < localThreads; ++i) {
        reslens[i] = new std::list<int>();
    }

    Debug(Debug::INFO) << "Starting prefiltering scores calculation (step " << (split + 1) << " of " << splits << ")\n";
    Debug(Debug::INFO) << "Query db start " << (queryFrom + 1) << " to " << queryFrom + querySize << "\n";
    Debug(Debug::INFO) << "Target db start " << (dbFrom + 1) << " to " << dbFrom + dbSize << "\n";
    Debug::Progress progress(querySize);

#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Sequence seq(qdbr->getMaxSeqLen(), querySeqType, kmerSubMat, kmerSize, spacedKmer, aaBiasCorrection, true,
                     spacedKmerPattern);
        QueryMatcher matcher(indexTable, sequenceLookup, kmerSubMat, ungappedSubMat,
                             kmerThr, kmerSize, dbSize, std::max(tdbr->getMaxSeqLen(), qdbr->getMaxSeqLen()),
                             maxResListLen, aaBiasCorrection, aaBiasCorrectionScale,
                             diagonalScoring, minDiagScoreThr, takeOnlyBestKmer,
                             targetSeqType == Parameters::DBTYPE_NUCLEOTIDES);

        if (seq.profile_matrix != nullptr) {
            matcher.setProfileMatrix(seq.profile_matrix);
        } else if (_3merSubMatrix.isValid() && _2merSubMatrix.isValid()) {
            matcher.setSubstitutionMatrix(&_3merSubMatrix, &_2merSubMatrix);
        } else {
            matcher.setSubstitutionMatrix(nullptr, nullptr);
        }

        if (taxonomyHook != nullptr) {
            matcher.setQueryMatcherHook(taxonomyHook);
        }

        char buffer[128];
        std::string result;
        result.reserve(1000000);

#pragma omp for schedule(dynamic, 1) reduction (+: kmersPerPos, resSize, dbMatches, doubleMatches, querySeqLenSum, diagonalOverflow)
        for (size_t id = queryFrom; id < queryFrom + querySize; id++) {
            progress.updateProgress();
            // get query sequence
            char *seqData = qdbr->getData(id, thread_idx);
            unsigned int qKey = qdbr->getDbKey(id);
            seq.mapSequence(id, qKey, seqData, qdbr->getSeqLen(id));
            size_t targetSeqId = UINT_MAX;
            if (sameQTDB || includeIdentical) {
                targetSeqId = tdbr->getId(seq.getDbKey());
                // only the corresponding split should include the id (hack for the hack)
                if (targetSeqId >= dbFrom && targetSeqId < (dbFrom + dbSize) && targetSeqId != UINT_MAX) {
                    targetSeqId = targetSeqId - dbFrom;
                    if (targetSeqId > tdbr->getSize()) {
                        Debug(Debug::ERROR) << "targetSeqId: " << targetSeqId << " > target database size: "
                                            << tdbr->getSize() << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                } else {
                    targetSeqId = UINT_MAX;
                }
            }
            // calculate prefiltering results
            if (taxonomyHook != nullptr)
                taxonomyHook->setDbFrom(dbFrom);

            std::pair<hit_t *, size_t> prefResults = matcher.matchQuery(&seq, targetSeqId, targetSeqType ==
                                                                                           Parameters::DBTYPE_NUCLEOTIDES);
            size_t resultSize = prefResults.second;
            const auto queryLength = static_cast<float>(qdbr->getSeqLen(id));
            for (size_t i = 0; i < resultSize; i++) {
                hit_t *res = prefResults.first + i;
                // correct the 0 indexed sequence id again to its real identifier
                size_t targetSeqId1 = res->seqId + dbFrom;
                // replace id with key
                res->seqId = tdbr->getDbKey(targetSeqId1);
                if (UNLIKELY(targetSeqId1 >= tdbr->getSize())) {
                    Debug(Debug::WARNING) << "Wrong prefiltering result for query: " << qdbr->getDbKey(id) << " -> "
                                          << targetSeqId1 << "\t" << res->prefScore << "\n";
                }

                // TODO: check if this should happen when diagonalScoring == false
                if (covThr > 0.0 && (covMode == Parameters::COV_MODE_BIDIRECTIONAL
                                     || covMode == Parameters::COV_MODE_QUERY
                                     || covMode == Parameters::COV_MODE_LENGTH_SHORTER)) {
                    const auto targetLength = static_cast<float>(tdbr->getSeqLen(targetSeqId1));
                    if (!Util::canBeCovered(covThr, covMode, queryLength, targetLength)) {
                        continue;
                    }
                }

                // write prefiltering results to a string
                int len = QueryMatcher::prefilterHitToBuffer(buffer, *res);
                result.append(buffer, len);
            }
            tmpDbw.writeData(result.c_str(), result.length(), qKey, thread_idx);
            result.clear();

            // update statistics counters
            if (resultSize != 0) {
                notEmpty[id - queryFrom] = 1;
            }

            if (Debug::debugLevel >= Debug::INFO) {
                kmersPerPos += matcher.getStatistics()->kmersPerPos;
                dbMatches += matcher.getStatistics()->dbMatches;
                doubleMatches += matcher.getStatistics()->doubleMatches;
                querySeqLenSum += seq.L;
                diagonalOverflow += matcher.getStatistics()->diagonalOverflow;
                resSize += resultSize;
                realResSize += std::min(resultSize, maxResListLen);
                reslens[thread_idx]->emplace_back(resultSize);
            }
        } // step end
    }

    if (Debug::debugLevel >= Debug::INFO) {
        statistics_t stats(kmersPerPos / static_cast<double>(totalQueryDBSize),
                           dbMatches / totalQueryDBSize,
                           doubleMatches / totalQueryDBSize,
                           querySeqLenSum, diagonalOverflow,
                           resSize / totalQueryDBSize);

        size_t empty = 0;
        for (size_t id = 0; id < querySize; id++) {
            if (notEmpty[id] == 0) {
                empty++;
            }
        }

        printStatistics(stats, reslens, localThreads, empty, maxResListLen);
    }

    if (splitMode == Parameters::TARGET_DB_SPLIT && splits == 1) {
#ifdef HAVE_MPI
        // if a mpi rank processed a single split, it must have it merged before all ranks can be united
        tmpDbw.close(true);
#else
        tmpDbw.close(merge);
#endif
    } else {
        tmpDbw.close(merge);
    }

    // sort by ids
    // needed to speed up merge later on
    // sorts this datafile according to the index file
    if (splitMode == Parameters::TARGET_DB_SPLIT && splits > 1) {
        // free memory early since the merge might need quite a bit of memory
        if (indexTable != nullptr) {
            delete indexTable;
            indexTable = nullptr;
        }
        if (sequenceLookup != nullptr) {
            delete sequenceLookup;
            sequenceLookup = nullptr;
        }
        DBReader<unsigned int> resultReader(tmpDbw.getDataFileName(), tmpDbw.getIndexFileName(), threads,
                                            DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        resultReader.open(DBReader<unsigned int>::NOSORT);
        resultReader.readMmapedDataInMemory();
        const std::pair<std::string, std::string> tempDb = Util::databaseNames((resultDB + "_tmp"));
        DBWriter resultWriter(tempDb.first.c_str(), tempDb.second.c_str(), localThreads, compressed,
                              Parameters::DBTYPE_PREFILTER_RES);
        resultWriter.open();
        resultWriter.sortDatafileByIdOrder(resultReader);
        resultWriter.close(true);
        resultReader.close();
        DBReader<unsigned int>::removeDb(resultDB);
        DBReader<unsigned int>::moveDb(tempDb.first, resultDB);
    }

    for (size_t i = 0; i < localThreads; i++) {
        reslens[i]->clear();
        delete reslens[i];
    }
    delete[] reslens;
    delete[] notEmpty;

    return true;
}
