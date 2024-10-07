//
// Created by minghang on 10/7/24.
//

#ifndef FOLDSEEK_COMBINEDPREFILTERING_H
#define FOLDSEEK_COMBINEDPREFILTERING_H


#include "Prefiltering.h"

class CombinedPrefiltering: public Prefiltering {
public:
    using Prefiltering::Prefiltering;

    bool runSplit(const std::string &resultDB, const std::string &resultDBIndex, size_t split, bool merge) override;
};


#endif //FOLDSEEK_COMBINEDPREFILTERING_H
