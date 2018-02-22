//
//  coronagraph.cpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include <iostream>
#include <string>
#include <assert.h>
#include "telescope.hpp"
#include "../lib/csim_parser.hpp"

telescope *globalTelescope = NULL;

telescope::telescope(initCommandSet*& cmdBlocks) {
    init(cmdBlocks);
}

void telescope::init(initCommandSet*& cmdBlocks) {
    std::cout << "initing a telescope" << std::endl;
    
    for (int c=0; c<cmdBlocks->commandList.size(); c++) {
        set(cmdBlocks->commandList[c]->getCmdStr(),
            cmdBlocks->commandList[c]->getArgStr());
    }
    assert(primaryDiameter != 0 & primaryfRatio != 0);
    primaryfLength = primaryfRatio*primaryDiameter;
}

void telescope::set(std::string fieldName, const char *value) {
    if (fieldName == "telescope")
        ;
    else if (fieldName == "primaryDiameter")
        primaryDiameter = atof(value);
    else if (fieldName == "primaryFRatio" || fieldName == "primaryfRatio")
        primaryfRatio = atof(value);
    else
        std::cout << "!!!!! telescope set: unknown fieldName: " << fieldName << std::endl;
}

double telescope::get(std::string fieldName) {
    if (fieldName == "primaryDiameter")
        return primaryDiameter;
    else if (fieldName == "primaryFRatio" || fieldName == "primaryfRatio")
    return primaryfRatio;
    else if (fieldName == "primaryFLength" || fieldName == "primaryfLength")
    return primaryfLength;
    else
        std::cout << "!!!!! telescope get: unknown fieldName: " << fieldName << std::endl;
}

void telescope::print(const char *header) {
    std::cout << "I'm a telescope " << header << std::endl;
    std::cout << "primaryDiameter = " << primaryDiameter << std::endl;
    std::cout << "primaryfRatio = " << primaryfRatio << std::endl;
}
