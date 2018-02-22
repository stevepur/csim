//
//  celem.cpp
//  csim
//
//  Created by steve on 5/24/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include "breakExecution.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>
#include "../lib/csim_lib.hpp"

breakExecution::breakExecution() {
}

breakExecution::breakExecution(initCommandSet*& cmdBlock) {
    init(cmdBlock);
}

efield* breakExecution::execute(efield* E, celem* prev, celem* next, double time) {
    assert(NULL);
    
    return E;
}

void breakExecution::init(initCommandSet*& cmdBlock) {
    std::cout << "initing a breakExecution" << std::endl;
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
    post_init();
}

void breakExecution::set(std::string fieldName, const char *arg) {
    bool found = celem::set(fieldName, arg);
    if (fieldName == "breakExecution")
        ;
    else if (!found)
        std::cout << "!!! breakExecution bad set field name: " << fieldName << std::endl;
}

