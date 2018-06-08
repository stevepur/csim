//
//  mirror.cpp
//  csim
//
//  Created by steve on 5/24/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include "mirror.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>
#include "../lib/csim_lib.hpp"

mirror::mirror() {
}

mirror::mirror(initCommandSet*& cmdBlock) {
    init(cmdBlock);
}

efield* mirror::execute(efield* E, celem* prev, celem* next, double time) {
//    std::cout << "executed a mirror " << name << std::endl;
    std::complex<double> i1(0, 1);

    pre_execute(E, prev, next, time);
    
    for (int s=0; s<E->E.size(); s++) {
        for (int p=0; p<E->E[s].size(); p++) {
            for (int k=0; k<E->E[s][p]->n_slices; k++) {
                double lambda = E->lambdaData[k].lambda;
                E->E[s][p]->slice(k) = E->E[s][p]->slice(k) % exp((-mirrorSign*2*2*M_PI/lambda)*i1*mirrorMat); // element-wise multiplication
            }
        }
    }
    
    post_execute(E, prev, next, time);
    
    return E;
}

void mirror::init(const char *filename) {
    load_mat(filename, mirrorMat);
}

void mirror::init(initCommandSet*& cmdBlock) {
//    std::cout << "initing a mirror" << std::endl;
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
    
    post_init();
}

void mirror::set(std::string fieldName, const char *arg) {

    bool found = celem::set(fieldName, arg);
    if (fieldName == "mirror")
        ;
    else if (fieldName == "filename") {
        // arg is the .fits filename that contains the mirror definition
        init(arg);
    }
    else if (fieldName == "sign") {
        double mSign = atof(arg);
        if (mSign == -1 | mSign == 1)
            mirrorSign = mSign;
        else
            std::cout << "!!!! illegal propagationSign, ignoring" << std::endl;
    }
    else if (!found)
        std::cout << "!!! mirror bad set field name: " << fieldName << std::endl;
}

void mirror::draw(const char *title) {
    char str[200];
    sprintf(str, "%s mirror %s", title, name);
    draw_mat(mirrorMat, str, "gray");
}

