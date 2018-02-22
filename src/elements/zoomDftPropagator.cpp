//
//  zoomDftPropagator.cpp
//  csim
//
//  Created by steve on 5/24/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include "zoomDftPropagator.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>

zoomDftPropagator::zoomDftPropagator() {
}

zoomDftPropagator::zoomDftPropagator(initCommandSet*& cmdBlock) {
    init(cmdBlock);
}

efield* zoomDftPropagator::execute(efield* E, celem* prev, celem* next, double time) {
    std::cout << "executed zoomDftPropagator " << name << std::endl;
    double z;
    
    print("execute parameters:");
    pre_execute(E, prev, next, time);

    propZoomDft.init(E->arrayGeometry.pixelX.n_elem, E->arrayGeometry.pixelY.n_elem, zoomFactor, dir);
    
    for (int s=0; s<E->E.size(); s++) {
        for (int p=0; p<E->E[s].size(); p++) {
            *(E->E[s][p]) = propZoomDft.execute(*(E->E[s][p]));
        }
    }
    
    post_execute(E, prev, next, time);
    
    return E;
}

void zoomDftPropagator::init(initCommandSet*& cmdBlock) {
    std::cout << "initing a zoomDftPropagator" << std::endl;
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
    post_init();
}

void zoomDftPropagator::set(std::string fieldName, const char *arg) {
    bool found = celem::set(fieldName, arg);
    if (fieldName == "zoomDftPropagator")
    ;
    else if (fieldName == "dir") {
        dir = atoi(arg);
    }
    else if (fieldName == "zoomFactor") {
        zoomFactor = atof(arg);
    }
    else if (!found)
        std::cout << "!!! zoomDftPropagator bad set field name: " << fieldName << std::endl;
}

void zoomDftPropagator::print(const char *hdr) {
    std::cout << "zoomDftPropagator " << hdr << std::endl;
    std::cout << "dir = " << dir << std::endl;
    std::cout << "zoomFactor = " << zoomFactor << std::endl;
}

