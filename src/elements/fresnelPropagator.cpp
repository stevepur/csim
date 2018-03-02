//
//  fresnelPropagator.cpp
//  csim
//
//  Created by steve on 5/24/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include "fresnelPropagator.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>
#include "../lib/csim_lib.hpp"

fresnelPropagator::fresnelPropagator() {
}

fresnelPropagator::fresnelPropagator(initCommandSet*& cmdBlock) {
    init(cmdBlock);
}

efield* fresnelPropagator::execute(efield* E, celem* prev, celem* next, double time) {
//    std::cout << "executed fresnelPropagator " << name << std::endl;
    std::complex<double> i1(0, 1);
    double z;
    
    assert(prev != NULL & next != NULL);
    
    pre_execute(E, prev, next, time);
    
    z = propagationSign*abs(next->position - prev->position);
//    std::cout << "propagation distance " << z << std::endl;

    double *lambda = new double[E->E[0][0]->n_slices];
    for (int i=0; i<E->E[0][0]->n_slices; i++)
        lambda[i] = E->lambdaData[i].lambda;

    for (int s=0; s<E->E.size(); s++) {
        for (int p=0; p<E->E[s].size(); p++) {
            propagator.execute(*(E->E[s][p]), lambda, E->arrayRadiusPhysical, z);
        }
    }
    
    post_execute(E, prev, next, time);
    
    return E;
}

void fresnelPropagator::init(initCommandSet*& cmdBlock) {
//    std::cout << "initing a fresnelPropagator" << std::endl;
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
    post_init();
}

void fresnelPropagator::set(std::string fieldName, const char *arg) {
    bool found = celem::set(fieldName, arg);
    if (fieldName == "fresnelPropagator")
        ;
    else if (fieldName == "propagationSign") {
        double pSign = atof(arg);
        if (pSign == -1 | pSign == 1)
            propagationSign = pSign;
        else
            std::cout << "!!!! illegal propagationSign, ignoring" << std::endl;
    }
    else if (!found)
        std::cout << "!!! fresnelPropagator bad set field name: " << fieldName << std::endl;
}

