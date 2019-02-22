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
            if (padFactor == 1)
                propagator.execute(*(E->E[s][p]), lambda, E->arrayRadiusPhysical, z);
            else {
//                std::cout << "==================== padding fresnel ====================" << std::endl;
                int nRowsE = E->E[s][p]->n_rows;
                int nColsE = E->E[s][p]->n_cols;
                arma::cx_cube paddedE;
                paddedE = arma::zeros<arma::cx_cube>(2*nRowsE, 2*nColsE, E->E[0][0]->n_slices);
                // embed E in the zero-padded cube
                paddedE(arma::span(nRowsE/2, 3*nRowsE/2-1), arma::span(nColsE/2, 3*nColsE/2-1), arma::span::all) = *(E->E[s][p]);
                propagator.execute(paddedE, lambda, padFactor*E->arrayRadiusPhysical, z);
                for (int sl=0; sl<E->E[s][p]->n_slices; sl++)
                    E->E[s][p]->slice(sl) = paddedE(arma::span(nRowsE/2, 3*nRowsE/2-1), arma::span(nColsE/2, 3*nColsE/2-1), arma::span(sl, sl));
            }
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
    else if (fieldName == "padFactor") {
        padFactor = atoi(arg);
        std::cout << "========= set padFactor to " << padFactor << std::endl;
    }
    else if (!found)
        std::cout << "!!! fresnelPropagator bad set field name: " << fieldName << std::endl;
}

