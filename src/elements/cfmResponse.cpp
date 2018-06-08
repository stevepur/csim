//
//  celem.cpp
//  csim
//
//  Created by steve on 5/24/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include "cfmResponse.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>
#include "../lib/csim_lib.hpp"

cfmResponse::cfmResponse() {
}

cfmResponse::cfmResponse(initCommandSet*& cmdBlock) {
    init(cmdBlock);
}

efield* cfmResponse::execute(efield* E, celem* prev, celem* next, double time) {
    std::complex<double> i1(0, 1);
    int sliceSize = cfmResponseData->cfmResponseData[0][0]->n_rows*cfmResponseData->cfmResponseData[0][0]->n_cols;
//    std::cout << "executed binary cfmResponse " << name << std::endl;

    pre_execute(E, prev, next, time);
    
    for (int s=0; s<cfmResponseData->cfmResponseData.size(); s++) {
        for (int p=0; p<cfmResponseData->cfmResponseData[s].size(); p++) {
            // first collect ns slices into one big nr*ns x nc matrix
            for (int k=0; k<cfmResponseData->cfmResponseData[s][p]->n_slices; k++) {
                respMat(arm::span::all, arma::span(k*sliceSize,(k+1)*slizeSize - 1)) = cfmResponseData->cfmResponseData[s][p]->slice(k);
                
                // compute the phase rotation
                phaseMat(arma::span(k*sliceSize,(k+1)*slizeSize - 1)) = exp(-4*i1*pi*sagVec/E->lambdaData[i].get_wavelength(k));
            }
            
            
            // compute E for each column as the matrix product
            for (int k=0; k<E->E[s][p]->n_slices; k++) {
                E->E[s][p]->slice(k) = real(E->E[s][p]->slice(k)) % cfmResponseMat + i1*(imag(E->E[s][p]->slice(k)) % cfmResponseMat); // element-wise multiplication
            }
            // put back in cube order
            for (int k=0; k<E->E[s][p]->n_slices; k++) {
                E->E[s][p]->slice(k) = real(E->E[s][p]->slice(k)) % cfmResponseMat + i1*(imag(E->E[s][p]->slice(k)) % cfmResponseMat); // element-wise multiplication
            }
        }
    }
    
    post_execute(E, prev, next, time);
    
    return E;
}

void cfmResponse::init(const char *filename) {
    load_mat(filename, cfmResponseMat);
}

void cfmResponse::init(initCommandSet*& cmdBlock) {
//    std::cout << "initing a cfmResponse" << std::endl;
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
    post_init();
}

void cfmResponse::set(std::string fieldName, const char *arg) {
    bool found = celem::set(fieldName, arg);
    if (fieldName == "cfmResponse")
        ;
    else if (fieldName == "filename") {
        // arg is the .fits filename that contains the cfmResponse definition
        init(arg);
    } else if (!found)
        std::cout << "!!! cfmResponse bad set field name: " << fieldName << std::endl;
}

void cfmResponse::draw(const char *title) {
    char str[200];
    sprintf(str, "%s cfmResponse %s", title, name);
    draw_mat(cfmResponseMat, str, "rainbow");
}

