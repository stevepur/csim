//
//  celem.cpp
//  csim
//
//  Created by steve on 5/24/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include "mask.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>
#include "../lib/csim_lib.hpp"

mask::mask() {
}

mask::mask(initCommandSet*& cmdBlock) {
    init(cmdBlock);
}

efield* mask::execute(efield* E, celem* prev, celem* next, double time) {
    std::complex<double> i1(0, 1);
//    std::cout << "executed binary mask " << name << std::endl;

    pre_execute(E, prev, next, time);
    if (maskMat.n_rows != E->E[0][0]->n_rows) {
//        draw_mat(maskMat, "mask before downsampling");
        downsample_via_convolution(maskMat, maskMat, E->initArrayGeometry, E->arrayGeometry);
//        draw_mat(maskMat, "mask after downsampling");
        cxMaskMat.set_size(size(maskMat));
        cxMaskMat.set_real(maskMat);
        cxMaskMat.set_imag(maskMat);
    }
//    save_mat("resampled mask.fits", maskMat);
    
    for (int s=0; s<E->E.size(); s++) {
        for (int p=0; p<E->E[s].size(); p++) {
            #pragma omp parallel for
            for (int k=0; k<E->E[s][p]->n_slices; k++) {
//                E->E[s][p]->slice(k) = real(E->E[s][p]->slice(k)) % maskMat + i1*(imag(E->E[s][p]->slice(k)) % maskMat); // element-wise multiplication
                E->E[s][p]->slice(k) = E->E[s][p]->slice(k) % cxMaskMat; // element-wise multiplication
            }
        }
    }
    
    post_execute(E, prev, next, time);
    
    return E;
}

void mask::init(const char *filename) {
    load_mat(filename, maskMat);
    cxMaskMat.set_size(size(maskMat));
    cxMaskMat.set_real(maskMat);
    cxMaskMat.set_imag(maskMat);
}

void mask::init(initCommandSet*& cmdBlock) {
//    std::cout << "initing a mask" << std::endl;
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
    post_init();
}

void mask::set(std::string fieldName, const char *arg) {
    bool found = celem::set(fieldName, arg);
    if (fieldName == "mask")
        ;
    else if (fieldName == "filename") {
        // arg is the .fits filename that contains the mask definition
        init(arg);
    } else if (!found)
        std::cout << "!!! mask bad set field name: " << fieldName << std::endl;
}

void mask::draw(const char *title) {
    char str[200];
    sprintf(str, "%s mask %s", title, name);
    draw_mat(maskMat, str, "rainbow");
}

