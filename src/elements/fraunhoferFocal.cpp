//
//  fraunhoferFocal.cpp
//  csim
//
//  Created by steve on 5/24/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include "fraunhoferFocal.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>
#include "../telescope/telescope.hpp"

fraunhoferFocal::fraunhoferFocal() {
}

fraunhoferFocal::fraunhoferFocal(initCommandSet*& cmdBlock) {
    init(cmdBlock);
}

efield* fraunhoferFocal::execute(efield* E, celem* prev, celem* next, double time) {
//    std::cout << "executed fraunhoferFocal " << name << std::endl;
    double z;
    
    pre_execute(E, prev, next, time);
    
    int Nx = arrayGeometry.pixelX.n_elem;
    int Ny = arrayGeometry.pixelY.n_elem;
    
    assert(E->beamRadiusPhysical > 0);
    
//    E->print("in fraunhoferFocal::execute: ");
    double *lambdaFocalLength = new double[E->E[0][0]->n_slices];
//    std::cout << "lambdaFocalLength :" << std::endl;
    for (int i=0; i<E->E[0][0]->n_slices; i++) {
        lambdaFocalLength[i] = E->lambdaData[i].get_wavelength()*2*E->beamRadiusPhysical*globalTelescope->get("primaryfRatio");
//        std::cout << " lambda = " << E->lambdaData[i].get_wavelength() << ", lFl = " << lambdaFocalLength[i] << "; ";
    }
//    std::cout << std::endl;
//    E->arrayGeometry.print("fraunhoferFocal input E geometry: ");
//    arrayGeometry.print("fraunhoferFocal output E geometry: ");
    
    propZoomFft.init(E->arrayGeometry, arrayGeometry, lambdaFocalLength, E->E[0][0]->n_slices);
    
    for (int s=0; s<E->E.size(); s++) {
        for (int p=0; p<E->E[s].size(); p++) {
            *(E->E[s][p]) = propZoomFft.execute(*(E->E[s][p]));
        }
    }

    // E has now changed geometry to this array geometry
    E->arrayGeometry = arrayGeometry;
    
    post_execute(E, prev, next, time);
    
    return E;
}

void fraunhoferFocal::set_geometry(void) {
    
//    scaleFlD = globalTelescope->get("primaryfLength")*referenceLambda/globalTelescope->get("primaryDiameter");
    scaleFlD = globalTelescope->compute_loD(referenceLambda);
//    arrayGeometry.set_geometry(scaleFlD, samplesPerflD, FOVflD);
    arrayGeometry.set_geometry(ceil(FOVflD*samplesPerflD), scaleFlD/samplesPerflD);

}

void fraunhoferFocal::init(initCommandSet*& cmdBlock) {
    std::cout << "initing a fraunhoferFocal" << std::endl;
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
    set_geometry();
    post_init();
}

void fraunhoferFocal::set(std::string fieldName, const char *arg) {
    bool found = celem::set(fieldName, arg);
    if (fieldName == "fraunhoferFocal")
        ;
    else if (fieldName == "referenceLambda") {
        referenceLambda = atof(arg);
    }
    else if (fieldName == "samplesPerflD") {
        samplesPerflD = atof(arg);
    }
    else if (fieldName == "FOVflD") {
        FOVflD = atof(arg);
    }
    else if (!found)
        std::cout << "!!! fraunhoferFocal bad set field name: " << fieldName << std::endl;
}

void fraunhoferFocal::print(const char *hdr) {
    std::cout << "fraunhoferFocal " << hdr << std::endl;
    std::cout << "referenceLambda = " << referenceLambda << std::endl;
    std::cout << "samplesPerflD = " << samplesPerflD << std::endl;
    std::cout << "FOVflD = " << FOVflD << std::endl;
    std::cout << "scaleFlD = " << scaleFlD << std::endl;
    std::cout << "focalLength = " << focalLength << std::endl;
}

