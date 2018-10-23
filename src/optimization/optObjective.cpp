//
//  optObjective_hpp.cpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>

#include "armadillo"
#include "optObjective.hpp"

optHexRespObjective::optHexRespObjective() {
}

optHexRespObjective::optHexRespObjective(initCommandSet*& cmdBlocks) {
    init(cmdBlocks);
}

void optHexRespObjective::execute(arma::mat& optData, arma::rowvec& objVal) {
    std::complex<double> i1(0, 1);
    
//    std::cout << "starting optHexRespObjective::execute" << std::endl;
    // need to add a column of zeros to the data matrix
    arma::mat zmat = arma::zeros<arma::mat>(1, optData.n_cols);
//    std::cout << "made zmat: " << size(zmat) << std::endl;
    arma::mat fullOptData = join_cols(optData, zmat);
//    std::cout << "made fullOptData: " << size(fullOptData) << std::endl;
    arma::mat pixEnergy = arma::zeros(response->M[0][0]->slice(0).n_rows, fullOptData.n_cols);
//    std::cout << "made pixEnergy: " << size(pixEnergy) << std::endl;
    rotatedPhases.set_size(size(fullOptData));
    for (int s=0; s<response->M.size(); s++) {
        for (int p=0; p<response->M[s].size(); p++) {
            for (int sl=0; sl<response->M[s][p]->n_slices; sl++) {
                double a = phaseSign*4*M_PI/response->wavelengths[sl];
                #pragma omp parallel for
                for (int i=0; i<fullOptData.n_cols; i++)
                    rotatedPhases(arma::span::all, i) = exp(i1*a*fullOptData(arma::span::all, i));
//                std::cout << "computed rotatedPhases: " << size(rotatedPhases) << std::endl;
//                std::cout << "mean(mean(rotatedPhases)): " << mean(mean(rotatedPhases)) << std::endl;

                arma::cx_mat E = response->M[s][p]->slice(sl) * rotatedPhases;
//                std::cout << "mean(mean(E)): " << mean(mean(E)) << std::endl;
//                std::cout << "computed E: " << size(E) << std::endl;
                pixEnergy += real(E % arma::conj(E));
//                std::cout << "mean(mean(pixEnergy)): " << mean(mean(pixEnergy)) << std::endl;
//                std::cout << "computed pixEnergy: " << size(pixEnergy) << std::endl;
            }
        }
    }
//    std::cout << "completed computing pix energy" << std::endl;
    switch (objectiveType) {
        case OBJECTIVE_MEAN:
            objVal = mean(pixEnergy)/response->calibIntensity;
            break;
            
        case OBJECTIVE_MAX:
            objVal = max(pixEnergy)/response->calibIntensity;
            break;
            
        default:
            assert(NULL);
    }
//    std::cout << "computed objVal: " << size(objVal) << std::endl;
}

void optHexRespObjective::init(initCommandSet*& cmdBlock) {
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
}

void optHexRespObjective::set(std::string fieldName, const char *arg) {
    if (fieldName == "optHexRespObjective")
        ;
    else if (fieldName == "responseDirectory") {
        // arg is two strings separated by a comma,
        // the first string is the name of the component to be optimized
        // the second string is the name of the data to be optimized
        response = new responseData;
        response->set((std::string) "outputDirectory", arg);
        response->load();
    }
    else if (fieldName == "objectiveType") {
        // arg is a string
        if (!strcmp(arg, "mean")) {
            objectiveType = OBJECTIVE_MEAN;
        }
        else if (!strcmp(arg, "max")) {
            objectiveType = OBJECTIVE_MAX;
        }
        else
            std::cout << "!!! differentialEvolutionOptimizer initType bad arg: " << arg << std::endl;
    }
    else if (fieldName == "phaseSign") {
        // arg is a single double
        phaseSign = atof(arg);
    }
    else
        std::cout << "!!! optHexRespObjective bad set field name: " << fieldName << std::endl;
    
}

int optHexRespObjective::get_opt_vector_size(void) {
    assert(response);
    
    return response->M[0][0]->slice(0).n_cols;
}

void optHexRespObjective::print(const char *hdr) {
    if (hdr != NULL)
        std::cout << hdr << ":" << std::endl;
    std::cout << "phaseSign = " << phaseSign << std::endl;
    response->print("response data");
    
}
