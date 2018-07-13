//
//  efiled.cpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//


#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "responseData.hpp"
#include "../lib/csim_lib.hpp"
#include "../telescope/telescope.hpp"


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

responseData::responseData() {
}

responseData::responseData(int nRows, int nColumns, int nSlices, int nPolarizations, int nSources) {
    init(nRows, nColumns, nSlices, nPolarizations, nSources);
}

responseData::responseData(responseData& in) {
    
    if (in.name != NULL) {
        name = new char[strlen(in.name)+1];
        strcpy(name, in.name);
    }
    
    verbose = in.verbose;
    
    outputDirectory = new char[strlen(in.outputDirectory)+1];
    strcpy(outputDirectory, in.outputDirectory);
    
    for (int s=0; s<in.M.size(); s++) {
        for (int p=0; p<in.M[s].size(); p++) {
            arma::cx_cube *newM = new arma::cx_cube;
            *newM = *(in.M[s][p]);
            polarizations tmpP;
            tmpP.push_back(newM);
            M.push_back(tmpP);
        }
    }
}

responseData::responseData(char *dirName) {
}

responseData::responseData(initCommandSet*& cmdBlock) {
    init(cmdBlock);
}

responseData::~responseData(void) {
    for (int s=0; s<M.size(); s++) {
        for (int p=0; p<M[s].size(); p++) {
            delete M[s][p];
        }
    }
}


responseData& responseData::operator=(const responseData& in) {
    
    if (in.name != NULL) {
        name = new char[strlen(in.name)+1];
        strcpy(name, in.name);
    }
    
    verbose = in.verbose;
    
    outputDirectory = new char[strlen(in.outputDirectory)+1];
    strcpy(outputDirectory, in.outputDirectory);
    
    for (int s=0; s<in.M.size(); s++) {
        for (int p=0; p<in.M[s].size(); p++) {
            arma::cx_cube *newM = new arma::cx_cube;
            *newM = *(in.M[s][p]);
            polarizations tmpP;
            tmpP.push_back(newM);
            M.push_back(tmpP);
        }
    }
    
    return *this;
}

void responseData::init(int nRows, int nColumns, int nSlices, int nPolarizations, int nSources) {
    // fill in all the sources with one polarization (for now)
    for (int i=0; i<nSources; i++) {
        arma::cx_cube *newM = new arma::cx_cube(nRows, nColumns, nSlices, arma::fill::zeros);
        polarizations tmpP;
        tmpP.push_back(newM);
        M.push_back(tmpP);
    }
}

void responseData::init(initCommandSet*& cmdBlock) {
    std::cout << "initing a responseData" << std::endl;
    
    const char *arg = "";
    outputDirectory = new char[strlen(arg)+1];
    strcpy(outputDirectory, arg);

    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
    print();
}

void responseData::set(std::string fieldName, const char *arg) {
    if (fieldName == "responseData")
        ;
    else if (fieldName == "setSize") {
        // arg is three integers, comma separated
        int nRows;
        int nCols;
        int nSlices;
        int nPolarizations;
        int nSources;
        sscanf(arg, "%d, %d, %d, %d, %d", &nRows, &nCols, &nSlices, &nPolarizations, &nSources);
        init(nRows, nCols, nSlices, nPolarizations, nSources);
    }
    else if (fieldName == "outputDirectory") {
        // arg is a string
        outputDirectory = new char[strlen(arg)+1];
        strcpy(outputDirectory, arg);
        if (access(outputDirectory, 0) == -1) {
            std::cout << "making " << outputDirectory << std::endl;
            mkdir(outputDirectory, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
        }
    }
    else if (fieldName == "name") {
        // arg is a string
        name = new char[strlen(arg)+1];
        strcpy(name, arg);
    }
    else
        std::cout << "!!!!! responseData set: unknown fieldName: " << fieldName << std::endl;
}

void responseData::set(std::string fieldName, double arg) {
    if (fieldName == "calibIntensity") {
        // arg is a float
        calibIntensity = arg;
    }
    else
        std::cout << "!!!!! responseData set: unknown fieldName: " << fieldName << std::endl;
}

void responseData::set_size(int nRows, int nColumns, int nSlices) {
    for (int s=0; s<M.size(); s++) {
        for (int p=0; p<M[s].size(); p++) {
            M[s][p]->set_size(nRows, nColumns, nSlices);
        }
    }
}

void responseData::set_wavelengths(lambdaDataClass *lambdaData) {
    assert(M[0][0]->n_slices > 0);
    
    wavelengths.set_size(M[0][0]->n_slices);
    for (int i=0; i<M[0][0]->n_slices; i++)
        wavelengths[i] = lambdaData[i].get_wavelength();
}

void responseData::save(void) {
    
    std::string fname = (std::string)outputDirectory + "/responseParams.txt";
    FILE *fid = fopen(fname.c_str(), "w");
    fprintf(fid, "nRows = %lu, nCols = %lu, nSlices = %lu, nPolarizations = %lu, nSources = %lu\n",
            M[0][0]->n_rows, M[0][0]->n_cols, M[0][0]->n_slices, M[0].size(), M.size());
    fclose(fid);
    for (int s=0; s<M.size(); s++) {
        for (int p=0; p<M[s].size(); p++) {
            fname = (std::string)outputDirectory + "/responseData_s" + std::to_string(s) + "_p" + std::to_string(p);
            save_cube(fname.c_str(), *(M[s][p]), false);
        }
    }
    save_problem_params();
}

void responseData::load(void) {
    int nRows = 0;
    int nColumns = 0;
    int nSlices = 0;
    int nPolarizations = 0;
    int nSources = 0;
    
    std::string fname = (std::string)outputDirectory + "/responseParams.txt";
    FILE *fid = fopen(fname.c_str(), "r");
    int d = fscanf(fid, "nRows = %d, nCols = %d, nSlices = %d, nPolarizations = %d, nSources = %d\n",
            &nRows, &nColumns, &nSlices, &nPolarizations, &nSources);
    fclose(fid);
    
//    printf("nRows = %d, nCols = %d, nSlices = %d, nPolarizations = %d, nSources = %d\n",
//           nRows, nColumns, nSlices, nPolarizations, nSources);
    init(nRows, nColumns, nSlices, nPolarizations, nSources);

    for (int s=0; s<M.size(); s++) {
        for (int p=0; p<M[s].size(); p++) {
            fname = (std::string)outputDirectory + "/responseData_s" + std::to_string(s) + "_p" + std::to_string(p) + ".fits";
            load_cube(fname.c_str(), *(M[s][p]), false);
        }
    }
    
    read_problem_params();
}

void responseData::save_problem_params(void) {
    std::string fname = (std::string)outputDirectory + "/problemParams.txt";
    FILE *fid = fopen(fname.c_str(), "w");
    
    fprintf(fid, "calibMaxIntensity: %lf\n", calibIntensity);
    fprintf(fid, "wavelengths: ");
    for (int i=0; i<wavelengths.n_elem; i++)
        fprintf(fid, "%g ", wavelengths[i]);
    fprintf(fid, "\n");
    
    fclose(fid);
}

void responseData::read_problem_params(void) {
    char wavelengthStr[1000];
    char singleWavelengthStr[100];
    int retVal = 1;
    
    std::string fname = (std::string)outputDirectory + "/problemParams.txt";
    FILE *fid = fopen(fname.c_str(), "r");

    double ci;
    int d = fscanf(fid, "calibMaxIntensity: %lf\n", &ci);
    calibIntensity = ci;
    retVal = fscanf(fid, "wavelengths: %s", wavelengthStr);
//    std::cout << "wavelengthStr: " << wavelengthStr << std::endl;
    while (retVal != EOF) {
        retVal = fscanf(fid, "%s ", singleWavelengthStr);
        strcat(wavelengthStr, " ");
        strcat(wavelengthStr, singleWavelengthStr);
//        std::cout << "wavelengthStr: " << wavelengthStr << std::endl;
    }
    fclose(fid);
    
//    std::cout << "wavelengthStr: " << wavelengthStr << std::endl;

    retVal = 1;
    int wCount = 0;
    int cPos = 0;
    double waveLengthTemp[100];
    while (retVal != EOF & cPos < strlen(wavelengthStr)) {
        int ncPos;
        retVal = sscanf(wavelengthStr+cPos, "%lf %n", &waveLengthTemp[wCount], &ncPos);
        cPos += ncPos;
//        std::cout << "wavelength " << wCount << ": " << waveLengthTemp[wCount] << ", retVal = " << retVal << ", cPos = " << cPos << std::endl;
        wCount++;
        if (wCount >= 1000)
            assert(NULL);
    }
    if (wCount > 1)
        wCount--;
    std::cout << wCount << " wavelengths read" << std::endl;
    wavelengths.set_size(wCount);
    for (int i=0; i<wCount; i++)
        wavelengths(i) = waveLengthTemp[i];
}


void responseData::print(const char *header) {
    if (header != NULL)
        std::cout << header << ":" << std::endl;
    std::cout << "num sources = " << M.size() << std::endl;
    std::cout << "num polarizations = " << M[0].size() << std::endl;
    std::cout << "size for each source and polarization = " << size(*M[0][0]) << std::endl;
    std::cout << "outputDirectory = " << outputDirectory << std::endl;
    std::cout << "calibIntensity: " << calibIntensity << std::endl;
    std::cout << "wavelengths:" << std::endl;
    for (int i=0; i<wavelengths.n_elem; i++)
        std::cout << wavelengths(i) << std::endl;
}

