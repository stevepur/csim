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
#include "efield.hpp"
#include "../lib/csim_lib.hpp"
#include "../telescope/telescope.hpp"

#define MAXNWAVELENGTHS 100

efield *initialEfield = NULL;

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

lambdaDataClass::lambdaDataClass(void) {
}

lambdaDataClass::lambdaDataClass(double wavelength) {
    set_wavelength(wavelength);
}

void lambdaDataClass::set_wavelength(double wavelength) {
    assert(globalTelescope->get("primaryfLength") != 0 & globalTelescope->get("primaryDiameter") != 0);
    
    lambda = wavelength;
//    focalLengthLambdaOverD = globalTelescope->get("primaryfLength") * lambda / globalTelescope->get("primaryDiameter");
    focalLengthLambdaOverD = globalTelescope->compute_loD(lambda);
}

double lambdaDataClass::get_wavelength(void) {
    return lambda;
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

pointSourceClass::pointSourceClass(void) {
}

pointSourceClass::pointSourceClass(double f, double tx, double ty) {
    flux = f;
    tipX = tx;
    tiltY = ty;
}

void pointSourceClass::print(const char *hdr) {
    std::cout << hdr << "flux=" << flux << ", tipX=" << tipX << ", tiltY=" << tiltY << std::endl;
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

efield::efield() {
}

efield::efield(int nRows, int nColumns, int nLambda) {
    init();
}

efield::efield(efield& in) {
    
    if (in.name != NULL) {
        name = new char[strlen(in.name)+1];
        strcpy(name, in.name);
    }
    
    verbose = in.verbose;
    
    outputDirectory = new char[strlen(in.outputDirectory)+1];
    strcpy(outputDirectory, in.outputDirectory);
    
    for (int s=0; s<in.E.size(); s++) {
        for (int p=0; p<in.E[s].size(); p++) {
            arma::cx_cube *newE = new arma::cx_cube;
            *newE = *(in.E[s][p]);
            polarizations tmpP;
            tmpP.push_back(newE);
            E.push_back(tmpP);
        }
    }
    
    lambdaData = new lambdaDataClass[E[0][0]->n_slices];
    for (int s=0; s<E[0][0]->n_slices; s++)
        lambdaData[s] = in.lambdaData[s];
    
    pointSourceList = in.pointSourceList;
    
    beamRadiusPhysical = in.beamRadiusPhysical;
    pixelScale = in.pixelScale;
    arraySizePixels = in.arraySizePixels;
    beamSizePixels = in.beamSizePixels;
    arrayRadiusPhysical = in.arrayRadiusPhysical;
    
    arrayGeometry = in.arrayGeometry;
    initArrayGeometry = in.initArrayGeometry;
}

efield::efield(char *reFilename, char *imFilename) {
    init(reFilename, imFilename);
}

efield::efield(initCommandSet*& cmdBlock) {
    arma::cx_cube *newE = new arma::cx_cube;
    polarizations tmpP;
    tmpP.push_back(newE);
    E.push_back(tmpP);

    init(cmdBlock);
}

efield::~efield(void) {
    for (int s=0; s<E.size(); s++) {
        for (int p=0; p<E[s].size(); p++) {
            delete E[s][p];
        }
    }
    delete lambdaData;
}


efield& efield::operator=(const efield& in) {
    
    if (in.name != NULL) {
        name = new char[strlen(in.name)+1];
        strcpy(name, in.name);
    }
    
    verbose = in.verbose;
    
    outputDirectory = new char[strlen(in.outputDirectory)+1];
    strcpy(outputDirectory, in.outputDirectory);
    
    for (int s=0; s<in.E.size(); s++) {
        for (int p=0; p<in.E[s].size(); p++) {
            arma::cx_cube *newE = new arma::cx_cube;
            *newE = *(in.E[s][p]);
            polarizations tmpP;
            tmpP.push_back(newE);
            E.push_back(tmpP);
        }
    }
    
    lambdaData = new lambdaDataClass[E[0][0]->n_slices];
    for (int s=0; s<E[0][0]->n_slices; s++)
        lambdaData[s] = in.lambdaData[s];
    
    beamRadiusPhysical = in.beamRadiusPhysical;
    pixelScale = in.pixelScale;
    arraySizePixels = in.arraySizePixels;
    beamSizePixels = in.beamSizePixels;
    arrayRadiusPhysical = in.arrayRadiusPhysical;
    
    arrayGeometry = in.arrayGeometry;
    initArrayGeometry = in.initArrayGeometry;

    return *this;
}

void efield::init(double initValue) {
    for (int s=0; s<E.size(); s++) {
        for (int p=0; p<E[s].size(); p++) {
            E[s][p]->ones(size(*E[s][p]));
            *E[s][p] = (*E[s][p])*initValue;
        }
    }
}

// input flux is relative to central star flux = 1
void efield::init_point_sources(void) {
    std::complex<double> i1(0, 1);
    int startInitedSources = nInitedSources;
    for (int s=nInitedSources; s<pointSourceList.size()+startInitedSources; s++) {
        // does a slot for this source exist?
        if (s >= E.size()) {
            // add the source
            arma::cx_cube *newE = new arma::cx_cube;
            newE->ones(size(*E[0][0]));
            polarizations tmpP;
            tmpP.push_back(newE);
            E.push_back(tmpP);
        }
        
        // initialize the source for each polarization
        for (int p=0; p<E[s].size(); p++) {
            // set it to constant flux
            E[s][p]->ones(size(*E[s][p]));
            *E[s][p] *= pointSourceList[s].flux;
            for (int sl=0; sl<E[s][p]->n_slices; sl++) {
                E[s][p]->slice(sl) %= exp(i1*2.0*M_PI
                            * (2*arrayGeometry.pixelXX*pointSourceList[s].tipX/arrayGeometry.physicalSizeX
                               + 2*arrayGeometry.pixelYY*pointSourceList[s].tiltY/arrayGeometry.physicalSizeY));
            }
        }
        nInitedSources++;
    }
}

void efield::add_point_source(double flux, double tipX, double tiltY) {
    pointSourceClass *tmpS = new pointSourceClass(flux, tipX, tiltY);
    pointSourceList.push_back(*tmpS);
    pointSourceList[pointSourceList.size()-1].print("added point source ");
}

void efield::init(char *reFilename, char *imFilename) {
    void *data = NULL;
    int nDims = 0;
    long *dims = NULL;
    int dataType;
    arma::cube reMat;
    arma::cube imMat;
    
    load_cube(reFilename, reMat);
    load_cube(imFilename, imMat);

    arma::cx_cube tmpE(reMat.n_rows, reMat.n_cols, reMat.n_slices);
    tmpE.set_real(reMat);
    tmpE.set_imag(imMat);
    
    set_size(reMat.n_rows, reMat.n_cols, reMat.n_slices);

    // does a slot for this source exist?
    if (nInitedSources >= E.size()) {
        // add the source
        arma::cx_cube *newE = new arma::cx_cube;
        newE->ones(size(*E[0][0]));
        polarizations tmpP;
        tmpP.push_back(newE);
        E.push_back(tmpP);
    }
    // set the source
    E[nInitedSources][0]->set_size(reMat.n_rows, reMat.n_cols, reMat.n_slices);
    *E[nInitedSources][0] = tmpE;
    nInitedSources++;
}

void efield::init(initCommandSet*& cmdBlock) {
    std::cout << "initing an efield" << std::endl;
    
    const char *arg = "";
    outputDirectory = new char[strlen(arg)+1];
    strcpy(outputDirectory, arg);

    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
    set_optical_parameters();
    set_array_geometry();
    init_point_sources();
//    print();
}

void efield::set(std::string fieldName, const char *arg) {
    if (fieldName == "efield")
        ;
    else if (fieldName == "setSize") {
        // arg is three integers, comma separated
        int nRows;
        int nCols;
        int nLambda;
        sscanf(arg, "%d, %d, %d", &nRows, &nCols, &nLambda);
        set_size(nRows, nCols, nLambda);
    }
    else if (fieldName == "filename") {
        // arg is two filenames separated by a comma, each giving
        // the .fits amplitude and phase filename that contains the complexMask definition
        const char *cPtr = strchr(arg, ',');
        int strAmpLen = cPtr - arg;
        int strPhLen = strlen(arg) - strAmpLen + 1;
        std::cout << "strAmpLen " << strAmpLen << ", strPhLen " << strPhLen << std::endl;
        char *strAmp = new char[strAmpLen + 1];
        char *strPh = new char[strPhLen + 1];
        strncpy(strAmp, arg, strAmpLen);
        strAmp[strAmpLen] = '\0';
        strncpy(strPh, cPtr + 1, strPhLen);
        strPh[strPhLen] = '\0';
        std::cout << "loading " << strAmp << " and " << strPh << std::endl;
        init(strAmp, strPh);
        //        draw("initial FPM");
    }
    else if (fieldName == "constantValue") {
        // arg is one double value
        init(atof(arg));
    }
    else if (fieldName == "pointSource") {
        // arg is three doubles, comma separated
        double flux;
        double tipX;
        double tiltY;
        sscanf(arg, "%lf, %lf, %lf", &flux, &tipX, &tiltY);
        tipX = tipX*globalTelescope->get("magnification");
        tiltY = tiltY*globalTelescope->get("magnification");
        add_point_source(sqrt(flux), tipX, tiltY);
    }
    else if (fieldName == "referenceWavelength") {
        // arg is one double value
        referenceWavelength = atof(arg);
    }
    else if (fieldName == "wavelength") {
        // arg is one double value
        lambdaData = new lambdaDataClass(atof(arg));
    }
    else if (fieldName == "wavelengthRange") {
        // args are start wavelength, end wavelength, delta wavelength
        set_wavelength_range(arg);
    }
    else if (fieldName == "wavelengthBandwidth") {
        // args are central wavelength, bandwidth fraction, n wavelengths
        set_wavelength_bandwidth(arg);
    }
    else if (fieldName == "wavelengthList") {
        // args are a space-separated list of wavelengths of arbitrary length
        set_wavelength_list(arg);
    }
    else if (fieldName == "outputDirectory") {
        // arg is a string
        outputDirectory = new char[strlen(arg)+1];
        strcpy(outputDirectory, arg);
    }
    else if (fieldName == "name") {
        // arg is a string
        name = new char[strlen(arg)+1];
        strcpy(name, arg);
    }
    else if (fieldName == "beamRadiusPhysical") {
        // arg is one double value
        beamRadiusPhysical = atof(arg);
    }
    else if (fieldName == "pixelScale") {
        // arg is one double value
        pixelScale = atof(arg);
    }
    else if (fieldName == "arraySizePixels") {
        // arg is one integer value
        arraySizePixels = atoi(arg);
    }
    else
        std::cout << "!!!!! efield set: unknown fieldName: " << fieldName << std::endl;
}

void efield::set_wavelength_range(const char *arg) {
    // args are start wavelength, end wavelength, delta wavelength
    std::cout << "set_wavelength_range arg: " << arg << std::endl;
    sscanf(arg, "%lf, %lf, %lf", &wavelengthStart, &wavelengthEnd, &deltaWavelength);
    std::cout << "read wavelengthRange " << wavelengthStart << ", " << wavelengthEnd << ", " << deltaWavelength << std::endl;
    int nWavelengths = round((wavelengthEnd - wavelengthStart)/deltaWavelength) + 1;
    std::cout << "nWavelengths = " << nWavelengths << std::endl;
    lambdaData = new lambdaDataClass[nWavelengths];
    std::cout << "wavelength list:" << std::endl;
    for (int i=0; i<nWavelengths; i++) {
        lambdaData[i].set_wavelength(wavelengthStart + i*deltaWavelength);
        std::cout << lambdaData[i].get_wavelength() << ", ";
    }
    std::cout << std::endl;
    set_size(E[0][0]->n_rows, E[0][0]->n_cols, nWavelengths);
}

void efield::set_wavelength_bandwidth(const char *arg) {
    // args are central wavelength, bandwidth fraction, n wavelengths
    double centralWavelength;
    double bandwidthFraction;
    int nWavelengths;
    
    std::cout << "set_wavelength_bandwidth arg: " << arg << std::endl;
    sscanf(arg, "%lf, %lf, %d", &centralWavelength, &bandwidthFraction, &nWavelengths);
    std::cout << "read wavelengthBandwidth " << centralWavelength << ", " << bandwidthFraction << ", " << nWavelengths << std::endl;
    wavelengthStart = centralWavelength - centralWavelength*bandwidthFraction/2;
    wavelengthEnd = centralWavelength + centralWavelength*bandwidthFraction/2;
    deltaWavelength = (wavelengthEnd - wavelengthStart)/((double) nWavelengths - 1);
    std::cout << "wavelengthStart = " << wavelengthStart << ", wavelengthEnd = " << wavelengthEnd << ", deltaWavelength = " << deltaWavelength << std::endl;
    lambdaData = new lambdaDataClass[nWavelengths];
    std::cout << "wavelength list:" << std::endl;
    for (int i=0; i<nWavelengths; i++) {
        lambdaData[i].set_wavelength(wavelengthStart + i*deltaWavelength);
        std::cout << lambdaData[i].get_wavelength() << ", ";
    }
    std::cout << std::endl;
    set_size(E[0][0]->n_rows, E[0][0]->n_cols, nWavelengths);
}

void efield::set_wavelength_list(const char *arg) {
    // arg is a variable length listof wavelengths
    std::cout << "set_wavelength_list arg: " << arg << std::endl;
    int retVal = 1;
    int wCount = 0;
    int cPos = 0;
    double waveLengths[MAXNWAVELENGTHS];
    while (retVal != EOF) {
        int ncPos;
        retVal = sscanf(arg+cPos, "%lf %n", &waveLengths[wCount], &ncPos);
        cPos += ncPos;
        std::cout << "wavelength " << wCount << ": " << waveLengths[wCount] << ", retVal = " << retVal << ", cPos = " << cPos << std::endl;
        wCount++;
        if (wCount >= MAXNWAVELENGTHS)
            assert(NULL);
    }
    wCount--;
    std::cout << wCount << " wavelengths read" << std::endl;
    lambdaData = new lambdaDataClass[wCount];
    std::cout << "wavelength list:" << std::endl;
    for (int i=0; i<wCount; i++) {
        lambdaData[i].set_wavelength(waveLengths[i]);
        std::cout << lambdaData[i].get_wavelength() << ", ";
    }
    std::cout << std::endl;
    set_size(E[0][0]->n_rows, E[0][0]->n_cols, wCount);
}

void efield::set_optical_parameters(void) {
    assert(beamRadiusPhysical != 0 & pixelScale != 0 & arraySizePixels != 0);
    
    std::cout << "set_optical_parameters: beamRadiusPhysical = " << beamRadiusPhysical << std::endl;
    std::cout << "set_optical_parameters: pixelScale = " << pixelScale << std::endl;
    beamSizePixels = ceil(2*beamRadiusPhysical/pixelScale); // mablab: NBeam
    std::cout << "set_optical_parameters: beamSizePixels = " << beamSizePixels << std::endl;
    std::cout << "set_optical_parameters: arraySizePixels = " << arraySizePixels << std::endl;
    arrayRadiusPhysical = beamRadiusPhysical*arraySizePixels/beamSizePixels; // matlab: prad
    std::cout << "set_optical_parameters: arrayRadiusPhysical = " << arrayRadiusPhysical << std::endl;

}

void efield::set_array_geometry(void) {
    std::cout << arrayRadiusPhysical << E[0][0]->n_rows << E[0][0]->n_cols << std::endl;
    assert(arrayRadiusPhysical != 0 & E[0][0]->n_rows != 0 & E[0][0]->n_cols != 0);
    
    arrayGeometry.set_geometry(E[0][0], pixelScale);
    initArrayGeometry.set_geometry(E[0][0], pixelScale);

//    arrayGeometry.print("Efield");
}

void efield::set_size(int nRows, int nColumns, int nLambda) {
//    if (E[0][0]->n_slices != nLambda)
//        lambdaData = new lambdaDataClass[nLambda];
    std::cout << "set_size: nRows = " << nRows << ", nColumns = " << nColumns << ", nLambda = " << nLambda << std::endl;
    for (int s=0; s<E.size(); s++) {
        for (int p=0; p<E[s].size(); p++) {
            E[s][p]->set_size(nRows, nColumns, nLambda);
        }
    }
    arraySizePixels = nRows;
}

inline std::complex<double> efield::getV(int r, int c, int lambda, int p, int s) {
    return (*E[s][p])(r, c, lambda);
}

inline void efield::setV(std::complex<double> val, int r, int c, int lambda, int p, int s) {
    (*E[s][p])(r, c, lambda) = val;
}

void efield::save(const char *coreName) {
    for (int s=0; s<E.size(); s++) {
        for (int p=0; p<E[s].size(); p++) {
            
            arma::cube reE = real(*E[s][p]);
            arma::cube imE = imag(*E[s][p]);
            
            // fits is transposed so...
            for (int s=0; s<reE.n_slices; s++) {
                inplace_strans(reE.slice(s));
                inplace_strans(imE.slice(s));
            }
            
            int nDims = 3;
            long dims[3];
            int dataType = TDOUBLE;
            
            dims[0] = E[s][p]->n_rows;
            dims[1] = E[s][p]->n_cols;
            dims[2] = E[s][p]->n_slices;
            
            std::string fname = "!" + (std::string)outputDirectory + "E_s" + std::to_string(s) + "_p" + std::to_string(p) + "_" + (std::string)coreName + "_re.fits";
            write_fits_array(fname.c_str(), (double(*)[3])&reE(0,0,0), nDims, dims, dataType);
            
            fname = "!" + (std::string)outputDirectory + "E_s" + std::to_string(s) + "_p" + std::to_string(p) + "_" + (std::string)coreName + "_im.fits";
            write_fits_array(fname.c_str(), (double(*)[3])&imE(0,0,0), nDims, dims, dataType);

            arma::cube absE = abs(*E[s][p]);
            arma::cube argE = arg(*E[s][p]);
            
            // fits is transposed so...
            for (int s=0; s<absE.n_slices; s++) {
                inplace_strans(absE.slice(s));
                inplace_strans(argE.slice(s));
            }
            
            fname = "!" + (std::string)outputDirectory + "E_s" + std::to_string(s) + "_p" + std::to_string(p) + "_" + (std::string)coreName + "_amp.fits";
            write_fits_array(fname.c_str(), (double(*)[3])&absE(0,0,0), nDims, dims, dataType);
            
            fname = "!" + (std::string)outputDirectory + "E_s" + std::to_string(s) + "_p" + std::to_string(p) + "_" + (std::string)coreName + "_phase.fits";
            write_fits_array(fname.c_str(), (double(*)[3])&argE(0,0,0), nDims, dims, dataType);

        }
    }
}

void efield::draw(const char *title, const char *drawType) {
    for (int s=0; s<E.size(); s++) {
        for (int p=0; p<E[s].size(); p++) {
            arma::cube drawE;
            
            if (!strcmp(drawType, "amp"))
                drawE = abs(*E[s][p]);
            else if (!strcmp(drawType, "logamp"))
                drawE = log10(abs(*E[s][p]));
            else if (!strcmp(drawType, "phase"))
                drawE = arg(*E[s][p]);
            else
                std::cout << "!!!!! efield draw: unknown type: " << drawType << std::endl;
            
            for (int k = 0; k < drawE.n_slices; k++) {
                char str[200];
                sprintf(str, "%s w%d p%d %s", title, k, p, drawType);
                draw_mat(drawE.slice(k), str, "matlab");
            }
        }
    }
}

void efield::print(const char *header) {
    std::cout << "I'm an efield: " << header << std::endl;
    std::cout << "num sources = " << E.size() << std::endl;
    std::cout << "num polarizations = " << E[0].size() << std::endl;
    std::cout << "size for each source and polarization = " << size(*E[0][0]) << std::endl;
    std::cout << "beamRadiusPhysical = " << beamRadiusPhysical << std::endl;
    std::cout << "pixelScale = " << pixelScale << std::endl;
    std::cout << "arraySizePixels = " << arraySizePixels << std::endl;
    std::cout << "beamSizePixels = " << beamSizePixels << std::endl;
    std::cout << "arrayRadiusPhysical = " << arrayRadiusPhysical << std::endl;
    std::cout << "referenceWavelength = " << referenceWavelength << std::endl;
    std::cout << "wavelengths:" << std::endl;
    for (int i=0; i<E[0][0]->n_slices; i++) {
        std::cout << lambdaData[i].get_wavelength() << ", ";
    }
    std::cout << std::endl;
    std::cout << pointSourceList.size() << " point sources:" << std::endl;
    for (int i=0; i<pointSourceList.size(); i++) {
        pointSourceList[i].print();
    }
}

