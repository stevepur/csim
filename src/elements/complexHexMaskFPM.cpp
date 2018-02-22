//
//  complexHexMaskFPM.cpp
//  csim
//
//  Created by steve on 5/24/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include "complexHexMaskFPM.hpp"
#include <iostream>
#include <string>
#include <cstring>
#include <assert.h>
#include "../lib/csim_lib.hpp"
#include "../coronagraph/coronagraph.hpp"

complexHexMaskFPM::complexHexMaskFPM() {
}

complexHexMaskFPM::complexHexMaskFPM(initCommandSet*& cmdBlock, double z) {
    zoomFactor = z;
    init(cmdBlock);
}

void complexHexMaskFPM::execute(efield* E, arma::cx_mat& tcMat, int sl, double time) {
    std::complex<double> i1(0, 1);
    
    // make the interpolated complex mask
    arma::cx_mat complexIntMask = make_complex_intpolated_mask(E->lambdaData[sl].lambda, time);
    
    // apply the complex mask
    int row0 = tcMat.n_rows/2 - complexIntMask.n_rows/2;
    int col0 = tcMat.n_cols/2 - complexIntMask.n_cols/2;
#pragma omp parallel
    {
#pragma omp for
        for (int c=0; c<complexIntMask.n_cols; ++c) {
            for (int r=0; r<complexIntMask.n_rows; ++r)
                tcMat(row0 + r, col0 + c) *= complexIntMask(r,c);
        }
    }
}

arma::cx_mat complexHexMaskFPM::make_complex_intpolated_mask(double lambda, double time) {
    std::complex<double> i1(0, 1);

    // set the sag values
#pragma omp parallel
    {
#pragma omp for
        for (int r=0; r<interpHexNum.n_rows; r++)
            for (int c=0; c<interpHexNum.n_cols; c++)
                for (int s=0; s<interpHexNum.n_slices; s++) {
                    if (interpHexNum(r,c,s) > -1)
                        interpSagVals(r,c,s) = fpmSags((int)interpHexNum(r,c,s));
                    else
                        interpSagVals(r,c,s) = 0.0;
                }
    }
    
    // make the interpolated complex mask
    arma::cx_mat complexIntMask = sum(interpNSubPix % exp(-2*M_PI*i1*2.0*interpSagVals/lambda), 2)/(nSubPix*nSubPix);
//    save_mat("used_complexIntMask.fits", complexIntMask);
    
    return complexIntMask;
}

void complexHexMaskFPM::init(initCommandSet*& cmdBlock) {
    std::cout << "initing a complexHexMaskFPM" << std::endl;
    
    for (int c=0; c<cmdBlock->commandList.size(); c++) {
        set(cmdBlock->commandList[c]->getCmdStr(),
            cmdBlock->commandList[c]->getArgStr());
    }
    init_mask();
}

void complexHexMaskFPM::set(std::string fieldName, const char *arg) {
    bool found = celem::set(fieldName, arg);
    if (fieldName == "complexHexMaskFPM")
        ;
    else if (fieldName == "maskNRings") {
        // arg is a single integer
        maskNRings = atoi(arg);
    } else if (fieldName == "fpmArraySize") {
        // arg is a single integer
        fpmArraySize = atoi(arg);
    } else if (fieldName == "nSubPix") {
        // arg is a single integer
        nSubPix = atoi(arg);
    } else if (fieldName == "hexStep") {
        // arg is a single double
        hexStep = atof(arg);
    } else if (fieldName == "fpmScaleFactor") {
        // arg is a single double
        fpmScaleFactor = atof(arg);
    } else if (fieldName == "hexGap") {
        // arg is a single double
        hexGap = atof(arg);
    } else if (fieldName == "lambdaRef") {
        // arg is a single double
        lambdaRef = atof(arg);
    } else if (fieldName == "fpmFRatio") {
        // arg is a single double
        fpmFRatio = atof(arg);
    } else if (fieldName == "fpmRadius") {
        // arg is a single double
        fpmRadius = atof(arg);
    } else if (fieldName == "sagFilename") {
        // arg is one filename
        load_sags(arg);
    } else if (fieldName == "maskFilenameRoot") {
        // arg is one filename root, which will get additions for the mask data files
        maskFilenameRoot = arg;
    } else if (!found)
        std::cout << "!!! complexHexMaskFPM bad set field name: " << fieldName << std::endl;
}

void complexHexMaskFPM::load_sags(const char *filename){
    load_vec(filename, fpmSags);
}

void complexHexMaskFPM::init_mask(void) {
    int fileStatus = 0;
    
    std::string s1 = maskFilenameRoot+"_nSubPix.fits";
    fileStatus += load_cube(s1.c_str(), interpNSubPix);
    std::string s2 = maskFilenameRoot+"_hexNum.fits";
    fileStatus += load_cube(s2.c_str(), interpHexNum);
    
    if (fileStatus) {
        std::cout << "did not find mask interpolation data files, creating" << std::endl;
        compute_hex_centers();
        make_hex_array();
        make_interpolation_data();
        arma::mat testSagArray = sum(interpNSubPix % interpSagVals, 2)/(nSubPix*nSubPix);
        
        save_cube(s1.c_str(), interpNSubPix);
        save_cube(s2.c_str(), interpHexNum);
    } else {
        fpmScale = compute_fpmScale(interpNSubPix.n_rows);
        std::cout << "fpmScale = " << fpmScale << std::endl;
        interpSagVals = arma::zeros<arma::cube>(interpNSubPix.n_rows, interpNSubPix.n_cols, interpNSubPix.n_slices);
    }
}

void complexHexMaskFPM::compute_hex_centers(void) {
    int ii1Max = (int) (maskNRings/3 + 2);
    int jj1Max = (int) (maskNRings/sqrt(3) + 2);
    
    int hIndex = 0;
    double badHexNum = 1e8;
    
    hexX = arma::zeros<arma::vec>(2000);
    hexY = arma::zeros<arma::vec>(2000);
    hexRing.set_size(2000);
    hexNum.set_size(2000);
    
    for (int ii1 = -ii1Max; ii1 < ii1Max; ii1++) {
        for (int jj1 = -jj1Max; jj1 < jj1Max; jj1++) {
            // place hex center on even row
            double hx = hexStep*ii1*3.0;
            double hy = hexStep*jj1*sqrt(3.0);
            int ring = (int) sqrt(hx*hx + hy*hy);
            hexX(hIndex) = hx;
            hexY(hIndex) = hy;
            hexRing(hIndex) = ring;
            if (ring < maskNRings)
                hexNum(hIndex) = 0;
            else
                hexNum(hIndex) = badHexNum + hIndex;
            hIndex = hIndex + 1;
            
            // place hex center on odd row
            hx += hexStep*1.5;
            hy += hexStep*sqrt(3.0)/2.0;
            ring = (int) sqrt(hx*hx + hy*hy);
            hexX(hIndex) = hx;
            hexY(hIndex) = hy;
            hexRing(hIndex) = ring;
            if (ring < maskNRings)
                hexNum(hIndex) = 0;
            else
                hexNum(hIndex) = badHexNum + hIndex;
            hIndex = hIndex + 1;
        }
    }
    hexNum.resize(hIndex);
    hexX.resize(hIndex);
    hexY.resize(hIndex);
    hexRing.resize(hIndex);
    std::cout << "# of hexes: " << hexNum.n_elem << std::endl;
    
    
    // now organize the hex list by ring, assigning each hex a number
    int hCount = 0;
    for (int r=0; r<=maskNRings; r++) {
        arma::uvec inRing = (hexRing == r) % (hexNum == 0);
        arma::uvec inRingIdx = arma::find(inRing);
        for (int rri = 0; rri < inRingIdx.n_elem; rri++) {
            hexNum(inRingIdx(rri)) = hCount;
            hCount = hCount + 1;
        }
    }
    
    // sort in order of increasing hexNum
    arma::uvec sortIdx = sort_index(hexNum);
    hexNum = hexNum(sortIdx);
    hexX = hexX(sortIdx);
    hexY = hexY(sortIdx);
    hexRing = hexRing(sortIdx);
    
    // trim to the legal hexes
    arma::uvec goodRingIdx = arma::find(hexNum < badHexNum);
    hexNum.resize(goodRingIdx.n_elem);
    hexX.resize(goodRingIdx.n_elem);
    hexY.resize(goodRingIdx.n_elem);
    hexRing.resize(goodRingIdx.n_elem);
    
    
    FILE *outFile = fopen("hex_positions.txt", "w");
    for (int i=0; i<hexNum.n_elem; i++)
        fprintf(outFile, "%d %f %f %d\n", hexNum(i), hexX(i), hexY(i), hexRing(i));
    fclose(outFile);
    
}

void complexHexMaskFPM::make_hex_array(void) {
    
    // for vectors vp = (sqrt(3)/2, 0.5), vm = (sqrt(3)/2, -0.5)
    double vx = sqrt(3)/2;
    double vy = 0.5;
    
    double hexStepPix = (0.5*fpmArraySize)/maskNRings * fpmScaleFactor;
    arma::vec hexXFpmIndexCoords = 0.5*fpmArraySize + hexX * hexStepPix;
    arma::vec hexYFpmIndexCoords = 0.5*fpmArraySize + hexY * hexStepPix;
    double hexRadius = hexStepPix*(1.0-hexGap)*(sqrt(3.0)/2.0); // inner radius
    double hexOuterRadius = 2*hexRadius/sqrt(3.0);
    double closeHexSq = hexOuterRadius*hexOuterRadius;
    
    fpmDesignArray = -1.0*arma::ones<arma::mat>(fpmArraySize, fpmArraySize);
    fpmDesignSagArray = arma::zeros<arma::mat>(fpmArraySize, fpmArraySize);
    for (int h=0; h<hexNum.n_elem; h++) {
        for (int ii=0; ii<fpmArraySize; ii++) {
            for (int jj=0; jj<fpmArraySize; jj++) {
                double dx = (double)jj - hexXFpmIndexCoords(h);
                double dy = (double)ii - hexYFpmIndexCoords(h);
                double r2 = dx*dx + dy*dy;
                if (r2 > closeHexSq)
                    continue;
                
                // compare the vector dv=(dx, dy) from this point to the
                // vp and vm by taking the absolute magnitude of the dot
                // products dv.vp and dv.vm
                // demand that abs(dv.vp) and abs(dv.vm) < hexRadius
                // and dy < hexRadius
                double dpx = dx*vx;  // help out the compiler optimizer
                double dpy = dy*vy;
                if (fabs(dy) < hexRadius
                    & fabs(dpx + dpy) < hexRadius
                    & fabs(dpx - dpy) < hexRadius) {
                    fpmDesignArray(ii, jj) = (double) hexNum(h);
                    fpmDesignSagArray(ii, jj) = fpmSags(h);
                }
            }
        }
    }
    
    FILE *outFile = fopen("hex_positions_fpmIndex.txt", "w");
    for (int i=0; i<hexNum.n_elem; i++)
        fprintf(outFile, "%d %f %f %d\n", hexNum(i), hexXFpmIndexCoords(i), hexYFpmIndexCoords(i), hexRing(i));
    fclose(outFile);
}

double complexHexMaskFPM::compute_fpmScale(int fpmSize) {
    return 2*(initialEfield->beamRadiusPhysical/(initialEfield->pixelScale*fpmSize*zoomFactor))*lambdaRef*fpmFRatio;
}

void complexHexMaskFPM::make_interpolation_data(void) {
    int fpmSize = initialEfield->E[0][0]->n_rows;
    arma::ivec hexList;
    
    // compute the pixel scale for the fpm array
    fpmScale = compute_fpmScale(fpmSize);
    std::cout << "fpmScale = " << fpmScale << std::endl;
    
    // precompute some things to help out the compiler optimizer
    double ac1 = -fpmSize/2.0 - 0.5;
    double xc1 = (0.5/fpmRadius)*fpmScaleFactor;
    
    interpNSubPix = arma::zeros<arma::cube>(fpmSize, fpmSize, 3);
    interpSagVals = arma::zeros<arma::cube>(fpmSize, fpmSize, 3);
//    interpHexNum.set_size(fpmSize, fpmSize, 3);
    interpHexNum =  -1.0*arma::ones<arma::cube>(fpmSize, fpmSize, 3);
    
    hexList.set_size(3);
    for (int i=0; i<fpmSize; i++) {
        std::cout << i << " ";
        fflush(stdout);
        // check midpoint of the cell to find if it's on the design array
        int iii = (int) (nSubPix/2.0);
        // compute the physical coordinate of the point
        double x = (i + ac1 + (0.5 + iii)/((double) nSubPix)) * fpmScale;
        // compute the index of this physical point on the design array
        int ii1 = (int)((0.5 + x*xc1)*fpmArraySize + 0.5);
        
        for (int j=0; j<fpmSize; j++) {
            // check midpoint of the cell to find if it's on the design array
            int jjj = (int) (nSubPix/2.0);
            // compute the physical coordinate of the point
            double y = (j + ac1 + (0.5 + jjj)/((double) nSubPix)) * fpmScale;
            // compute the index of this physical point on the design array
            int jj1 = (int)((0.5 + y*xc1)*fpmArraySize + 0.5);
            
            if (ii1 < 0 | ii1 >= fpmArraySize | jj1 < 0 | jj1 >= fpmArraySize)
                continue;
            
            // now collect sub-sample data
            hexList.reset();
            for (int iii = 0; iii < nSubPix; iii++) {
                for (int jjj = 0; jjj < nSubPix; jjj++) {
                    // physical coordinates of this subpixel
                    x = (i + ac1 + (0.5 + iii)/((double) nSubPix)) * fpmScale;
                    y = (j + ac1 + (0.5 + jjj)/((double) nSubPix)) * fpmScale;
                    // coordinates of this subpixel on the design array
                    ii1 = (int)((0.5 + x*xc1)*fpmArraySize + 0.5);
                    jj1 = (int)((0.5 + y*xc1)*fpmArraySize + 0.5);
                    if (ii1 < 0 | ii1 >= fpmArraySize | jj1 < 0 | jj1 >= fpmArraySize) {
                        interpNSubPix(i,j,0)++;
                        continue;
                    }
                    
                    // get the hexnum for this subcell
                    int inHexNum = (int) fpmDesignArray(ii1, jj1);
                    if (inHexNum < 0) {
                        interpNSubPix(i,j,0)++;
                        continue;
                    }
                    // is this hexNum in this cell's hexList?
                    arma::uvec hlii = find(hexList == inHexNum);
                    if (hlii.is_empty()) {
                        // no, so add it to the list of hex numbers in this cell
                        int curSize = hexList.n_elem;
                        hexList.resize(curSize + 1);
                        hexList(curSize) = inHexNum;
                        assert(hexList.n_elem <= 3);
                    }
                    
                    // now find the slot with this hexNum
                    arma::uvec hexListIdx = find(hexList == inHexNum);
                    assert(hexListIdx.n_elem == 1);
                    // and set the interpolatin data for this cell
                    if (!hexListIdx.is_empty()) {
                        interpNSubPix(i,j,hexListIdx(0))++;
//                        interpSagVals(i,j,hexListIdx(0)) = fpmSags(inHexNum);
                        interpHexNum(i,j,hexListIdx(0)) = (double) inHexNum;
                    }
                } // jjj loop
            } // iii loop
        } // j loop
    } // i loop
    std::cout << std::endl;
    
    // fix zeros in interpNSubPix
    arma::mat ss = sum(interpNSubPix, 2);
    for (int i=0; i<ss.n_rows; i++)
        for (int j=0; j<ss.n_cols; j++)
            if (ss(i,j) == 0)
                interpNSubPix(i, j, 0) = nSubPix*nSubPix;
    
}

void complexHexMaskFPM::print(const char *hdr) {
    std::cout << "complexHexMaskFPM " << hdr << std::endl;
    std::cout << "maskNRings = " << maskNRings << std::endl;
    std::cout << "fpmArraySize = " << fpmArraySize << std::endl;
    std::cout << "nSubPix = " << nSubPix << std::endl;
    std::cout << "hexStep = " << hexStep << std::endl;
    std::cout << "fpmScaleFactor = " << fpmScaleFactor << std::endl;
    std::cout << "hexGap = " << hexGap << std::endl;
    std::cout << "lambdaRef = " << lambdaRef << std::endl;
    std::cout << "fpmFRatio = " << fpmFRatio << std::endl;
    std::cout << "fpmRadius = " << fpmRadius << std::endl;
}

void complexHexMaskFPM::draw(const char *title) {
    char str[200];
    sprintf(str, "%s complexHexMaskFPM %s", title, name);
}

