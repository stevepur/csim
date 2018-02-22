//
//  csim_fits.cpp
//  csim
//
//  Created by steve on 4/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include <iostream>
#include <assert.h>

#include "csim_geom.hpp"
#include "csim_plot.hpp"

// set up array geometry including pixel radius and angle
void arrayGeom::set_geometry(int nRows, int nCols, double physicalRadius, bool display) {
    
    set_xy(nRows, nCols, physicalRadius, display);
    
    set_mesh();
    
    if (display) {
        print();
        draw();
    }
}

// set up just the x and y coordinates of each pixel
void arrayGeom::set_xy(int nRows, int nCols, double physicalRadius, bool display) {
    
    physicalSize = 2*physicalRadius;
    physicalSizeX = physicalSize;
    physicalSizeY = physicalSize;
    pixelSizeX = physicalSizeX/nCols;
    pixelSizeY = physicalSizeY/nRows;
    
    // matlab: plane.x = ((1:plane.N) - 1 - plane.N/2)/(plane.N)*plane.D;
    pixelX = arma::zeros<arma::vec>(nCols);
    for (int i=0; i<nCols; i++)
        pixelX[i] = ((double) i - ((double)nCols)/2.)*pixelSizeX;
    pixelY = arma::zeros<arma::vec>(nRows);
    for (int i=0; i<nRows; i++)
        pixelY[i] = ((double) i - ((double)nRows)/2.)*pixelSizeY;
    
    if (display) {
        print();
        draw();
    }
}

void arrayGeom::set_xy_m1(int nRows, int nCols, double physicalRadius, bool display) {
    
    physicalSize = 2*physicalRadius;
    physicalSizeX = physicalSize;
    physicalSizeY = physicalSize;
    pixelSizeX = physicalSizeX/nCols;
    pixelSizeY = physicalSizeY/nRows;
    
    // sci.x = -(sci.gridsize - sci.dx)/2 : sci.dx : (sci.gridsize - sci.dx)/2;
    pixelX = arma::zeros<arma::vec>(nCols);
    for (int i=0; i<nCols; i++)
        pixelX[i] = -(physicalSizeX - pixelSizeX)/2.0 + i*pixelSizeX;
    pixelY = arma::zeros<arma::vec>(nRows);
    for (int i=0; i<nRows; i++)
        pixelY[i] = -(physicalSizeY - pixelSizeY)/2.0 + i*pixelSizeY;

    if (display) {
        print();
        draw();
    }
}

// set up array geometry including pixel radius and angle
void arrayGeom::set_geometry(double flD, double samplesPerFld, double fovInFld, bool display) {
    
    // matlab: sci.N = ceil(FOVflD * samplesPerflD);
    int N = ceil(fovInFld*samplesPerFld);
    
    // matlab: sci.dx = flD / samplesPerflD;
    pixelSizeX = flD/samplesPerFld;
    pixelSizeY = flD/samplesPerFld;
    
    // matlab: sci.gridsize = sci.N*sci.dx;
    physicalSize = N*pixelSizeX;
    physicalSizeX = physicalSize;
    physicalSizeY = physicalSize;
    
    // sci.x = -(sci.gridsize - sci.dx)/2 : sci.dx : (sci.gridsize - sci.dx)/2;
    pixelX = arma::zeros<arma::vec>(N);
    for (int i=0; i<N; i++)
        pixelX[i] = -(physicalSizeX - pixelSizeX)/2.0 + i*pixelSizeX;
    pixelY = arma::zeros<arma::vec>(N);
    for (int i=0; i<N; i++)
        pixelY[i] = -(physicalSizeY - pixelSizeY)/2.0 + i*pixelSizeY;
    
    set_mesh();
    
    if (display) {
        print();
        draw();
    }
}

void arrayGeom::set_mesh(void) {
    
    // matlab: [plane.xx plane.yy] = meshgrid(plane.x, plane.y);
    pixelXX = arma::repmat(pixelX.t(), pixelX.n_elem, 1);
    pixelYY = arma::repmat(pixelY, 1, pixelY.n_elem);
    
    // matlab: plane.rr = sqrt(plane.xx.^2 + plane.yy.^2);
    pixelRR = sqrt(pow(pixelXX, 2) + pow(pixelYY, 2));
    // matlab: plane.ttheta = atan2(plane.yy, plane.xx);
    pixelTT = atan2(pixelYY, pixelXX);
    
}

void arrayGeom::print(const char *hdr) {
    std::cout << "arrayGeom " << hdr << std::endl;
    printf("physicalSize = %0.25f\n", physicalSize);
    printf("physicalSizeX = %0.25f\n", physicalSizeX);
    printf("physicalSizeY = %0.25f\n", physicalSizeY);
    printf("pixelSizeX = %0.25f\n", pixelSizeX);
    printf("pixelSizeY = %0.25f\n", pixelSizeY);
    std::cout << "N elements X = " << pixelX.n_elem << std::endl;
    std::cout << "N elements Y = " << pixelY.n_elem << std::endl;
    for (int i=0; i<5; i++)
        printf("pixel %d: X = %0.23f, Y = %0.23f\n", i, pixelX[i], pixelY[i]);
    for (int i=pixelX.n_elem - 5; i<pixelX.n_elem; i++)
        printf("pixel %d: X = %0.23f, Y = %0.23f\n", i, pixelX[i], pixelY[i]);
}

void arrayGeom::draw(void) {
    draw_mat(pixelXX, pixelX(0), pixelX(pixelX.n_elem-1), pixelY(0), pixelY(pixelY.n_elem-1), "arrayGeom pixelXX", "matlab");
    draw_mat(pixelYY, pixelX(0), pixelX(pixelX.n_elem-1), pixelY(0), pixelY(pixelY.n_elem-1), "arrayGeom pixelYY", "matlab");
    draw_mat(pixelRR, pixelX(0), pixelX(pixelX.n_elem-1), pixelY(0), pixelY(pixelY.n_elem-1), "arrayGeom pixelRR", "matlab");
    draw_mat(pixelTT, pixelX(0), pixelX(pixelX.n_elem-1), pixelY(0), pixelY(pixelY.n_elem-1), "arrayGeom pixelTT", "matlab");
}




