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
void arrayGeom::set_geometry(arma::cx_cube *E, double pixelScale, double origin, bool display) {
    set_geometry(E->n_cols, E->n_rows, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_geometry(arma::cx_mat *E, double pixelScale, double origin, bool display) {
    set_geometry(E->n_cols, E->n_rows, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_geometry(arma::cube *E, double pixelScale, double origin, bool display) {
    set_geometry(E->n_cols, E->n_rows, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_geometry(arma::mat *E, double pixelScale, double origin, bool display) {
    set_geometry(E->n_cols, E->n_rows, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_xy(arma::cx_cube *E, double pixelScale, double origin, bool display) {
    set_xy(E->n_cols, E->n_rows, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_xy(arma::cx_mat *E, double pixelScale, double origin, bool display) {
    set_xy(E->n_cols, E->n_rows, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_xy(arma::cube *E, double pixelScale, double origin, bool display) {
    set_xy(E->n_cols, E->n_rows, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_xy(arma::mat *E, double pixelScale, double origin, bool display) {
    set_xy(E->n_cols, E->n_rows, pixelScale, pixelScale, origin, origin, display);
}

//

void arrayGeom::set_geometry(arma::cx_cube& E, double pixelScale, double origin, bool display) {
    set_geometry(E.n_cols, E.n_rows, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_geometry(arma::cx_mat& E, double pixelScale, double origin, bool display) {
    set_geometry(E.n_cols, E.n_rows, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_geometry(arma::cube& E, double pixelScale, double origin, bool display) {
    set_geometry(E.n_cols, E.n_rows, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_geometry(arma::mat& E, double pixelScale, double origin, bool display) {
    set_geometry(E.n_cols, E.n_rows, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_xy(arma::cx_cube& E, double pixelScale, double origin, bool display) {
    set_xy(E.n_cols, E.n_rows, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_xy(arma::cx_mat& E, double pixelScale, double origin, bool display) {
    set_xy(E.n_cols, E.n_rows, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_xy(arma::cube& E, double pixelScale, double origin, bool display) {
    set_xy(E.n_cols, E.n_rows, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_xy(arma::mat& E, double pixelScale, double origin, bool display) {
    set_xy(E.n_cols, E.n_rows, pixelScale, pixelScale, origin, origin, display);
}

//

void arrayGeom::set_geometry(int N, double pixelScale, double origin, bool display) {
    set_geometry(N, N, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_geometry(int nX, int nY, double pixelScale, double origin, bool display) {
    set_geometry(nX, nY, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_geometry(int nX, int nY, double pixelScaleX, double pixelScaleY, double oX, double oY, bool display) {

    set_xy(nX, nY, pixelScaleX, pixelScaleY, oX, oY, display);
    
    set_mesh();
    
    if (display) {
        print();
        draw();
    }
}

// set up just the x and y coordinates of each pixel
void arrayGeom::set_xy(int N, double pixelScale, double origin, bool display) {
    set_xy(N, N, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_xy(int nX, int nY, double pixelScale, double origin, bool display) {
    set_xy(nX, nY, pixelScale, pixelScale, origin, origin, display);
}

void arrayGeom::set_xy(int nX, int nY, double pixelScaleX, double pixelScaleY, double oX, double oY, bool display) {

    physicalSizeX = pixelScaleX*nX;
    physicalSizeY = pixelScaleY*nY;
    pixelSizeX = pixelScaleX;
    pixelSizeY = pixelScaleY;
    originX = oX;
    originY = oY;

    physicalSize = physicalSizeX;

    // sci.x = -(sci.gridsize - sci.dx)/2 : sci.dx : (sci.gridsize - sci.dx)/2;
    // this is symmetric:
    // i=0 => pixelX = -(physicalSizeX - pixelSizeX)/2.0;
    // i=nCols-1 => pixelX[ = -(physicalSizeX - pixelSizeX)/2.0 + (nCols - 1)*pixelSizeX
    //          = -physicalSizeX/2.0 + pixelSizeX/2.0 + nCols*pixelSizeX - pixelSizeX
    //          = -physicalSizeX/2.0 + pixelSizeX/2.0 + physicalSizeX - pixelSizeX
    //          = (physicalSizeX - pixelSizeX)/2.0
    pixelX = arma::zeros<arma::vec>(nX);
    for (int i=0; i<nX; i++)
        pixelX[i] = -(physicalSizeX - pixelSizeX)/2.0 + i*pixelSizeX - originX;
    pixelY = arma::zeros<arma::vec>(nY);
    for (int i=0; i<nY; i++)
        pixelY[i] = -(physicalSizeY - pixelSizeY)/2.0 + i*pixelSizeY - originY;
    
    if (display) {
        print();
        draw();
    }
}

void arrayGeom::set_xy_offset_m1(int nX, int nY, double pixelScale, bool display) {
    
    physicalSizeX = pixelScale*nX;
    physicalSizeY = pixelScale*nY;
    pixelSizeX = pixelScale;
    pixelSizeY = pixelScale;

    physicalSize = physicalSizeX;

    // matlab: plane.x = ((1:plane.N) - 1 - plane.N/2)/(plane.N)*plane.D;
    // this is not symmetric:
    // i=0 => pixelX = - (nCols/2.)*pixelSizeX
    // i=nCols-1 => pixelX = (nCols - 1 - nCols/2.)*pixelSizeX
    //                  = (nCols/2. - 1)*pixelSizeX
    pixelX = arma::zeros<arma::vec>(nX);
    for (int i=0; i<nX; i++)
        pixelX[i] = ((double) i - ((double)nX)/2.)*pixelSizeX;
    pixelY = arma::zeros<arma::vec>(nY);
    for (int i=0; i<nY; i++)
        pixelY[i] = ((double) i - ((double)nY)/2.)*pixelSizeY;
    
    if (display) {
        print();
        draw();
    }
}

/*
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
*/
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
    printf("...\n");
    for (int i=pixelX.n_elem - 5; i<pixelX.n_elem; i++)
        printf("pixel %d: X = %0.23f, Y = %0.23f\n", i, pixelX[i], pixelY[i]);
}

void arrayGeom::draw(void) {
    draw_mat(pixelXX, pixelX(0), pixelX(pixelX.n_elem-1), pixelY(0), pixelY(pixelY.n_elem-1), "arrayGeom pixelXX", "matlab");
    draw_mat(pixelYY, pixelX(0), pixelX(pixelX.n_elem-1), pixelY(0), pixelY(pixelY.n_elem-1), "arrayGeom pixelYY", "matlab");
    draw_mat(pixelRR, pixelX(0), pixelX(pixelX.n_elem-1), pixelY(0), pixelY(pixelY.n_elem-1), "arrayGeom pixelRR", "matlab");
    draw_mat(pixelTT, pixelX(0), pixelX(pixelX.n_elem-1), pixelY(0), pixelY(pixelY.n_elem-1), "arrayGeom pixelTT", "matlab");
}




