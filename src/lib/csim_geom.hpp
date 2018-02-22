//
//  csim_geom.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef csim_geom_hpp
#define csim_geom_hpp

#include "armadillo"

class arrayGeom {
public:
    // Physical size of the array
    double physicalSize = 0;  //  (.D in matlab code) only makes sense for square array
    double physicalSizeX = 0;  // (.Dx in matlab code)
    double physicalSizeY = 0;  // (.Dy in matlab code)
    // physical pixel size
    double pixelSizeX = 0;  // (.dx in matlab code)
    double pixelSizeY = 0;  // (.dy in matlab code)
    // physical pixel locations
    arma::vec pixelX;  // (.x in matlab code)
    arma::vec pixelY;  // (.y in matlab code)
    // physical pixel mesh locations
    arma::mat pixelXX;  // (.xx in matlab code)
    arma::mat pixelYY;  // (.yy in matlab code)
    arma::mat pixelRR;  // (.rr in matlab code) distance of pixel from array center
    arma::mat pixelTT;  // (.ttheta in matlab code) angle of pixel from +x axis

    
    arrayGeom() {}
    void set_geometry(int nRows, int nCols, double physicalRadius, bool display = false);
    void set_xy(int nRows, int nCols, double physicalRadius, bool display = false);
    void set_xy_m1(int nRows, int nCols, double physicalRadius, bool display = false);
    void set_geometry(double referenceLambda, double samplesPerFld, double fovInFld, bool display = false);
    void set_mesh(void);
    void print(const char *hdr = "");
    void draw(void);
};


#endif /* csim_geom_hpp */
