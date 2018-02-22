//
//  doubleFraunhoferFPM.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef doubleFraunhoferFPM_hpp
#define doubleFraunhoferFPM_hpp

#include "celem.hpp"
#include "armadillo"
#include "../lib/csim_lib.hpp"

class doubleFraunhoferFPM : public celem {
    zoomFft propZoomFftIn;
    zoomFft propZoomFftOut;
    double zoomFactor = 1.0;
    double padFactor = 1.0;

    arma::cube complexMaskMatAmp;
    arma::cube complexMaskMatPh;
    arma::cx_cube complexMaskCube;
    int maskIndex = -1;

public:
    doubleFraunhoferFPM();
    doubleFraunhoferFPM(initCommandSet*& cmdBlock);
    
    void init(initCommandSet*& cmdBlock);
    void initMask(const char *filenameRe, const char *filenameIm);
    
    void set(std::string fieldName, const char *arg);
    
    efield* execute(efield* E, celem* prev, celem* next, double time);
    
    void print(const char *hdr = "");
    void draw(const char *title = "");
};

#endif /* doubleFraunhoferFPM_hpp */
