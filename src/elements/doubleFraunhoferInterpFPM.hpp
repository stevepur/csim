//
//  doubleFraunhoferInterpFPM.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef doubleFraunhoferInterpFPM_hpp
#define doubleFraunhoferInterpFPM_hpp

#include "celem.hpp"
#include "armadillo"
#include "../lib/csim_lib.hpp"
#include "complexHexMaskFPM.hpp"

class doubleFraunhoferInterpFPM : public celem {
    zoomFft propZoomFftIn;
    zoomFft propZoomFftOut;
    double zoomFactor = 1.0;
    double padFactor = 1.0;
    
    double fpmScale = -1;
    
    complexHexMaskFPM *hexFPM;

public:
    doubleFraunhoferInterpFPM();
    doubleFraunhoferInterpFPM(initCommandSet*& cmdBlock);
    
    void init(initCommandSet*& cmdBlock);
    void init_geom(void);
    
    void set(std::string fieldName, const char *arg);
    
    efield* execute(efield* E, celem* prev, celem* next, double time);
    
    void print(const char *hdr = "");
    void draw(const char *title = "");
};

#endif /* doubleFraunhoferInterpFPM_hpp */
