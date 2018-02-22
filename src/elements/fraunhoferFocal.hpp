//
//  fraunhoferFocal.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef fraunhoferFocal_hpp
#define fraunhoferFocal_hpp

#include "celem.hpp"
#include "armadillo"
#include "../lib/csim_lib.hpp"

class fraunhoferFocal : public celem {
    zoomFft propZoomFft;
    arrayGeom arrayGeometry;
    double referenceLambda = 0;
    double samplesPerflD = 0;
    double FOVflD = 0;
    double scaleFlD = 0;
    double focalLength = 0;
    
public:
    fraunhoferFocal();
    fraunhoferFocal(initCommandSet*& cmdBlock);
    
    void init(initCommandSet*& cmdBlock);
    
    void set(std::string fieldName, const char *arg);
    void set_geometry(void);
    
    efield* execute(efield* E, celem* prev, celem* next, double time);
    
    void print(const char *hdr = "");

};

#endif /* fraunhoferFocal_hpp */
