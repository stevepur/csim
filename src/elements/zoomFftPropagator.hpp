//
//  zoomFftPropagator.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef zoomFftPropagator_hpp
#define zoomFftPropagator_hpp

#include "celem.hpp"
#include "armadillo"
#include "../lib/csim_lib.hpp"

class zoomFftPropagator : public celem {
    zoomFft propZoomFft;
    int dir = 1;
    double zoomFactor = 1.0;
    double padFactor = 1.0;
    
public:
    zoomFftPropagator();
    zoomFftPropagator(initCommandSet*& cmdBlock);
    
    void init(initCommandSet*& cmdBlock);
    
    void set(std::string fieldName, const char *arg);
    
    efield* execute(efield* E, celem* prev, celem* next, double time);
    
    void print(const char *hdr = "");

};

#endif /* zoomFftPropagator_hpp */
