//
//  celem.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef fresnelPropagator_hpp
#define fresnelPropagator_hpp

#include "celem.hpp"
#include "armadillo"
#include "../lib/csim_lib.hpp"

class fresnelPropagator : public celem {
    fresnelPropagateAS propagator;
    double propagationSign = 1;
    
public:
    fresnelPropagator();
    fresnelPropagator(initCommandSet*& cmdBlock);
    
    void init(initCommandSet*& cmdBlock);
    
    void set(std::string fieldName, const char *arg);
    
    efield* execute(efield* E, celem* prev, celem* next, double time);
    
};

#endif /* fresnelPropagator_hpp */
