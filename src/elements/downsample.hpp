//
//  celem.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef downsample_hpp
#define downsample_hpp

#include "celem.hpp"
#include "armadillo"

class downsample : public celem {
    int sampleFactor = 1;
   
public:
    downsample();
    downsample(initCommandSet*& cmdBlock);
    
    void init(initCommandSet*& cmdBlock);
    
    void set(std::string fieldName, const char *arg);
    
    efield* execute(efield* E, celem* prev, celem* next, double time);
};

#endif /* downsample_hpp */
