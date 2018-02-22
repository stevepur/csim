//
//  coronagraph.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef telescope_hpp
#define telescope_hpp

// #include celem.hpp
#include <iostream>
#include "../lib/csim_parser.hpp"



class telescope {
    double primaryDiameter = 0;
    double primaryfRatio = 0;
    double primaryfLength = 0;
    
public:
    telescope(initCommandSet*& cmdBlocks);
    
    void init(initCommandSet*& cmdBlocks);
    
    void set(std::string fieldName, const char *value);
    double get(std::string fieldName);
    
    void print(const char *header = "");
};

extern telescope *globalTelescope;

#endif /* telescope_hpp */
