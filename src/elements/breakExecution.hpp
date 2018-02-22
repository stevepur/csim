//
//  celem.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef binary_breakExecution_hpp
#define binary_breakExecution_hpp

#include "celem.hpp"
#include "armadillo"

class breakExecution : public celem {

public:
    breakExecution();
    breakExecution(initCommandSet*& cmdBlock);
    
    void init(initCommandSet*& cmdBlock);
    
    void set(std::string fieldName, const char *arg);
    
    efield* execute(efield* E, celem* prev, celem* next, double time);
    };

#endif /* binary_breakExecution_hpp */
