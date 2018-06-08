//
//  makeHexCFpmResponse.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef makeHexCFpmResponse_hpp
#define makeHexCFpmResponse_hpp

#include "../data/efield.hpp"
#include "regionContrast.hpp"

class makeHexCFpmResponse {
    double referenceLambda = -1.0;
    char *componentName = NULL;
    char *dataName = NULL;
    char *filename = NULL;
    bool draw = false;
    regionContrast *region = NULL;
    char *outputDirectory = NULL;
    
public:

    makeHexCFpmResponse();
    makeHexCFpmResponse(initCommandSet*& cmdBlocks);
    
    void init(initCommandSet*& cmdBlocks);
    void init_block(initCommandSet*& cmdBlock);
    void set(std::string fieldName, const char *arg);
    
    void execute(efield *inE, double time);
    void compute_response(void);
};


#endif /* makeHexCFpmResponse_hpp */
