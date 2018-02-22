//
//  celem.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef complexMask_hpp
#define complexMask_hpp

#include "celem.hpp"
#include "armadillo"

class complexMask : public celem {
    arma::cube complexMaskMatAmp;
    arma::cube complexMaskMatPh;
    arma::cx_cube complexMaskCube;
    int maskIndex = -1;
    
public:
    complexMask();
    complexMask(const char *inName);
    complexMask(initCommandSet*& cmdBlock);
    
    void init(initCommandSet*& cmdBlock);
    void init(const char *filenameRe, const char *filenameIm);
    
    void set(std::string fieldName, const char *arg);
    
    efield* execute(efield* E, celem* prev, celem* next, double time);
    
    void draw(const char *title = "");
};

#endif /* binary_complexMask_hpp */
