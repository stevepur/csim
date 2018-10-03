//
//  celem.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef binary_mask_hpp
#define binary_mask_hpp

#include "celem.hpp"
#include "armadillo"

class mask : public celem {
    arma::mat maskMat;
    arma::cx_mat cxMaskMat;
    
public:
    mask();
    mask(const char *inName);
    mask(initCommandSet*& cmdBlock);
    
    void init(initCommandSet*& cmdBlock);
    void init(const char *filename);
    
    void set(std::string fieldName, const char *arg);
    
    efield* execute(efield* E, celem* prev, celem* next, double time);
    
    void draw(const char *title = "");
};

#endif /* binary_mask_hpp */
