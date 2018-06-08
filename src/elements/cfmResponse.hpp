//
//  celem.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef cfmResponse_hpp
#define cfmResponse_hpp

#include "celem.hpp"
#include "armadillo"
#include "responseData.hpp"

class cfmResponse : public celem {
    responseData cfmResponseData;
    arma::mat sagMat;
    arma::cx_mat phaseMat;
    arma::cx_mat respMat;
    
public:
    cfmResponse();
    cfmResponse(const char *inName);
    cfmResponse(initCommandSet*& cmdBlock);
    
    void init(initCommandSet*& cmdBlock);
    void init(const char *filename);
    
    void set(std::string fieldName, const char *arg);
    
    efield* execute(efield* E, celem* prev, celem* next, double time);
    
    void draw(const char *title = "");
};

#endif /* cfmResponse_hpp */
