//
//  celem.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef mirror_hpp
#define mirror_hpp

#include "celem.hpp"
#include "armadillo"

class mirror : public celem {
    arma::mat mirrorMat;
    arma::cx_cube mirrorPhase;
    double mirrorSign = 1;
   
public:
    mirror();
    mirror(const char *inName);
    mirror(initCommandSet*& cmdBlock);
    
    void init(initCommandSet*& cmdBlock);
    void init(const char *filename);
    
    void set(std::string fieldName, const char *arg);
    void set_mirrorPhase(void);

    efield* execute(efield* E, celem* prev, celem* next, double time);
    
    void draw(const char *title = "");
};

#endif /* mirror_hpp */
