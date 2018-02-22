//
//  contrastCurve.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef contrastCurve_hpp
#define contrastCurve_hpp

#include "../efield/efield.hpp"

class contrastCurve {
    int pixelSampling = 1;
    char *filename = NULL;
    bool draw = false;
    double drawTo = -1;
    double specialNormalization = 0;
    
public:

    contrastCurve();
    contrastCurve(initCommandSet*& cmdBlocks);
    
    void init(initCommandSet*& cmdBlocks);
    void set(std::string fieldName, const char *arg);
    
    void execute(efield *inE, double time);
    void make_contrast_curve(void);
};


#endif /* contrastCurve_hpp */
