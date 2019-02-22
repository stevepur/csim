//
//  celem.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef celem_hpp
#define celem_hpp

#include "../data/efield.hpp"
#include "../lib/csim_parser.hpp"

class celem {
public:
    char *name = NULL;
    double position = 0;
    bool disableForCalibration = false;
    bool calibrating = false;
    
    bool verbose = false;
    bool onInitDraw = false;
    bool preExecuteDraw = false;
    bool postExecuteDraw = false;
    bool preExecuteSave = false;
    bool postExecuteSave = false;
    bool preExecuteBreak = false;
    bool postExecuteBreak = false;
    
    celem();
    celem(initCommandSet*& cmdBlock) {}
    
    virtual efield* execute(efield* E, celem* prev, celem* next, double time);
    void pre_execute(efield* E, celem* prev, celem* next, double time);
    void post_execute(efield* E, celem* prev, celem* next, double time);
    void post_init(void);
    bool set(std::string fieldName, const char *arg);
    void get_optimization_data(const char *dataName, arma::vec& data);
    void set_optimization_data(const char *dataName, arma::vec& data);
    void save_optimization_data(const char *dataName, char *outputDirectory = NULL);
    void print(const char *hdr = "");
    virtual void draw(const char *title);
};

#endif /* celem_hpp */
