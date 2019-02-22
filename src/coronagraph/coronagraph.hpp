//
//  coronagraph.hpp
//  csim
//
// The coronagraph class is the manager and container for the coronagraph's
// components. This class is responsible for
//  - parsing input text to create celem coronagraph component objects (init method)
//  - holding the list of celem components as a doubly linked list (elemList)
//  - propagating an efield object through the coronagraph (execute mehtod)
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef coronagraph_hpp
#define coronagraph_hpp

#include <list>
#include "../elements/celem.hpp"
#include "../data/efield.hpp"
#include "../lib/csim_parser.hpp"

// optimization data class
class optData {
public:
    char *componentName = NULL;
    char *dataName = NULL;
    int dataLength = 0;
    double lowerBound = 0.0;
    double upperBound = 0.0;
    
    optData(const char *componentName, const char *dataName, int dataLength, double lb, double ub);
    void print(const char *hdr = "");
};

// The coronagraph class is a sub-class of celem, so coronagraphs can be nested.
class coronagraph : public celem {
    // Holder for the propagated efield, set by execute method.
    efield *E = NULL;
    // Double-linked list of celem objects, defining the coronagraph.
    // This list is created by the init method based on the input text file.
    // This list is traversed in order by the execute method to execute
    // the coronagraph.
    std::list<celem*> elemList;
    // State flag indicating whether or not the coronagraph is being executed
    // in calibration mode.  Different celem objects respond to this flag
    // differently.
    bool calibrating = false;
    
    // list of optData objects, allowing multiple data to be subject to joint optimization
    std::list<optData*> optDataList;
    int optDataListDataLength = 0;

public:
    // Constructor that accepts an array of cmdBlocks objects, each of which
    // contains a keyword string and value string.
    // Passes the cmdBlocks argument to the init method.
    coronagraph(initCommandSet*& cmdBlocks);
    // Create the coronagraph by filling in elemList based on the inputs in
    // cmdBlocks.
    void init(initCommandSet*& cmdBlocks);

    // setters/getters for the calibrationState flag
    void set_calibration_state(bool calibrationState) {calibrating = calibrationState;}
    bool get_calibration_state(void) {return calibrating;}
    
    // setters/getters for the optimization data
    void add_optimization_data(const char *componentName, const char *dataName, double lb, double ub);
    void get_optimization_bounds(arma::vec& lowerBoundData, arma::vec& upperBoundData);
    void get_optimization_data(arma::vec& allOptData);
    void set_optimization_data(arma::vec& allOptData);
    void get_optimization_data(const char *componentName, const char *dataName, arma::vec& data);
    void set_optimization_data(const char *componentName, const char *dataName, arma::vec& data);
    void save_optimization_data(const char *componentName, const char *dataName, char *outputDirectory = NULL);

    // Execute the coronagraph by calling the execute method
    // of each celem object in elemList, passing the efield pointer E.
    // Modifies the efield pointed to by E.
    void execute(efield *inE, double time, bool showTimes = false);
    
    void print(void);
};

// global pointer to the top-level coronagraph
extern coronagraph *globalCoronagraph;

#endif /* coronagraph_hpp */
