//
//  csim.cpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//


#include <iostream>
#include <string>
#include <cstring>
#include <unistd.h>
#include "csim.hpp"
#include "efield.hpp"
#include "coronagraph.hpp"
#include "telescope.hpp"
#include "csim_parser.hpp"
#include "csim_fits.hpp"
#include "csim_fft.hpp"
#include "contrastCurve.hpp"
#include "regionContrast.hpp"
#include "makeHexCFpmResponse.hpp"
#include "differentialEvolutionOptimizer.hpp"

int main(int argn, char **argv) {
    arma::wall_clock timer;
    bool computeFftWisdon = false;
    bool showTimes = false;
    int c;
    
    while ((c = getopt(argn, argv, "wt")) != -1) {
        switch (c) {
            case 'w':
                computeFftWisdon = true;
                break;
                
            case 't':
                showTimes = true;
                break;
                
            default:
                ;
        }
    }
    
    init_fft_lib(computeFftWisdon); // initialize the FFTW library
    initCommandSet initCommands(argv[optind]); // parse the input script into command block sets

    std::vector<initCommandSet*> cmdBlocks = initCommands.find_command_blocks(); // extract the top level of command blocks
//    std::cout << "found " << cmdBlocks.size() << " command blocks" << std::endl;
    for (int i=0; i<cmdBlocks.size(); i++) { // for each command block, respond to the first command
//        cmdBlocks[i]->print();
        if (!strcmp(cmdBlocks[i]->commandList[0]->getCmdStr(), "telescope")) {
            // create the global telescope
            globalTelescope = new telescope(cmdBlocks[i]);
            globalTelescope -> print();
        }
        if (!strcmp(cmdBlocks[i]->commandList[0]->getCmdStr(), "coronagraph")) {
            // create the global coronagraph.  This initializes the entire coronagraph
            globalCoronagraph = new coronagraph(cmdBlocks[i]);
        }
        if (!strcmp(cmdBlocks[i]->commandList[0]->getCmdStr(), "efield")) {
            // create and initialize the initial e-field
            initialEfield = new efield(cmdBlocks[i]);
        }
        if (!strcmp(cmdBlocks[i]->commandList[0]->getCmdStr(), "makeContrastCurve")) {
            // create and execute the contrast curve tool
            contrastCurve *myContrastCurve = new contrastCurve(cmdBlocks[i]);
            myContrastCurve->make_contrast_curve();
        }
        if (!strcmp(cmdBlocks[i]->commandList[0]->getCmdStr(), "regionContrast")) {
            // create and execute the contrast curve tool
            regionContrast *myRegionContrast = new regionContrast(cmdBlocks[i]);
            myRegionContrast->compute_contrast();
        }
        if (!strcmp(cmdBlocks[i]->commandList[0]->getCmdStr(), "makeHexCFpmResponse")) {
            // create and execute the contrast curve tool
            makeHexCFpmResponse *myMakeHexCFpmResponse = new makeHexCFpmResponse(cmdBlocks[i]);
            myMakeHexCFpmResponse->compute_response();
        }
        if (!strcmp(cmdBlocks[i]->commandList[0]->getCmdStr(), "differentialEvolutionOptimizer")) {
            // create and execute the contrast curve tool
            differentialEvolutionOptimizer *deOptimizer = new differentialEvolutionOptimizer(cmdBlocks[i]);
            deOptimizer->optimize();
        }
        if (!strcmp(cmdBlocks[i]->commandList[0]->getCmdStr(), "execute")) {
            // do a single execution of the global coronagraph
            timer.tic();
            globalCoronagraph->execute(initialEfield, 0, showTimes);
            std::cout << "csim execution time: " << timer.toc() << " seconds" << std::endl;
        }
    }
    if (computeFftWisdon)
        save_fft_wisdom();
}
