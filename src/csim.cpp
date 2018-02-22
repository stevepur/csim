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
#include "csim.hpp"
#include "efield.hpp"
#include "coronagraph.hpp"
#include "telescope.hpp"
#include "csim_parser.hpp"
#include "csim_fits.hpp"
#include "csim_fft.hpp"
#include "contrastCurve.hpp"

int main(int argn, char **argv) {
    arma::wall_clock timer;
    
    init_fft_lib(); // initialize the FFTW library
    initCommandSet initCommands(argv[1]); // parse the input script into command block sets
    
    std::vector<initCommandSet*> cmdBlocks = initCommands.find_command_blocks(); // extract the top level of command blocks
    for (int i=0; i<cmdBlocks.size(); i++) { // for each command block, respond to the first command
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
        if (!strcmp(cmdBlocks[i]->commandList[0]->getCmdStr(), "execute")) {
            // do a single execution of the global coronagraph
            timer.tic();
            globalCoronagraph->execute(initialEfield, 0);
            std::cout << "csim execution time: " << timer.toc() << " seconds" << std::endl;
        }
    }
//    save_fft_wisdom();
}
