//
//  coronagraph.cpp
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

#include <iostream>
#include <string>
#include <cstring>

#include "coronagraph.hpp"
#include "../elements/breakExecution.hpp"
#include "../elements/mask.hpp"
#include "../elements/mirror.hpp"
#include "../elements/deformableMirror.hpp"
#include "../elements/complexMask.hpp"
#include "../elements/fresnelPropagator.hpp"
#include "../elements/fraunhoferFocal.hpp"
#include "../elements/doubleFraunhoferFPM.hpp"
#include "../elements/doubleFraunhoferInterpFPM.hpp"
#include "../elements/fpmPupToLyot.hpp"
#include "../elements/zoomDftPropagator.hpp"
#include "../elements/zoomFftPropagator.hpp"
#include "../elements/downsample.hpp"
#include "../elements/upsample.hpp"

// global pointer to the top-level coronagraph
coronagraph *globalCoronagraph = NULL;

// Constructor that accepts an array of cmdBlocks objects, each of which
// contains a keyword string and value string.
// Passes the cmdBlocks argument to the init method.
coronagraph::coronagraph(initCommandSet*& cmdBlocks) {
    init(cmdBlocks);
}

// Create the coronagraph by filling in elemList based on the inputs in
// cmdBlocks.  Each command block contains an array of commandList objects,
// each of which has getCmdStr method, which returns
// that command's keyword. At this level the only keywords seen should be
// names of celem classes, which will be the first command in the cmdBlock list.
void coronagraph::init(initCommandSet*& cmdBlocks) {
//    std::cout << "initing a coronagraph" << std::endl;

    // Find the command blocks in the input command blocks.
    // At this point there will be several command blocks, each of which
    // contains a list of commands defining a celem: the first command
    // specifies the celem subclass, in which case this command block is
    // passed to the constructor for that celem.  The celem then parses the
    // commands in teh command block which define parameters of that celem.

    // an initCommandSet contains an array commandList of commands.
    std::cout << "coronagraph::init" << std::endl;
    std::vector<initCommandSet*> subBlocks = cmdBlocks->find_command_blocks();
//    std::cout << "found " << subBlocks.size() << " subblocks" << std::endl;

    // walk through all the command bocks, creating celem classes as specified
    // by the first command in each block
    for (int j=0; j<subBlocks.size(); j++) {
//        char str[200];
//        sprintf(str, "-------------------- block %d", j);
//        subBlocks[j]->print(str);
        std::cout << "creating a " << subBlocks[j]->commandList[0]->getCmdStr() << std::endl;
        // empty celem
        if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "coronagraph")) {
            ;
        }
        // break execution celem (stops coronagraph execution)
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "optElem")) {
            celem *newE = new celem();
            elemList.push_back(newE);
        }
        // break execution celem (stops coronagraph execution)
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "breakExecution")) {
            celem *newE = new breakExecution(subBlocks[j]);
            elemList.push_back(newE);
        }
        // fresnel propagator celem
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "fresnelPropagator")) {
            celem *newE = new fresnelPropagator(subBlocks[j]);
            elemList.push_back(newE);
        }
        // grey-scale (multiplicative) mask celem
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "mask")) {
            celem *newE = new mask(subBlocks[j]);
            elemList.push_back(newE);
        }
        // grey-scale (multiplicative) mask celem
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "complexMask")) {
            celem *newE = new complexMask(subBlocks[j]);
            elemList.push_back(newE);
        }
        // prefect reflector mirror celem
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "mirror")) {
            celem *newE = new mirror(subBlocks[j]);
            elemList.push_back(newE);
        }
        // deformable prefect reflector mirror celem
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "deformableMirror")) {
            celem *newE = new deformableMirror(subBlocks[j]);
            elemList.push_back(newE);
        }
        // fraunhoferFocal, which propagates to a focal plane via fraunhofer
        // propagation
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "fraunhoferFocal")) {
            celem *newE = new fraunhoferFocal(subBlocks[j]);
            elemList.push_back(newE);
        }
        // zoomDftPropagator, which propagates via zoomDft
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "zoomDftPropagator")) {
            celem *newE = new zoomDftPropagator(subBlocks[j]);
            elemList.push_back(newE);
        }
        // zoomFftPropagator, which propagates via zoomFft
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "zoomFftPropagator")) {
            celem *newE = new zoomFftPropagator(subBlocks[j]);
            elemList.push_back(newE);
        }
        // doubleFraunhoferFPM, which does zoomFft, complex mask, zoomFft
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "doubleFraunhoferFPM")) {
            celem *newE = new doubleFraunhoferFPM(subBlocks[j]);
            elemList.push_back(newE);
        }
        // doubleFraunhoferFPM, which does zoomFft, on-the-fly-intperolated complex mask, zoomFft
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "doubleFraunhoferInterpFPM")) {
            celem *newE = new doubleFraunhoferInterpFPM(subBlocks[j]);
            elemList.push_back(newE);
        }
        // fpmPupToLyot, which propagates from a pupil through a focal plane
        // mask to a pupil (typically a Lyot stop).  Uses Fraunhofer propagation.
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "fpmPupToLyot")) {
            celem *newE = new fpmPupToLyot(subBlocks[j]);
            elemList.push_back(newE);
        }
        // downsample, which downsamples the E-field
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "downsample")) {
            celem *newE = new downsample(subBlocks[j]);
            elemList.push_back(newE);
        }
        // downsample, which downsamples the E-field
        else if (!strcmp(subBlocks[j]->commandList[0]->getCmdStr(), "upsample")) {
            celem *newE = new upsample(subBlocks[j]);
            elemList.push_back(newE);
        }
        else {
            std::cout << "!!! coronagraph init unknown object name: " << subBlocks[j]->commandList[0]->getCmdStr() << std::endl;
        }
#if 0
        // dump the remaining commands for debugging
        else {
            for (int sc=0; sc<subBlocks[j]->commandList.size(); sc++) {
                std::cout << "block " << j << " command " << sc << " "
                    << subBlocks[j]->commandList[sc]->getCmdStr() << ":"
                    << subBlocks[j]->commandList[sc]->getArgStr() << std::endl;
            }
        }
#endif
    }
//    std::cout << "elemList contains " << elemList.size() << " entries" << std::endl;
}

// Execute the coronagraph by calling the execute function of each celem in elemList.
// The celem method accepts (efield pointer, previous celem, next celem, time parameter).
// The previous and next celems are supplied for celems like propagatorss that need to know
// the physical distance between the previous and next celems.
// MODIFIES THE CONTENTS OF THE EFIELD POINTED TO BY inE
void coronagraph::execute(efield *inE, double time, bool showTimes) {
    arma::wall_clock timer;

    E = inE;

    // traverse elemList
    for (std::list<celem*>::iterator it = elemList.begin(); it != elemList.end(); ++it) {
//        std::cout << "======================= executing " << (*it)->name <<  std::endl;
        timer.tic();
        // if we're the first celem in elemList...
        if (it == elemList.begin()) {
//            std::cout << "first element" << std::endl;
            E = (*it)->execute(E, NULL, *std::next(it), time);
        }
        // if we're the last celem in elemList...
        else if (it == std::prev(elemList.end())){
//            std::cout << "last element" << std::endl;
            E = (*it)->execute(E, *std::prev(it), NULL, time);
        }
        else // we're inside the list
            E = (*it)->execute(E, *std::prev(it), *std::next(it), time);
        
        if (showTimes)
            std::cout << "======================= executing " << (*it)->name << " time: " << timer.toc() << " seconds" << std::endl;
    }
}

void coronagraph::add_optimization_data(const char *componentName, const char *dataName, double lb, double ub) {
    arma::vec tempVec;
    
    get_optimization_data(componentName, dataName, tempVec);
    optData *newOptData = new optData(componentName, dataName, tempVec.n_elem, lb, ub);
    optDataListDataLength += newOptData->dataLength;
    std::cout << "optDataListDataLength = " << optDataListDataLength << std::endl;
    
    optDataList.push_back(newOptData);
    
    for (std::list<optData*>::iterator it = optDataList.begin(); it != optDataList.end(); ++it) {
        (*it)->print();
    }
}

void coronagraph::get_optimization_bounds(arma::vec& lowerBoundData, arma::vec& upperBoundData) {
    lowerBoundData.set_size(optDataListDataLength);
    upperBoundData.set_size(optDataListDataLength);
    int dataCount = 0;
    arma::vec tempVec;
    for (std::list<optData*>::iterator it = optDataList.begin(); it != optDataList.end(); ++it) {
        get_optimization_data((*it)->componentName, (*it)->dataName, tempVec);
        lowerBoundData.subvec(dataCount, dataCount + (*it)->dataLength - 1) = (*it)->lowerBound*arma::ones<arma::vec>((*it)->dataLength);
        upperBoundData.subvec(dataCount, dataCount + (*it)->dataLength - 1) = (*it)->upperBound*arma::ones<arma::vec>((*it)->dataLength);
        dataCount += (*it)->dataLength;
    }
}

void coronagraph::get_optimization_data(arma::vec& allOptData) {
    allOptData.set_size(optDataListDataLength);
    int dataCount = 0;
    arma::vec tempVec;
    for (std::list<optData*>::iterator it = optDataList.begin(); it != optDataList.end(); ++it) {
        get_optimization_data((*it)->componentName, (*it)->dataName, tempVec);
        allOptData.subvec(dataCount, dataCount + (*it)->dataLength - 1) = tempVec;
        dataCount += (*it)->dataLength;
    }
}

void coronagraph::set_optimization_data(arma::vec& allOptData) {
    int dataCount = 0;
    arma::vec tempVec;
    for (std::list<optData*>::iterator it = optDataList.begin(); it != optDataList.end(); ++it) {
        tempVec = allOptData.subvec(dataCount, dataCount + (*it)->dataLength - 1);
        set_optimization_data((*it)->componentName, (*it)->dataName, tempVec);
        dataCount += (*it)->dataLength;
    }
}


// get optimization data dataName from the component componentName
void coronagraph::get_optimization_data(const char *componentName, const char *dataName, arma::vec& data) {
    for (std::list<celem*>::iterator it = elemList.begin(); it != elemList.end(); ++it) {
        if (!strcmp((*it)->name, componentName))
            (*it)->get_optimization_data(dataName, data);
    }
}

// set optimization data dataName from the component componentName
void coronagraph::set_optimization_data(const char *componentName, const char *dataName, arma::vec& data) {
    for (std::list<celem*>::iterator it = elemList.begin(); it != elemList.end(); ++it) {
        if (!strcmp((*it)->name, componentName)) {
            (*it)->set_optimization_data(dataName, data);
        }
    }
}

// set optimization data dataName from the component componentName
void coronagraph::save_optimization_data(const char *componentName, const char *dataName, char *outputDirectory) {
    for (std::list<celem*>::iterator it = elemList.begin(); it != elemList.end(); ++it) {
        if (!strcmp((*it)->name, componentName)) {
            (*it)->save_optimization_data(dataName, outputDirectory);
        }
    }
}

// Print the coronagraph
void coronagraph::print(void) {
    int elementCount = 0;
    
    for (std::list<celem*>::iterator it = elemList.begin(); it != elemList.end(); ++it) {
        std::cout << "=============== coronagraph element " << elementCount << std::endl;
        (*it)->print();
        elementCount++;
    }
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

optData::optData(const char *cName, const char *dName, int nElem, double lb, double ub) {
    if (componentName != NULL) {
        delete[] componentName;
        componentName = NULL;
    }
    componentName = new char[strlen(cName)+1];
    strcpy(componentName, cName);

    if (dataName != NULL) {
        delete[] dataName;
        dataName = NULL;
    }
    dataName = new char[strlen(dName)+1];
    strcpy(dataName, dName);
    
    dataLength = nElem;
    lowerBound = lb;
    upperBound = ub;
}

void optData::print(const char *hdr) {
    std::cout << "optData " << hdr << ":" << std::endl;
    std::cout << "componentName: " << componentName << std::endl;
    std::cout << "dataName: " << dataName << std::endl;
    std::cout << "dataLength: " << dataLength << std::endl;
    std::cout << "lowerBound: " << lowerBound << std::endl;
    std::cout << "upperbound: " << upperBound << std::endl;
}
