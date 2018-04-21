//
//  csim_parser.cpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>

#include "csim_parser.hpp"

initCommand::initCommand() {
    cmdStr = new char[1];
    cmdStr[0] = '\0';
    argStr = new char[1];
    argStr[0] = '\0';
}

initCommand::initCommand(const char *cStr, const char *aStr) {
    setCommand(cStr, aStr);
}

initCommand::initCommand(const initCommand& c) {
    setCommand(c);
}

initCommand::~initCommand() {
    delete []cmdStr;
    delete []argStr;
}

void initCommand::setCommand(const char *cStr, const char *aStr) {
    delete []cmdStr;
    cmdStr = new char[strlen(cStr) + 1];
    std::strcpy(cmdStr, cStr);
    delete []argStr;
    argStr = new char[strlen(aStr) + 1];
    std::strcpy(argStr, aStr);
}

void initCommand::setCommand(const initCommand& c) {
    
    cmdStr = new char[strlen(c.cmdStr) + 1];
    std::strcpy(cmdStr, c.cmdStr);
    argStr = new char[strlen(c.argStr) + 1];
    std::strcpy(argStr, c.argStr);
}

void initCommand::print(const char *header) {
    std::cout << header << cmdStr << " : " << argStr << "\n";
}

///////////////////////////////

initCommandSet::initCommandSet(char *filename) {
    
    std::string line;
    
    std::ifstream cmdFile(filename);
    
    if (cmdFile.is_open()) {
        while (std::getline(cmdFile, line)) {
//            std::cout << ">>>>>>>>> " << line << std::endl;
            // remove leading spaces or tabs from cStr
            int contentPos = line.find_first_not_of(" \t");
            if (contentPos > 0)
                line.erase(0, contentPos);
            // add the command if it is not a comment
            if (line[0] != '#' && line[0] != '\n' && line[0] != '\0') {
                add_command(line);
            }
        }
        cmdFile.close();
    } else {
        std::cout << "unable to open input file\n";
    }
}

initCommandSet::initCommandSet() {
}

initCommandSet::~initCommandSet() {
}


void initCommandSet::add_command(std::string cmdArgStr) {
    
//    std::cout << "adding " << cmdArgStr << std::endl;
    int colonPos = cmdArgStr.find(":");
    
    char cStr[300];
    char aStr[300];
    if (colonPos != -1) {
        int nCopied = cmdArgStr.copy(cStr, colonPos, 0);
        cStr[nCopied] = '\0';
        
        nCopied = cmdArgStr.copy(aStr, cmdArgStr.length()-colonPos, colonPos+1);
        aStr[nCopied] = '\0';
        
        add_command(cStr, aStr);
        
    } else { // assume this is just a command string
        
        
        std::strcpy(cStr, cmdArgStr.c_str());
        aStr[0] = '\0';
        add_command(cStr, aStr);
        
    }
}

void initCommandSet::add_command(initCommand *cmd) {
//    cmd->print();
    add_command(cmd->getCmdStr(), cmd->getArgStr());
    
}

void initCommandSet::add_command(const char *cStr, const char *aStr) {
    initCommand *newC = new initCommand;
    
//    std::cout << "adding " << cStr << " and " << aStr << std::endl;
    newC->setCommand(cStr, aStr);
    commandList.push_back(newC);
}

std::vector<initCommandSet*> initCommandSet::find_command_blocks(void) {
    int startCount = 0;
    std::vector<initCommandSet*> commandBlockList;
    initCommandSet *newBlock = NULL;
    initCommandSet *cmdBlock = NULL; // for commands not inside a block
    
    // Blocks are defined as what's between the start and end keywords.
    // Blocks may contain sub-blocks.
//    std::cout << "subblock search command list size: " << commandList.size() << std::endl;
    for (unsigned int i=0; i<commandList.size(); i++) {
        bool madeBlock = false;
        
        if (!strcmp(commandList[i]->getCmdStr(), "start")) {

            if (!startCount) {
                newBlock = new initCommandSet;
                commandBlockList.push_back(newBlock);
//                std::cout << "======= started a subblock" << std::endl;
            } else
                newBlock->add_command(commandList[i]);
            startCount++;
//            std::cout << "startCount = " << startCount << std::endl;
            madeBlock = true;
        } else if (!strcmp(commandList[i]->getCmdStr(), "end")) {
            startCount--;
//            std::cout << "startCount = " << startCount << std::endl;
            if (startCount)
                newBlock->add_command(commandList[i]);
//            else
//                std::cout << "======= ended a subblock" << std::endl;
        } else if (!startCount & !madeBlock) { // we've found a command not between a start and end
            if (cmdBlock == NULL) // init the cmdBlock if it's not already inited
                cmdBlock = new initCommandSet;
            cmdBlock->add_command(commandList[i]); // add this command to the non-block commands
        } else if (newBlock != NULL) {
            newBlock->add_command(commandList[i]);
            madeBlock = false;
        }
    }
    if (cmdBlock != NULL) { // if there are commands not between start and stop
        commandBlockList.insert(commandBlockList.begin(), cmdBlock);
        std::cout << "inserted a non-delimited command block at the beginning" << std::endl;
    }
    return commandBlockList;
}

void initCommandSet::print(const char *header) {
    std::cout << header << std::endl;
    for (unsigned int i=0; i<commandList.size(); i++)
        commandList[i]->print();
}
