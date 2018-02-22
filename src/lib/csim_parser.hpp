//
//  csim_parser.hpp
//  csim
//
//  Created by steve on 2/20/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#ifndef csim_parser_hpp
#define csim_parser_hpp

#define nCommandsDefault 20

#include <vector>

class initCommand {
    char *cmdStr;
    char *argStr;
    
public:
    initCommand();
    initCommand(const char *cStr, const char *aStr);
    initCommand(const initCommand& c);
    ~initCommand();
    void setCommand(const char *cStr, const char *aStr);
    char *getCmdStr() { return cmdStr; }
    char *getArgStr() { return argStr; }
    void setCommand(const initCommand& c);
    void print(const char *header = "");
};

class initCommandSet {
public:
    std::vector<initCommand*> commandList;
    
    initCommandSet();
    initCommandSet(char *filename);
    ~initCommandSet();
    
    void add_command(std::string cmdArgStr);
    void add_command(initCommand *cmd);
    void add_command(const char *cStr, const char *aStr);
    std::vector<initCommandSet*> find_command_blocks(void);
    void print(const char *header = "");
};

#endif /* csim_parser_hpp */
