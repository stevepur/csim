//
//  celem.cpp
//  csim
//
//  Created by steve on 5/24/17.
//  Copyright © 2017 NASA. All rights reserved.
//

#include "celem.hpp"
#include <assert.h>

celem::celem() {
}

efield* celem::execute(efield* E, celem* prev, celem* next, double time) {
    std::cout << "executed a celem" << std::endl;
    return E;
}

void celem::pre_execute(efield* E, celem* prev, celem* next, double time) {
    if (preExecuteDraw) {
        char str[200];
        sprintf(str, "%s input efield", name);
        E->draw(str);
    }
    if (preExecuteSave) {
        char str[200];
        sprintf(str, "%s_in", name);
        E->save(str);
    }
    if (preExecuteBreak) {
        assert(NULL);
    }
}

void celem::post_execute(efield* E, celem* prev, celem* next, double time) {
    if (postExecuteDraw) {
        char str[200];
        sprintf(str, "%s output efield", name);
        E->draw(str);
    }
    if (postExecuteSave) {
        char str[200];
        sprintf(str, "%s_out", name);
        E->save(str);
    }
    if (postExecuteBreak) {
        assert(NULL);
    }
}

void celem::post_init(void) {
    if (onInitDraw) {
        char str[200];
        sprintf(str, "%s initial", name);
        draw(str);
    }
}


bool celem::set(std::string fieldName, const char *arg) {
    
    if (fieldName == "position") {
        // arg is one double value
        position = atof(arg);
    }
    else if (fieldName == "name") {
        // arg is one double value
        name = new char[strlen(arg)+1];
        strcpy(name, arg);
    }
    else if (fieldName == "disableForCalibration") {
        // arg is a string
        if (!strcmp(arg, "on"))
            disableForCalibration = true;
        else if (!strcmp(arg, "off")) {
            disableForCalibration = false;
        } else
            std::cout << "!!! bad disableForCalibration arg" << std::endl;
    }
    else if (fieldName == "draw") {
        // arg is a string
        if (!strcmp(arg, "preExecute"))
            preExecuteDraw = true;
        else if (!strcmp(arg, "postExecute"))
            postExecuteDraw = true;
        else if (!strcmp(arg, "onInit"))
            onInitDraw = true;
        else if (!strcmp(arg, "off")) {
            onInitDraw = false;
            preExecuteDraw = false;
            postExecuteDraw = false;
        } else
            std::cout << "!!! bad draw arg" << std::endl;
    }
    else if (fieldName == "save") {
        // arg is a string
        if (!strcmp(arg, "preExecute"))
        preExecuteSave = true;
        else if (!strcmp(arg, "postExecute"))
        postExecuteSave = true;
        else if (!strcmp(arg, "off")) {
            preExecuteSave = false;
            postExecuteSave = false;
        } else
        std::cout << "!!! bad save arg" << std::endl;
    }
    else if (fieldName == "break") {
        // arg is a string
        if (!strcmp(arg, "preExecute"))
        preExecuteBreak = true;
        else if (!strcmp(arg, "postExecute"))
        postExecuteBreak = true;
        else if (!strcmp(arg, "off")) {
            postExecuteBreak = false;
        } else
        std::cout << "!!! bad break arg" << std::endl;
    }
    else
        return false;
    
    return true;
}

void celem::print(const char *hdr) {
    std::cout << "celem: " << hdr << std::endl;
    std::cout << "name: " << name << std::endl;
    std::cout << "position: " << position << std::endl;
}

void celem::draw(const char *title) {
    std::cout << "can't draw a virtual celem" << std::endl;
}

