/*
 * cshell.h
 *
 *  Created on: May 2016
 *      Author: JML
 */

#pragma once

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <list>
#include <sstream>
#include <time.h>
using namespace std;

#include "io/CD3d_io.h"
#include "util/CD3d_util.h"
#include "CD3d_stat.h"
#include "util/model/objReader.h"
#include "util/gettime.h"

#ifndef _CORISELIB_
#include "dude3d_draw.h"
#endif

///////////////////////////////////////////////////////////////////////////////
// Prototypes
inline bool parseCmd(int argc,char ** argv);
inline void printUsage(char * arg0);
inline void createModels(const list<string>& filenames, ML& ml);
inline void model_COM_R(ML& models);
///////////////////////////////////////////////////////////////////////////////
