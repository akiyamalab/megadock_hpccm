/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : ControlPDB
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef ControlPDB_h
#define ControlPDB_h 1

#include "control.h"
#include "cpu_time.h"
#include "parameter_pdb.h"
#include "docking.h"

using namespace std;

class ControlPDB : public Control<ParameterPDB, DockingPDB>
{
private:
    CPUTime   *_cputime;
protected:
    virtual void  autogridr(const int &ngrid,vector<int> &ngrid_table);
    virtual void  autogridl(const int &ngrid,vector<int> &ngrid_table);
    virtual void  checkgridr();
    virtual void  checkgridl();
public:
    ControlPDB(CPUTime *pcputime,Parallel *pparallel)
        : _cputime(pcputime),Control<ParameterPDB, DockingPDB>(pparallel) {
#ifdef DEBUG
        cout << "Constructing ControlPDB.\n";
#endif
    }
    virtual ~ControlPDB() {
#ifdef DEBUG
        cout << "Destructing ControlPDB.\n";
#endif
    }
    virtual void  initialize(int argc,char *argv[]);
    virtual void  execute();
};

#endif
