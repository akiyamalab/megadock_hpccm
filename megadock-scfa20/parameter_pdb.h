/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : ParameterPDB
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef ParameterPDB_h
#define ParameterPDB_h 1

#include "parameter.h"

using namespace std;

class ParameterPDB : public Parameter
{
private:
    friend class          ControlPDB;
    friend class          DockingPDB;
    friend class          FFTProcessPDB;
    string            _RecPDB_file;
    string            _LigPDB_file;
    int               _IO_flag[3];

protected:
    virtual void          pdb_step();

public:
    ParameterPDB(Parallel *pparallel) : Parameter(pparallel) {
#ifdef DEBUG
        cout << "Constructing ParameterPDB.\n";
#endif
    }
    virtual           ~ParameterPDB() {
#ifdef DEBUG
        cout << "Destructing ParameterPDB.\n";
#endif
    }
    virtual void          process_args(int argc,char *argv[]);
};

#endif
