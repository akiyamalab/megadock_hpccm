/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : ParameterTable
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef ParameterTable_h
#define ParameterTable_h 1

#include "parameter.h"

using namespace std;

class ParameterTable : public Parameter
{
private:
    friend class          ControlTable;
    friend class          DockingTable;
    friend class          FFTProcessTable;
#ifdef MPI_DP
    string            _RecPDB_file;
    string            _LigPDB_file;
    int               _IO_flag[3];

protected:
    virtual void          pdb_step();
#else
    string            _Table_file;
#endif

public:
    ParameterTable(Parallel *pparallel) : Parameter(pparallel) {
#ifdef DEBUG
        cout << "Constructing ParameterTable.\n";
#endif
    }
    virtual           ~ParameterTable() {
#ifdef DEBUG
        cout << "Destructing ParameterTable.\n";
#endif
    }
    virtual void          process_args(int argc,char *argv[]);
    using                 Parameter::initialize;
    virtual void          initialize(ParameterTable *pparameter);
    virtual void          output_file_name(const string rec_file, const string lig_file);
};

#endif
