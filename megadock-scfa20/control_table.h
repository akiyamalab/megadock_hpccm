/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : ControlTable
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef ControlTable_h
#define ControlTable_h 1

#include "control.h"
#include "exec_logger.h"
#include "parameter_table.h"
#include "docking_table.h"

using namespace std;

class ControlTable : public Control<ParameterTable, DockingTable>
{
private:
    ExecLogger   *_exec_logger;
protected:
    virtual void  autogridr(const int &ngrid,vector<int> &ngrid_table);
    virtual void  autogridl(const int &ngrid,vector<int> &ngrid_table);
    virtual void  checkgridr();
    virtual void  checkgridl();
public:
    ControlTable(ExecLogger *pexec_logger,Parallel *pparallel,ParameterTable *pparameter)
        : _exec_logger(pexec_logger),Control<ParameterTable, DockingTable>(pparallel,pparameter) {
#ifdef DEBUG
        cout << "Constructing ControlTable.\n";
#endif
    }
    virtual ~ControlTable() {
#ifdef DEBUG
        cout << "Destructing ControlTable.\n";
#endif
    }
    virtual void  initialize(bool verbose);
#ifdef MPI_DP
    virtual void  prepare();
#else
    virtual void  prepare(string rec_file, string lig_file, string out_file);
#endif
    virtual void  execute();
#ifndef MPI_DP
    virtual string input_file() {
        return _parameter->_Table_file;
    }
#endif
};

#endif
