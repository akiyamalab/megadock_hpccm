/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : Control
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef Control_h
#define Control_h 1

#include "parallel.h"
#include "receptor.h"
#include "ligand.h"
#include "docking_pdb.h"
#include "docking_table.h"

using namespace std;

template<class P, class D> class Control
{
private:
    Control(Control &c) {}
    const Control & operator=(const Control &c);
protected:
    Parallel  *_parallel;
    P *_parameter;
    Receptor<P>  *_receptor;
    Ligand<P>    *_ligand;
    D            *_docking;
    virtual void  gridtable_11base_normal(int &ngrid,vector<int> &ngrid_table);
    virtual void  gridtable_13base_normal(int &ngrid,vector<int> &ngrid_table);
    virtual void  gridtable_07base_normal(int &ngrid,vector<int> &ngrid_table);
    virtual void  gridtable_fftw_custom(int &ngrid,vector<int> &ngrid_table);
    virtual void  gridtable_cufft_custom(int &ngrid,vector<int> &ngrid_table);
    virtual void  autogridr(const int &ngrid,vector<int> &ngrid_table) = 0;
    virtual void  autogridl(const int &ngrid,vector<int> &ngrid_table) = 0;
    virtual void  checkgridr() = 0;
    virtual void  checkgridl() = 0;
public:
    Control(Parallel *pparallel) : _parallel(pparallel) {
#ifdef DEBUG
        cout << "Constructing Control.\n";
#endif
    }
    Control(Parallel *pparallel, P *pparameter) : _parallel(pparallel), _parameter(pparameter) {
#ifdef DEBUG
        cout << "Constructing Control.\n";
#endif
    }
    virtual ~Control() {
#ifdef DEBUG
        cout << "Destructing Control.\n";
#endif
    }
    virtual void  execute() = 0;
};

#endif
