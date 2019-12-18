/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : Docking
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef Docking_h
#define Docking_h 1

#include "parallel.h"
#include "receptor.h"
#include "ligand.h"

using namespace std;

template<class P, class F> class Docking
{
private:
    Docking(Docking &c) {}
    const Docking & operator=(const Docking &c);
protected:
    Parallel  *_parallel;
    P *_parameter;
    Receptor<P>  *_receptor;
    Ligand<P>    *_ligand;
    F    *_fft_process;
    int       _Num_grid;
    float     *_Grid_coord;
    float     *_Fwork;
    int       *_Iwork;
    size_t       _Memfw;
    size_t       _Memiw;
    virtual void  maxsize_voxel() = 0;
    virtual void  alloc_array(const int &maxatom, const int &nag, const size_t &ng3) = 0;
public:
    Docking(Parallel *pparallel,P *pparameter,
            Receptor<P> *rreceptor,Ligand<P> *rligand)
        : _parallel(pparallel),_parameter(pparameter),
          _receptor(rreceptor),_ligand(rligand) {
#ifdef DEBUG
        cout << "Constructing Docking.\n";
#endif
    }
    virtual ~Docking() {
#ifdef DEBUG
        cout << "Destructing Docking.\n";
#endif
        delete [] _Grid_coord;
        delete [] _Fwork;
        delete [] _Iwork;
        delete _fft_process;
    }
    virtual void  initialize() = 0;
    virtual void  rec_init() = 0;
    virtual void  dockz() = 0;
    virtual void  dock_memory_free() = 0;
    virtual void  output() = 0;
    virtual void  output_detail() = 0; // for analysis
    virtual void  output_calc_time_log() = 0; // for analysis
};

#endif
