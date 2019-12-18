/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : DockingPDB
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef DockingPDB_h
#define DockingPDB_h 1

#include "docking.h"
#include "cpu_time.h"
#include "parameter_pdb.h"
#include "fft_process_pdb.h"

using namespace std;

class DockingPDB : public Docking<ParameterPDB, FFTProcessPDB>
{
private:
    CPUTime   *_cputime;
    float     **_Mol_coord;
protected:
    virtual void  maxsize_voxel();
    virtual void  alloc_array(const int &maxatom, const int &nag, const size_t &ng3);
    virtual void  create_voxel(Protein<ParameterPDB> *rprotein, size_t myid2);
    virtual void  ligand_rotationz(float *theta, size_t myid2);
public:
    DockingPDB(CPUTime *pcputime,Parallel *pparallel,ParameterPDB *pparameter,
            Receptor<ParameterPDB> *rreceptor,Ligand<ParameterPDB> *rligand)
        : _cputime(pcputime),Docking(pparallel,pparameter,rreceptor,rligand) {
#ifdef DEBUG
        cout << "Constructing DockingPDB.\n";
#endif
    }
    virtual ~DockingPDB() {
#ifdef DEBUG
        cout << "Destructing DockingPDB.\n";
#endif
        delete [] _Mol_coord;
    }
    virtual void  initialize();
    virtual void  rec_init();
    virtual void  dockz();
    virtual void  dock_memory_free();
    virtual void  output();
    virtual void  output_detail(); // for analysis
    virtual void  output_calc_time_log(); // for analysis
};

#endif
