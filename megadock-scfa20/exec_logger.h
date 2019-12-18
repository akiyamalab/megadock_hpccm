/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : ExecLogger
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef Exec_logger_h
#define Exec_logger_h 1

#include <stdio.h>
#include <string>
#include <iostream>

#include "cpu_time.h"

using namespace std;

class ExecLogger
{
private:
    ExecLogger(ExecLogger &c) {}
    const ExecLogger & operator=(const ExecLogger &c);

public:
    CPUTime *_cputime;
    long double mb;

#ifdef CUFFT
    size_t devmem_free, devmem_total, devmem_use;
#endif

    string _RLOut_file;
    int _Num_fft_flag;
    string rec_filename, lig_filename;
    float rec_max_size, lig_max_size;
    float rec_voxel_size, lig_voxel_size;
    int rec_num_grid, lig_num_grid;
    float grid_width;

    ExecLogger() {
#ifdef DEBUG
        cout << "Constructing ExecLogger.\n";
#endif
        _cputime = new CPUTime();
    }
    virtual   ~ExecLogger()      {
#ifdef DEBUG
        cout << "Destructing ExecLogger\n";
#endif
        delete _cputime;
    }
    virtual void      initialize();
    virtual void      output(const int myid2);
    virtual void      record_malloc(const int &size) {
        _cputime->record_malloc(size);
    }
    virtual void      record_free(const int &size) {
        _cputime->record_free(size);
    }
};

#endif
