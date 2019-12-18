/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : Application
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef Application_h
#define Application_h 1

#include "cpu_time.h"
#include "exec_logger.h"
#include "control_pdb.h"
#include "control_table.h"

class Application
{
private:
    Application(Application &c) {}
    const Application & operator=(const Application &c);

    int nproc2;
    int device_count_gpu;
    Parallel  **_parallels;
    ExecLogger   **_exec_loggers;
    ControlTable   **_controls;
    ParameterTable **_parameters;

public:
    Application(const int nproc2) : nproc2(nproc2) {}
    virtual ~Application() {
#pragma omp parallel for
        for (int i = 0; i < nproc2; i++) {
            delete _exec_loggers[i];
            delete _controls[i];
            delete _parallels[i];
            delete _parameters[i];
        }

        delete [] _parallels;
        delete [] _exec_loggers;
        delete [] _controls;
        delete [] _parameters;
    }
    virtual void initialize();
    virtual int application(int argc, char *argv[], int myid2);
};

#endif
