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

#include "application.h"

//============================================================================//
void Application::initialize()
//============================================================================//
{
#ifdef CUFFT
    checkCudaErrors( cudaGetDeviceCount(&device_count_gpu) );
    if (device_count_gpu == 0) {
        fprintf(stderr, "GPU Error: no devices supporting CUDA.\n");
        exit(-1);
    }

    cudaDeviceProp deviceProp;
    checkCudaErrors( cudaGetDeviceProperties(&deviceProp, 0));
    if (deviceProp.major < 1) {
        fprintf(stderr, "GPU Error: device does not support CUDA.\n");
        exit(-1);
    }

    cudaSetDeviceFlags(cudaDeviceMapHost);
    fprintf(stdout, "# Using CUDA device %d: %s\n", 0, deviceProp.name);
    cudaSetDevice(0);
    //fprintf(stdout, "# Init CUDA device OK.\n");

    int cufft_version;
    cufftGetVersion(&cufft_version);
    printf("# CUFFT version : %d\n", cufft_version);

    printf("# Number of available [threads / GPUs] : [%d / %d]\n",nproc2,device_count_gpu);
#endif

    _parallels = new Parallel*[nproc2];
    _exec_loggers = new ExecLogger*[nproc2];
    _controls = new ControlTable*[nproc2];
    _parameters = new ParameterTable*[nproc2];

    for (int i = 0; i < nproc2; i++) {
        _parallels[i] = new Parallel(nproc2);
        _parallels[i]->num_gpu(device_count_gpu);
        _exec_loggers[i] = new ExecLogger();

        // ParameterTable
        _parameters[i] = new ParameterTable(_parallels[i]);
        if (i == 0) {
            _parameters[i]->initialize();
        } else {
            _parameters[i]->initialize(_parameters[0]);
        }
        _exec_loggers[i]->record_malloc(_parameters[i]->allocate_size()); //Rotation angles[], Atom radius, charge, ACE[]

        _controls[i] = new ControlTable(_exec_loggers[i],_parallels[i],_parameters[i]);
        _controls[i]->initialize(i == 0);
    }
}

//============================================================================//
int Application::application(int argc, char *argv[], int myid2)
//============================================================================//
{
    struct timeval et1, et2;
    gettimeofday(&et1,NULL);
    _exec_loggers[myid2]->initialize();
#pragma omp critical (prepare)
    {
        _parameters[myid2]->process_args(argc, argv);
        _controls[myid2]->prepare();
    }
    _controls[myid2]->execute();

    gettimeofday(&et2,NULL);

    const float elapsed_time = (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
    printf("\n");

#pragma omp critical (output)
    {
        printf("# ========================================\n");
        _exec_loggers[myid2]->output(myid2);
        printf("Elapsed time                  = %8.2f sec.\n"
               "# ========================================\n"
               ,elapsed_time);
    }
    return 0;
}
