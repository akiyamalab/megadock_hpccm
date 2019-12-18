/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : (main)
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#include <string.h>
#include "cpu_time.h"
#include "exec_logger.h"
#include "control_pdb.h"
#include "control_table.h"

#ifdef CUFFT
#include <helper_cuda.h>
#define VERSION "4.1.0 for GPU & "
#else
#define VERSION "4.1.0 for CPU & "
#endif

#ifdef MPI_DP
#define VTEXT "multiple nodes"
#else
#define VTEXT "single node"
#endif

#define LASTUPDATED "26 March, 2019"

struct DockingPair {
    string rec_file, lig_file, out_file;
    DockingPair(string rec_file, string lig_file, string out_file) : rec_file(rec_file), lig_file(lig_file), out_file(out_file) {}
};

//============================================================================//
void get_pair(string line, string &rec_file, string &lig_file, string &out_file)
//============================================================================//
{
    int first_tab_index = line.find_first_of('\t');
    if (first_tab_index == string::npos) {
        cerr << "[Error] Ligand is not specified." << endl;
        exit(1);
    }
    rec_file = line.substr(0, first_tab_index);

    int second_tab_index = line.find_last_of('\t');
    if (first_tab_index == second_tab_index) {
        lig_file = line.substr(first_tab_index + 1, line.size() - 1 - first_tab_index);
        out_file = "";
    } else {
        lig_file = line.substr(first_tab_index + 1, second_tab_index - 1 - first_tab_index);
        out_file = line.substr(second_tab_index + 1, line.size() - 1 - second_tab_index);
    }
}

//============================================================================//
void initialize(int argc, char *argv[], int &nproc2, int &device_count_gpu)
//============================================================================//
{
    cout << " MEGADOCK ver. "<< VERSION << VTEXT <<  endl;
    cout << "      megadock@bi.c.titech.ac.jp   lastupdated: " << LASTUPDATED << endl;
    cout << endl;

#ifdef _OPENMP
    #pragma omp parallel
    {
        nproc2 = omp_get_num_threads();
        if(omp_get_thread_num() == 0) {
            cout << "# Using OpenMP parallelization: " << nproc2 << " threads." << endl;
        }
    }
    //printf("#OpenMP version %d\n", _OPENMP);
#else
    nproc2 = 1;
#endif //#ifdef _OPENMP

#ifdef CUFFT
    int nogpu_flag = 0;
    for (int num = 0; num < (argc-1); ++num) {
        if(!strncmp(argv[num], "-G", 2)) {
            if(argv[num+1] != NULL) {
                if(atoi(argv[num+1]) == 0) {
                    nogpu_flag = 1;
                }
            }
        }
    }

    if(nogpu_flag != 1) {
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
    }

    printf("# Number of available [threads / GPUs] : [%d / %d]\n",nproc2,device_count_gpu);
#endif
}

//============================================================================//
void main_pdb(int argc, char *argv[])
//============================================================================//
{
    Parallel  *_parallel;
    CPUTime   *_cputime;
    ControlPDB   *_control;

    struct timeval et1, et2;
    struct timeval et3, et4;
    int nproc2 = 0;
    int device_count_gpu = 0;

    gettimeofday(&et1,NULL);
    gettimeofday(&et3,NULL);

    initialize(argc, argv, nproc2, device_count_gpu);

    _cputime = new CPUTime();
    _cputime->initialize();

    _parallel = new Parallel(nproc2);
    _parallel->num_gpu(device_count_gpu); 

    gettimeofday(&et4,NULL);
    _cputime->t1_initialize += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));

    _control = new ControlPDB(_cputime,_parallel);
    _control->initialize(argc,argv);
    _control->execute();

    delete _control;
    delete _parallel;

    _cputime->output();

    delete _cputime;

    gettimeofday(&et2,NULL);

    const float elapsed_time = (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
    printf("\n");
    printf("Elapsed time                  = %8.2f sec.\n",elapsed_time);
}

//============================================================================//
void main_table(int argc, char *argv[])
//============================================================================//
{
    struct timeval et3, et4;
    int nproc2 = 0;
    int device_count_gpu = 0;

    gettimeofday(&et3,NULL);

    initialize(argc, argv, nproc2, device_count_gpu);

    struct timeval et1[nproc2], et2[nproc2];

    Parallel  *_parallels[nproc2];
    ExecLogger   *_exec_loggers[nproc2];
    ControlTable   *_controls[nproc2];
    ParameterTable *_parameters[nproc2];

    for (int i = 0; i < nproc2; i++) {
        _parallels[i] = new Parallel(nproc2);
        _parallels[i]->num_gpu(device_count_gpu);
        _exec_loggers[i] = new ExecLogger();

        // ParameterTable
        _parameters[i] = new ParameterTable(_parallels[i]);
        if (i == 0) {
            _parameters[i]->initialize(argc,argv);
        } else {
            _parameters[i]->initialize(_parameters[0]);
        }
        _exec_loggers[i]->record_malloc(_parameters[i]->allocate_size()); //Rotation angles[], Atom radius, charge, ACE[]

        _controls[i] = new ControlTable(_exec_loggers[i],_parallels[i],_parameters[i]);
        _controls[i]->initialize(i == 0);
    }


    ifstream input_stream(_controls[0]->input_file());
    if (!input_stream.is_open()) {
        cerr << "Unable to open input file." << endl;
        exit(1);
    }
    string line;
    vector<DockingPair> pairs;
    while (getline(input_stream, line)) {
        string rec_file, lig_file, out_file;
        get_pair(line, rec_file, lig_file, out_file);
        pairs.push_back(DockingPair(rec_file, lig_file, out_file));
    }

#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < pairs.size(); i++) {
        int myid2 = omp_get_thread_num();
        DockingPair _docking_pair = pairs[i];
        gettimeofday(&et1[myid2],NULL);
        _exec_loggers[myid2]->initialize();
#pragma omp critical (prepare)
        {
            _controls[myid2]->prepare(_docking_pair.rec_file, _docking_pair.lig_file, _docking_pair.out_file);
        }
        _controls[myid2]->execute();

        gettimeofday(&et2[myid2],NULL);

        const float elapsed_time = (et2[myid2].tv_sec-et1[myid2].tv_sec + (float)((et2[myid2].tv_usec-et1[myid2].tv_usec)*1e-6));
        printf("\n");

#pragma omp critical (output)
        {
            printf("# ========================================\n");
            _exec_loggers[myid2]->output(myid2);
            printf("Elapsed time                  = %8.2f sec.\n"
                   "# ========================================\n"
                   ,elapsed_time);
        }
    }

#pragma omp parallel for
    for (int i = 0; i < nproc2; i++) {
        delete _exec_loggers[i];
        delete _controls[i];
        delete _parallels[i];
        delete _parameters[i];
    }

    gettimeofday(&et4,NULL);

    const float total_time = (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));
    printf("\n");
    printf("Total time                    = %8.2f sec.\n",total_time);
}

//============================================================================//
#ifdef MPI_DP
int application(int argc,char *argv[])
#else
int main(int argc, char *argv[])
#endif
//============================================================================//
{
    bool table_input_flag = false, pdb_input_flag = false;
    for (int num = 0; num < argc; ++num) {
        if (!(strncmp(argv[num], "-R", 2) && strncmp(argv[num], "-L", 2) && strncmp(argv[num], "-o", 2))) {
            pdb_input_flag = true;
        } else if (!strncmp(argv[num], "-I", 2)) {
            table_input_flag = true;
        } else if (!strncmp(argv[num], "-h", 2)) {
            usage();
        }
    }
    if (pdb_input_flag) {
        if (table_input_flag) {
            fprintf(stderr, "[ERROR] A pair of PDB files and a docking pair list file cannot be specified simultaneously.\n");
            usage();
        } else {
            main_pdb(argc, argv);
        }
    } else {
        if (table_input_flag) {
            main_table(argc, argv);
        } else {
            fprintf(stderr, "[ERROR] A pair of PDB files or a docking pair list file has to be specified.\n");
            usage();
        }
    }
    return 0;
}
