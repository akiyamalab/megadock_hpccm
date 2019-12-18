/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : FFTProcess
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef FFTProcess_h
#define FFTProcess_h 1

#include <algorithm>
#include <functional>
#include <cstring>
#include <unistd.h>


#include "parallel.h"
#include "parameter_table.h"
#include "parameter_pdb.h"
#include "protein.h"
#include "ligand.h"
#include "receptor.h"

#ifdef CUFFT
#include "helper_cuda.h"
#include "cufft.h"
//#include <thrust/device_vector.h>
//#include <thrust/device_ptr.h>

#endif

#include "fftw3.h"

using namespace std;

typedef struct {
    float score;
    int   index[4];
} SortScore;

template<class P> class FFTProcess
{
private:
    FFTProcess(FFTProcess &c) {}
    const FFTProcess & operator=(const FFTProcess &c);

protected:
    Parallel          *_parallel;
    P                 *_parameter;
    Ligand<P>         *_ligand;
    Receptor<P>       *_receptor;
    vector<SortScore>     _Top;
    vector< vector<SortScore> >   _Select;
    int               _Num_fft;
    float             *_FFT_rec_r;
    float             *_FFT_rec_i;
    virtual void              alloc_fft() = 0;
public:
    FFTProcess(Parallel *pparallel,P *pparameter, Receptor<P> *preceptor, Ligand<P> *pligand)
        : _parallel(pparallel),_parameter(pparameter),_receptor(preceptor),_ligand(pligand) {

#ifdef DEBUG
        cout << "Constructing FFTProcess.\n";
#endif
        //cout << "FFTP const "<< _parameter->_Num_sort <<endl; cout.flush();
    }
    virtual ~FFTProcess() {
#ifdef DEBUG
        cout << "Destructing FFTProcess.\n";
#endif
        delete [] _FFT_rec_r;
        delete [] _FFT_rec_i;

        vector< vector<SortScore> > tmp1;
        vector<SortScore> tmp2;
        _Select.swap(tmp1);
        _Top.swap(tmp2);
    }
    virtual void      alloc_array(const int &num_fft) = 0;
    virtual void      receptor_fft(float *grid_r,float *grid_i) = 0;
#ifdef CUFFT
    virtual void      cuda_fft(float *grid_r,float *grid_i,float *grid_coord,float *atom_coord_rotated,float *theta, size_t myid2) = 0;
    virtual void      ligand_voxelization_on_gpu(float *theta, size_t myid2) = 0;
    virtual void      ligand_data_transfer_gpu(float *grid_coord) = 0;
#endif
    virtual void      fft_memory_free() = 0;
    virtual void      top_score_clean();
    virtual int       num_fft() {
        return _Num_fft;
    }
    virtual void      num_fft(const int &i) {
        _Num_fft = i;
    }
    virtual float     top_score(const int &j) {
        return _Top[j].score;
    }
    virtual int       top_index(const int &j,int k) {
        return _Top[j].index[k];
    }
    virtual void      sort_index(float *fwork,int *iwork);
};

#endif
