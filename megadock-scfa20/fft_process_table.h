/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : FFTProcessTable
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef FFTProcessTable_h
#define FFTProcessTable_h 1

#include "exec_logger.h"
#include "fft_process.h"

class FFTProcessTable : public FFTProcess<ParameterTable>
{
private:
#ifdef CUFFT
    //host side
    cufftHandle cufft_plan;
    cufftResult cufft_result;
    cufftComplex *CUFFTin_host;
    cufftComplex *CUFFTout_host;
    float *top_score_host;
    int   *top_index_host;

    cudaStream_t _cuda_stream;

    //device side
    cufftComplex *CUFFTin_gpu;
    cufftComplex *CUFFTout_gpu;
    float *_FFT_rec_r_gpu;
    float *_FFT_rec_i_gpu;
    float *top_score_gpu;
    int   *top_index_gpu;

    float *radius_core2_gpu;
    float *radius_surf2_gpu;
    float *_Charge_gpu;
    float *xd_gpu;
    float *yd_gpu;
    float *zd_gpu;
    float *grid_coord_gpu;
    float *atom_coord_rotated_gpu;
    float *grid_r_gpu;
    float *grid_i_gpu;
    float *atom_coord_orig_gpu;
    float *mole_center_coord_gpu;
    float *ligand_rotation_angle_gpu;


#else

    fftwf_plan plan_fftw_forward;
    fftwf_plan plan_fftw_inverse;

    fftwf_complex         *_FFTWin; 
    fftwf_complex         *_FFTWout;

#endif

    ExecLogger           *_exec_logger;
    vector<SortScore>   _Select;
    int               _Current_rot_angle_num;
protected:
    virtual void              alloc_fft();
public:
    virtual void      fft3d(const float &theta);
    FFTProcessTable(ExecLogger *pexec_logger,Parallel *pparallel,ParameterTable *pparameter, Receptor<ParameterTable> *preceptor, Ligand<ParameterTable> *pligand)
        : _exec_logger(pexec_logger),FFTProcess(pparallel, pparameter, preceptor, pligand) {

#ifdef DEBUG
        cout << "Constructing FFTProcessTable.\n";
#endif
        //cout << "FFTP const "<< _parameter->_Num_sort <<endl; cout.flush();
    }
    virtual ~FFTProcessTable() {
#ifdef DEBUG
        cout << "Destructing FFTProcessTable.\n";
#endif

#ifdef CUFFT
        delete [] CUFFTin_host;
        delete [] CUFFTout_host;
#else

        fftwf_free(_FFTWin);
        fftwf_free(_FFTWout);
        fftwf_cleanup();

#endif

        vector<SortScore> tmp1;
        vector<SortScore> tmp2;
        _Select.swap(tmp1);
        _Top.swap(tmp2);
    }
    virtual void      alloc_array(const int &num_fft);
    virtual void      receptor_fft(float *grid_r,float *grid_i);
#ifdef CUFFT
    virtual void      cuda_fft(float *grid_r,float *grid_i,float *grid_coord,float *atom_coord_rotated,float *theta, size_t myid2);
    virtual void      ligand_voxelization_on_gpu(float *theta, size_t myid2);
    virtual void      ligand_data_transfer_gpu(float *grid_coord);
#else
    virtual void      ligand_preparation(float *grid_r,float *grid_i);
    virtual void      convolution();
    virtual void      score_sort();
#endif
    virtual void      fft_memory_free();
    virtual void      rotation_index(const int &j) {
        _Current_rot_angle_num = j;
    }
};

#endif
