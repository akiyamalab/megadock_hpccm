/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : FFTProcessPDB
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef FFTProcessPDB_h
#define FFTProcessPDB_h 1

#include "cpu_time.h"
#include "fft_process.h"

class FFTProcessPDB : public FFTProcess<ParameterPDB>
{
private:
#ifdef CUFFT
    //host side
    cufftHandle *cufft_plan;
    cufftResult *cufft_result;
    cufftComplex *CUFFTin_host;
    cufftComplex *CUFFTout_host;
    float **top_score_host;
    int   **top_index_host;

    //device side
    cufftComplex **CUFFTin_gpu;
    cufftComplex **CUFFTout_gpu;
    float **_FFT_rec_r_gpu;
    float **_FFT_rec_i_gpu;
    float **top_score_gpu;
    int   **top_index_gpu;

    float **radius_core2_gpu;
    float **radius_surf2_gpu;
    float **_Charge_gpu;
    float **xd_gpu;
    float **yd_gpu;
    float **zd_gpu;
    float **grid_coord_gpu;
    float **atom_coord_rotated_gpu;
    float **grid_r_gpu;
    float **grid_i_gpu;
    float **atom_coord_orig_gpu;
    float **mole_center_coord_gpu;
    float **ligand_rotation_angle_gpu;


#endif /* CUFFT */

    fftwf_plan *plan_fftw_forward;
    fftwf_plan *plan_fftw_inverse;

    fftwf_complex         *_FFTWin; 
    fftwf_complex         *_FFTWout;

    CPUTime           *_cputime;
    vector< vector<SortScore> >   _Select;
    int               *_Current_rot_angle_num;
protected:
    virtual void              alloc_fft();
public:
    virtual void          fft3d(const float &theta, size_t myid2);
    FFTProcessPDB(CPUTime *pcputime,Parallel *pparallel,ParameterPDB *pparameter, Receptor<ParameterPDB> *preceptor, Ligand<ParameterPDB> *pligand)
        : _cputime(pcputime),FFTProcess(pparallel, pparameter, preceptor, pligand) {

#ifdef DEBUG
        cout << "Constructing FFTProcessPDB.\n";
#endif
        //cout << "FFTP const "<< _parameter->_Num_sort <<endl; cout.flush();
    }
    virtual ~FFTProcessPDB() {
#ifdef DEBUG
        cout << "Destructing FFTProcessPDB.\n";
#endif
        delete [] _Current_rot_angle_num;

#ifdef CUFFT
        //host side
        delete [] cufft_plan;
        delete [] cufft_result;
        delete [] CUFFTin_host;
        delete [] CUFFTout_host;
        delete [] top_score_host;
        delete [] top_index_host;

        //device side
        delete [] CUFFTin_gpu;
        delete [] CUFFTout_gpu;
        delete [] _FFT_rec_r_gpu;
        delete [] _FFT_rec_i_gpu;
        delete [] top_score_gpu;
        delete [] top_index_gpu;

        delete [] radius_core2_gpu;
        delete [] radius_surf2_gpu;
        delete [] _Charge_gpu;
        delete [] xd_gpu;
        delete [] yd_gpu;
        delete [] zd_gpu;
        delete [] grid_coord_gpu;
        delete [] atom_coord_rotated_gpu;
        delete [] grid_r_gpu;
        delete [] grid_i_gpu;

        delete [] atom_coord_orig_gpu;
        delete [] mole_center_coord_gpu;
        delete [] ligand_rotation_angle_gpu;


#endif

        delete [] plan_fftw_forward;
        delete [] plan_fftw_inverse;

        fftwf_free(_FFTWin);
        fftwf_free(_FFTWout);
        fftwf_cleanup();

        vector< vector<SortScore> > tmp1;
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
#endif
    virtual void      ligand_preparation(float *grid_r,float *grid_i, size_t myid2);
    virtual void      convolution(size_t myid2);
    virtual void      score_sort(size_t myid2);
    virtual void      fft_memory_free();
    virtual void      rotation_index(int i,const int &j) {
        _Current_rot_angle_num[i] = j;
    }
};

#endif
