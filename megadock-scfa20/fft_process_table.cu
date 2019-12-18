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

#include "fft_process_table.h"

#define NUM_THREADS 512 //should be power of 2

#ifdef CUFFT

#include "cuda_kernel.h"

#endif

//============================================================================//
void FFTProcessTable::alloc_array(const int &num_fft)
//============================================================================//
{
    //cout << "FFT::alloc_array |" <<num_fft<< endl; cout.flush();
    _Num_fft = num_fft;

    const size_t nf3       = _Num_fft * _Num_fft * _Num_fft;
    const int num_sort  = _parameter->_Num_sort;
    const int num_angle = _parameter->_Num_rot_angles;
    const int no        = _parameter->_Num_output;
    const size_t nproc2    = _parallel->nproc2();
    int   num_toprank;

    num_toprank = num_angle * num_sort;
    if( no > num_toprank ) num_toprank = no;

    alloc_fft();

    _Select.resize(num_sort);

    _Top.resize(num_toprank);

    //---------- memory allocation for _Current_rot_angle_num
    //_Current_rot_angle_num = new int[nproc2];

    _exec_logger->record_malloc( sizeof(float)*nf3*2*(1 + nproc2));

    //---------- memory allocation for _FFT_rec_r
    _FFT_rec_r = new float[nf3];
    if( !_FFT_rec_r ) {
        cerr << "[ERROR] Out of memory. Number of listed receptors = ("
             << nf3 << ") for (_FFT_rec_r) in fft_process.cpp!!\n";
        exit(1);
    }

    //---------- memory allocation for _FFT_rec_i
    _FFT_rec_i = new float[nf3];
    if( !_FFT_rec_i ) {
        cerr << "[ERROR] Out of memory. Number of listed receptors = ("
             << nf3 << ") for (_FFT_rec_i) in fft_process.cpp!!\n";
        exit(1);
    }

    return;
}

//============================================================================//
void FFTProcessTable::alloc_fft()
//============================================================================//
{
    const int nf1 = _Num_fft;
    const size_t nf3 = _Num_fft * _Num_fft * _Num_fft;
    const size_t nproc2  = _parallel->nproc2();
    const int num_gpu = _parallel->num_gpu();
    const int na = _ligand->num_atoms();

#ifdef CUFFT
    const int num_sort = _parameter->_Num_sort;
    const int ng1 = _Num_fft / 2;
    const int ng3 = ng1 * ng1 * ng1;
    const int nag = na * ng1;
    //for ligand voxelization on GPU
    const int nThreads = NUM_THREADS;
    const int nBlocks_nf3 = (nf3 + (nThreads-1)) / nThreads;

    CUFFTin_host  = new cufftComplex[nf3];
    CUFFTout_host = new cufftComplex[nf3];

    _exec_logger->record_malloc( sizeof(cufftComplex)*nf3*2 ); //_in/outBuf

    //printf(" start: %p\n",&CUFFTin_host[0].x);

    int lenCUFFTin_host = (int)(((long int)&CUFFTin_host[nf3-1].x) - ((long int)&CUFFTin_host[0].x) + sizeof(CUFFTin_host[nf3-1]))/sizeof(CUFFTin_host[nf3-1]);
    if(lenCUFFTin_host !=nf3) printf("# discontinuous memory allocation occurs\n");

    //printf("   end: %ld\n",(long long int)&CUFFTin_host[nf3-1].y - &CUFFTin_host[0].x);

    int myid2 = omp_get_thread_num();
    cudaSetDevice(myid2 % num_gpu);
    checkCudaErrors( cudaStreamCreate(&_cuda_stream));
    cufft_result = cufftPlan3d(&cufft_plan, nf1, nf1, nf1, CUFFT_C2C);
    cufftSetStream(cufft_plan, _cuda_stream);

    checkCudaErrors( cudaMalloc((void **)&CUFFTin_gpu,  sizeof(cufftComplex)*nf3) );
    checkCudaErrors( cudaMalloc((void **)&CUFFTout_gpu, sizeof(cufftComplex)*nf3) );
    checkCudaErrors( cudaMalloc((void **)&_FFT_rec_r_gpu, sizeof(float)*nf3) );
    checkCudaErrors( cudaMalloc((void **)&_FFT_rec_i_gpu, sizeof(float)*nf3) );

    checkCudaErrors( cudaMalloc((void **)&grid_r_gpu,  sizeof(float)*ng3));
    checkCudaErrors( cudaMalloc((void **)&grid_i_gpu,  sizeof(float)*ng3));
    checkCudaErrors( cudaMalloc((void **)&grid_coord_gpu,  sizeof(float)*ng1));
    checkCudaErrors( cudaMalloc((void **)&radius_core2_gpu,  sizeof(float)*na));
    checkCudaErrors( cudaMalloc((void **)&radius_surf2_gpu,  sizeof(float)*na));
    checkCudaErrors( cudaMalloc((void **)&_Charge_gpu,  sizeof(float)*na));
    checkCudaErrors( cudaMalloc((void **)&xd_gpu,  sizeof(float)*nag));
    checkCudaErrors( cudaMalloc((void **)&yd_gpu,  sizeof(float)*nag));
    checkCudaErrors( cudaMalloc((void **)&zd_gpu,  sizeof(float)*nag));
    checkCudaErrors( cudaMalloc((void **)&atom_coord_rotated_gpu,  sizeof(float)*na*3));
    checkCudaErrors( cudaMalloc((void **)&atom_coord_orig_gpu,  sizeof(float)*na*3));
    checkCudaErrors( cudaMalloc((void **)&mole_center_coord_gpu,  sizeof(float)*3));
    checkCudaErrors( cudaMalloc((void **)&ligand_rotation_angle_gpu,  sizeof(float)*3));
    checkCudaErrors( cudaMalloc((void **)&top_score_gpu, sizeof(float)*nBlocks_nf3*num_sort) );
    checkCudaErrors( cudaMalloc((void **)&top_index_gpu, sizeof(int)*nBlocks_nf3*num_sort) );
    top_score_host = new float[nBlocks_nf3];
    top_index_host = new int[nBlocks_nf3];

    _exec_logger->record_malloc( sizeof(float)*nBlocks_nf3 + sizeof(int)*nBlocks_nf3 );

    cudaMemGetInfo(&(_exec_logger->devmem_free), &(_exec_logger->devmem_total));
    _exec_logger->devmem_use = _exec_logger->devmem_total - _exec_logger->devmem_free;

#else

    _FFTWin  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*nf3);
    _FFTWout = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*nf3);

    plan_fftw_forward=fftwf_plan_dft_3d(nf1,nf1,nf1,_FFTWin,_FFTWout,FFTW_FORWARD,FFTW_ESTIMATE);
    plan_fftw_inverse=fftwf_plan_dft_3d(nf1,nf1,nf1,_FFTWin,_FFTWout,FFTW_BACKWARD,FFTW_ESTIMATE);

    _exec_logger->record_malloc( sizeof(fftwf_complex)*nf3*2 );

#endif
    return;
}

//============================================================================//
void FFTProcessTable::receptor_fft(float *grid_r,float *grid_i)
//============================================================================//
{
    const int num_grid= _Num_fft / 2;
    const size_t nf3 = _Num_fft * _Num_fft * _Num_fft;
    const int ndata   = ( _Num_fft - num_grid ) / 2;
    const float   theta   = -2.0 * PI / _Num_fft;

    const int num_gpu = _parallel->num_gpu();
    const int nproc2 = _parallel->nproc2();

    if(num_gpu > 0) {
#ifdef CUFFT
        int myid2 = omp_get_thread_num();
        struct timeval et1, et2;
        //memset(CUFFTin_host[0], make_cuComplex(0.0, 0.0), sizeof(cufftComplex)*nf3);
        for( int i = 0 ; i < nf3 ; i++ ) {
            CUFFTin_host[i] = make_cuComplex(0.0, 0.0);
        }

        for( int i = 0, m = 0 ; i < num_grid ; i++ ) {
            const int ic = _Num_fft*_Num_fft*(i+ndata);
            for( int j = 0 ; j < num_grid ; j++ ) {
                const int jc = ic + _Num_fft*(j+ndata);
                for( int k = 0 ; k < num_grid ; k++ ) {
                    CUFFTin_host[jc+k+ndata] = make_cuComplex(grid_r[m  ], grid_i[m]);
                    m++;
                }
            }
        }

        cudaSetDevice(myid2 % num_gpu); //CUFFTin_dev[0] : [0] means 0th GPU

        gettimeofday(&et1,NULL);
        checkCudaErrors( cudaMemcpyAsync(CUFFTin_gpu, CUFFTin_host, sizeof(cufftComplex)*nf3, cudaMemcpyHostToDevice, _cuda_stream) );
        gettimeofday(&et2,NULL);
        _exec_logger->_cputime->t6_data_transfer_rec += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

        fft3d(theta); // [0] means performed on 0th GPU

        gettimeofday(&et1,NULL);
        checkCudaErrors( cudaMemcpyAsync(CUFFTout_host,CUFFTout_gpu,sizeof(cufftComplex)*nf3,cudaMemcpyDeviceToHost, _cuda_stream) );
        gettimeofday(&et2,NULL);
        _exec_logger->_cputime->t6_data_transfer_rec += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

        for( int i = 0 ; i < nf3 ; i++ ) {
            _FFT_rec_r[i] = cuCrealf(CUFFTout_host[i]);
            _FFT_rec_i[i] = cuCimagf(CUFFTout_host[i]);
        }

        gettimeofday(&et1,NULL);

        checkCudaErrors( cudaMemcpyAsync(_FFT_rec_r_gpu, _FFT_rec_r, sizeof(float)*nf3, cudaMemcpyHostToDevice, _cuda_stream) );
        checkCudaErrors( cudaMemcpyAsync(_FFT_rec_i_gpu, _FFT_rec_i, sizeof(float)*nf3, cudaMemcpyHostToDevice, _cuda_stream) );

        gettimeofday(&et2,NULL);
        _exec_logger->_cputime->t6_data_transfer_rec += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
#endif
    } else {
#ifndef CUFFT
        memset(_FFTWin, 0.0, sizeof(fftwf_complex)*nf3);

        for( int i = 0, m = 0 ; i < num_grid ; i++ ) {
            const int ic = _Num_fft*_Num_fft*(i+ndata);

            for( int j = 0 ; j < num_grid ; j++ ) {
                const int jc = ic + _Num_fft*(j+ndata);

                for( int k = 0 ; k < num_grid ; k++ ) {
                    _FFTWin[jc+k+ndata][0] = grid_r[m  ];
                    _FFTWin[jc+k+ndata][1] = grid_i[m++];
                }
            }
        }

        fft3d(theta);

        for( int i = 0 ; i < nf3 ; i++ ) {
            _FFT_rec_r[i] = _FFTWout[i][0];
            _FFT_rec_i[i] = _FFTWout[i][1];
        }
#endif
    }


    return;
}

#ifndef CUFFT
//============================================================================//
void FFTProcessTable::ligand_preparation(float *grid_r,float *grid_i)
//============================================================================//
{
    const int ng1 = _Num_fft / 2;
    const int nf2 = _Num_fft * _Num_fft;
    const size_t nf3 = _Num_fft * _Num_fft * _Num_fft;
    const int ndata   = ( _Num_fft - ng1 ) / 2;
   
    memset(_FFTWin[0], 0.0, sizeof(fftwf_complex)*nf3);
        
    for( int i = 0, m = 0 ; i < ng1 ; i++ ) {
        const int ic = nf2*(i+ndata);

        for( int j = 0 ; j < ng1 ; j++ ) {
            int jc = ic + _Num_fft*(j+ndata);
            
            for( size_t k = 0, myijk=jc+ndata ; k < ng1 ; k++, myijk++ ) {
                _FFTWin[myijk][0] = grid_r[m  ];
                _FFTWin[myijk][1] = grid_i[m++];
            }
        }
    }
    
    return;
}

//============================================================================//
void FFTProcessTable::convolution()
//============================================================================//
{
    const int nf1 = _Num_fft;
    const int nf2 = nf1*nf1;
    const size_t nf3 = nf1*nf2;

    for( size_t i = 0, j=0 ; i < nf3 ; i++,j++ ) {
      _FFTWin[j][0] = _FFT_rec_r[i]*_FFTWout[j][0] + _FFT_rec_i[i]*_FFTWout[j][1];
      _FFTWin[j][1] = _FFT_rec_r[i]*_FFTWout[j][1] - _FFT_rec_i[i]*_FFTWout[j][0];
    }

    return;
}
#endif

//============================================================================//
void FFTProcessTable::fft3d(const float &theta)
//============================================================================//
{   
    const size_t nproc2  = _parallel->nproc2();
    const int num_gpu = _parallel->num_gpu();

#ifdef CUFFT
    const int nf1 = _Num_fft;
    cufftHandle plan;
    cufftResult res;

    res = cufftPlan3d(&plan, nf1, nf1, nf1, CUFFT_C2C);
    cufftSetStream(plan, _cuda_stream);
    if(!res == CUFFT_SUCCESS) {
        cout << "!fail to plan 3d FFT (DFT):" << res << endl;
        exit(-1);
    }

    if( theta < 0.0 ) {
        res = cufftExecC2C(plan, CUFFTin_gpu, CUFFTout_gpu, CUFFT_FORWARD);
    } else {
        res = cufftExecC2C(plan, CUFFTin_gpu, CUFFTout_gpu, CUFFT_INVERSE);
    }

    if(!res == CUFFT_SUCCESS) {
        cout << "!fail to exec 3d FFT(in fft3d()):" << res << endl;
        exit(-1);
    }

    res =  cufftDestroy(plan);
#else
    struct timeval et3, et4;
    gettimeofday(&et3,NULL);
    if( _parameter->fft_library_type == 2 ) {        
    } else {
        if( theta < 0.0 ) {
            fftwf_execute(plan_fftw_forward);
        } else {
            fftwf_execute(plan_fftw_inverse);
        }
    }
    gettimeofday(&et4,NULL);
    //printf(" [FFT(host),%s] %10.5f\n\n",((theta<0.0)?"Forward":"Inverse"),(et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6)));
#endif

    return;
}

#ifndef CUFFT
//============================================================================//
void FFTProcessTable::score_sort()
//============================================================================//
{
    const int num_sort  = _parameter->_Num_sort;
    const int nf2 = _Num_fft * _Num_fft;
    const int nf3 = _Num_fft * _Num_fft * _Num_fft;
    float temp_top_score;
    int temp_top_index;

    for( int i = 0 ; i < num_sort ; i++ ) {
        _Select[i].score = -99999.0;
    }

    fftwf_complex *fftout;
    fftout = _FFTWout;
    
    if(num_sort!=1) {
        for( size_t i = 0,myi= 0 ; i < nf3 ; i++,myi++ ) {
            const float raw = fftout[myi][0] / nf3;
            if( raw < _Select[num_sort-1].score) continue;
            for( int j = 0 ; j < num_sort ; j++ ) {
                if( raw > _Select[j].score ) {
                    for( int k = num_sort-1 ; k > j ; k-- ) {
                        _Select[k] = _Select[k-1];
                    }
                    _Select[j].score    = raw;
                    _Select[j].index[1] = i / nf2;
                    _Select[j].index[2] = (i / _Num_fft) % _Num_fft;
                    _Select[j].index[3] = i % _Num_fft;
                    break;
                }
            }
        }
    } else { // num_sort = 1, take only 1 score per angle
        temp_top_score = 0.0;
        temp_top_index = 0;
        for( size_t i = 0, myi=0 ; i < nf3 ; i++,myi++ ) {
            const float raw = fftout[myi][0];
            if (temp_top_score < raw) {
                temp_top_score = raw;
                temp_top_index = i;
            }
        }
        _Select[0].score    = temp_top_score / nf3;
        _Select[0].index[1] = temp_top_index / nf2;
        _Select[0].index[2] = (temp_top_index / _Num_fft) % _Num_fft;
        _Select[0].index[3] = temp_top_index % _Num_fft;
    }

    for( int i = 0 ; i < num_sort ; i++ ) {
        //printf(" top %d %f\n",i,_Select[i].score);
        _Select[i].index[0] = _Current_rot_angle_num;
    }

    for( int i = 0 ; i < num_sort ; i++ ) {
        _Top[_Current_rot_angle_num*num_sort+i] = _Select[i];
    }

    return;
}
#endif

#ifdef CUFFT
//============================================================================//
void FFTProcessTable::cuda_fft(float *grid_r,float *grid_i,float *grid_coord,float *atom_coord_rotated,float *theta, size_t myid2)
//============================================================================//
{
    const int nf1 = _Num_fft;
    const int nf2 = nf1 * nf1;
    const size_t nf3 = nf2 * nf1;
    const int num_gpu = _parallel->num_gpu();
    const size_t nproc2    = _parallel->nproc2();

    const int num_sort = _parameter->_Num_sort;
    const int na = _ligand->num_atoms();

    struct timeval et1, et2;
    struct timeval et3, et4;
    gettimeofday(&et1,NULL);

    float temp_top_score = -999999.0;
    int temp_top_index = -999999;

    const int nThreads = NUM_THREADS;
    const int nBlocks_nf3 = (nf3 + (nThreads-1)) / nThreads;
    if(nBlocks_nf3 * nThreads < nf3) {
        printf(" nf3:%d, nBlocks_nf3:%d, nThreads:%d , nf3=nBlocks_nf3*nThreads\n",nf3,nBlocks_nf3,nThreads);
        fprintf(stderr, " [ERROR] too large FFT size. nf3:%d, nBlocks_nf3:%d\n", nf3, nBlocks_nf3);
        exit(1);
    }

    cudaSetDevice(myid2 % num_gpu);
    //printf(" #p10 [myid=%d]\n",myid2);

    ligand_voxelization_on_gpu(theta,myid2);
    checkCudaErrors( cudaStreamSynchronize(_cuda_stream) );

    gettimeofday(&et2,NULL);
    _exec_logger->_cputime->t3_1_ligand_voxelization += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
    gettimeofday(&et1,NULL);

    cufft_result = cufftExecC2C(cufft_plan, CUFFTin_gpu, CUFFTout_gpu, CUFFT_FORWARD);
    if(!cufft_result == CUFFT_SUCCESS) {
        cout << "!fail to exec 3d FFT (DFT, Lig):" << cufft_result << endl;
        exit(-1);
    }

    //*/
    checkCudaErrors( cudaStreamSynchronize(_cuda_stream) );

    gettimeofday(&et2,NULL);
    _exec_logger->_cputime->t3_2_fftprocess_ligand_fft += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

    gettimeofday(&et1,NULL);
    convolution_gpu<<<nBlocks_nf3, nThreads, 0, _cuda_stream>>>(nf3, _FFT_rec_r_gpu, _FFT_rec_i_gpu, CUFFTout_gpu, CUFFTin_gpu);

    checkCudaErrors( cudaStreamSynchronize(_cuda_stream) );

    gettimeofday(&et2,NULL);
    _exec_logger->_cputime->t3_3_fftprocess_convolution += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
    gettimeofday(&et1,NULL);

    cufft_result = cufftExecC2C(cufft_plan, CUFFTin_gpu, CUFFTout_gpu, CUFFT_INVERSE);
    if(!(cufft_result == CUFFT_SUCCESS)) {
        cout << "!fail to exec 3d FFT (IDFT):" << cufft_result << endl;
        exit(-1);
    }
    //*
    checkCudaErrors( cudaStreamSynchronize(_cuda_stream) );
    gettimeofday(&et2,NULL);
    _exec_logger->_cputime->t3_4_fftprocess_fft_inverse += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
    gettimeofday(&et1,NULL);

    // Search max score translation position from CUFFTout_gpu[nf3]

    //printf(" t=%d per angle\n",num_sort);

    for( int i = 0 ; i < num_sort ; i++ ) {
        _Select[i].score = -99999.0;
    }

    max_pos_single<<<nBlocks_nf3, nThreads, sizeof(float)*nThreads, _cuda_stream>>>(nf3, CUFFTout_gpu,  top_score_gpu, top_index_gpu);
    checkCudaErrors( cudaStreamSynchronize(_cuda_stream) );

    gettimeofday(&et3,NULL);
    checkCudaErrors( cudaMemcpyAsync(top_score_host,top_score_gpu,sizeof(float)*nBlocks_nf3,cudaMemcpyDeviceToHost, _cuda_stream) );
    checkCudaErrors( cudaMemcpyAsync(top_index_host,top_index_gpu,sizeof(int)*nBlocks_nf3,cudaMemcpyDeviceToHost, _cuda_stream) );
    gettimeofday(&et4,NULL);
    _exec_logger->_cputime->t6_data_transfer_in_loop += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));
    checkCudaErrors( cudaStreamSynchronize(_cuda_stream) );

    if(num_sort!=1) {
        for(int i=0; i<nBlocks_nf3; i++) {
            if(top_index_host[i]/nf2 > nf1 || top_index_host[i] < 0){
                top_score_host[i] = -99999.99;
                //printf(" error, %d | score, %f \n", top_index_host[i]/nf2, top_score_host[i]);
            }
            const float raw = top_score_host[i];
            if( raw < _Select[num_sort-1].score) continue;
            for( int j = 0 ; j < num_sort ; j++ ) {
                if( raw > _Select[j].score ) {
                    for( int k = num_sort-1 ; k > j ; k-- ) {
                        _Select[k] = _Select[k-1];
                    }
                    const int index = top_index_host[i];
                    _Select[j].score    = raw;
                    _Select[j].index[1] = index / nf2;
                    _Select[j].index[2] = (index / _Num_fft) % _Num_fft;
                    _Select[j].index[3] = index % _Num_fft;
                    break;
                }
            }
        }

    } else { // num_sort = 1, select only 1 score per 1 ligand angle
        for(int i=0; i<nBlocks_nf3; i++) {
            if(top_index_host[i]/nf2 > nf1 || top_index_host[i] < 0){
                top_score_host[i] = -99999.99;
                //printf(" error, %d | score, %f \n", top_index_host[i]/nf2, top_score_host[i]);
            }
            if(temp_top_score < top_score_host[i]) {
                temp_top_score = top_score_host[i];
                temp_top_index = top_index_host[i];
            }
        }

        //printf("  m:%f\n\n",temp_top_score);
        //printf("%g (%d) [%d %d %d]\n", temp_top_score, _p, temp_top_index/(n*n),(temp_top_index/n)%n, temp_top_index%n );
        //printf("<%d> %g (%d/%d) %d\n", nBlocks,temp_top_score, temp_top_index, nf3, temp_top_index/nf2);

        _Select[0].score    = temp_top_score;
        _Select[0].index[1] = temp_top_index / nf2;
        _Select[0].index[2] = (temp_top_index / nf1) % nf1;
        _Select[0].index[3] = temp_top_index % nf1;
        /* / DEBUG
        printf("TEST,  %d\n", _Select[0].index[1]);
        if ( _Select[0].index[1] > nf1 ){
            printf(" error, %d\n", _Select[0].index[1]);
            }*/

    }

    //*** score_sort ***********************************************************

    for( int i = 0 ; i < num_sort ; i++ ) {
        _Select[i].index[0] = _Current_rot_angle_num;
        _Top[_Current_rot_angle_num*num_sort+i] = _Select[i];
    }

    //size_t devmem_use, devmem_free, devmem_total;
    //cudaMemGetInfo(&devmem_free, &devmem_total);
    //devmem_use = devmem_total - devmem_free;
    //printf(" [GPU (%d) memory] Use : %10u (%4.1f%%), Free : %10u (%4.1f%%), Total : %10u\n",myid2,devmem_use,(float)(100*devmem_use/devmem_total), devmem_free, (float)(100*devmem_free/devmem_total), devmem_total);


    checkCudaErrors( cudaStreamSynchronize(_cuda_stream) );
    gettimeofday(&et2,NULL);
    _exec_logger->_cputime->t3_5_fftprocess_score_sort += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

    return;
}


//============================================================================//
void FFTProcessTable::ligand_voxelization_on_gpu(float *theta, size_t myid2)
//============================================================================//
{
    const int ng1 = _Num_fft / 2;
    const int ng3 = ng1 * ng1 * ng1;
    const int nf1 = _Num_fft;
    const int nf2 = nf1 * nf1;
    const size_t nf3 = nf2 * nf1;

    const float delta = 1.0;
    const float surface = 1.0;
    const float grid_width = _parameter->grid_width;
    const int sr_half = (2.4 + grid_width - 0.01) / grid_width;
    const int sr = 2 * sr_half + 1;

    const int na = _ligand->num_atoms();
    const int nag = na * ng1;
    const int na_sr3 = na * sr * sr * sr;

    struct timeval et1, et2;
    struct timeval et3, et4;

    const int nThreads = NUM_THREADS;
    //const int nBlocks_na = (na + (nThreads-1)) / nThreads;
    const int nBlocks_nag = (nag + (nThreads-1)) / nThreads;
    const int nBlocks_na_sr3 = (na_sr3 + (nThreads-1)) / nThreads;
    const int nBlocks_ng3 = (ng3 + (nThreads-1)) / nThreads;
    const int nBlocks_nf3 = (nf3 + (nThreads-1)) / nThreads;
    if(nBlocks_nf3 * nThreads < nf3) {
        printf(" nf3:%d, nBlocks_nf3:%d, nThreads:%d , nf3=nBlocks_nf3*nThreads\n",nf3,nBlocks_nf3,nThreads);
        fprintf(stderr, " [ERROR] too large FFT size. nf3:%d, nBlocks_nf3:%d\n", nf3, nBlocks_nf3);
        exit(1);
    }

    //*
    //transfer ligand angle & calc xd,yd,zd,atom_coord_rotated
    gettimeofday(&et3,NULL);

    gettimeofday(&et1,NULL);
    checkCudaErrors( cudaMemcpyAsync(ligand_rotation_angle_gpu, theta, sizeof(float)*3, cudaMemcpyHostToDevice, _cuda_stream) );
    gettimeofday(&et2,NULL);
    _exec_logger->_cputime->t3_1_ligand_voxelization += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

    //lig_rotation<<<nBlocks_na, nThreads, 0, _cuda_stream>>>(na, ligand_rotation_angle_gpu,atom_coord_orig_gpu, mole_center_coord_gpu, atom_coord_rotated_gpu);
    //checkCudaErrors( cudaStreamSynchronize(_cuda_stream) );
    //lig_calc_dis_atomgrid<<<nBlocks_nag, nThreads, 0, _cuda_stream>>>(na, ng1, xd_gpu, yd_gpu, zd_gpu, grid_coord_gpu, atom_coord_rotated_gpu);
    //checkCudaErrors( cudaStreamSynchronize(_cuda_stream) );
    ligvoxgpu_copy_htod<<<nBlocks_nag, nThreads, 0, _cuda_stream>>>
        (na, ligand_rotation_angle_gpu, ng1, atom_coord_orig_gpu, mole_center_coord_gpu, atom_coord_rotated_gpu, xd_gpu, yd_gpu, zd_gpu, grid_coord_gpu);
    checkCudaErrors( cudaStreamSynchronize(_cuda_stream) );
    gettimeofday(&et4,NULL);
    _exec_logger->_cputime->t3_1_1_ligvoxgpu_copy_htod += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));

    //grid[] initialize
    gettimeofday(&et3,NULL);
    lig_vox_init<<<nBlocks_nf3, nThreads, 0, _cuda_stream>>>(ng3,nf3,grid_r_gpu,grid_i_gpu,CUFFTin_gpu);
    //lig_vox_init_fft<<<nBlocks_nf3, nThreads, 0, _cuda_stream>>>(nf3,CUFFTin_gpu);
    checkCudaErrors( cudaStreamSynchronize(_cuda_stream) );
    gettimeofday(&et4,NULL);
    _exec_logger->_cputime->t3_1_2_ligvoxgpu_kernel_init += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));

    //atom fill(core)
    gettimeofday(&et3,NULL);
    lig_vox_fill<<<nBlocks_na_sr3, nThreads, 0, _cuda_stream>>>
    (ng1,na,delta,radius_core2_gpu,xd_gpu,yd_gpu,zd_gpu,grid_coord_gpu,atom_coord_rotated_gpu,grid_r_gpu, grid_width);
    checkCudaErrors( cudaStreamSynchronize(_cuda_stream) );
    gettimeofday(&et4,NULL);
    _exec_logger->_cputime->t3_1_3_ligvoxgpu_kernel_fill_core += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));

    //surface cutting
    gettimeofday(&et3,NULL);
    lig_vox_surface_cut_CtoT<<<nBlocks_ng3, nThreads, 0, _cuda_stream>>>(ng1,delta,grid_r_gpu);
    checkCudaErrors( cudaStreamSynchronize(_cuda_stream) );
    gettimeofday(&et4,NULL);
    _exec_logger->_cputime->t3_1_4_ligvoxgpu_kernel_cut_surf += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));

    //atom fill(surf)
    gettimeofday(&et3,NULL);
    lig_vox_fill<<<nBlocks_na_sr3, nThreads, 0, _cuda_stream>>>
    (ng1,na,surface,radius_surf2_gpu,xd_gpu,yd_gpu,zd_gpu,grid_coord_gpu,atom_coord_rotated_gpu,grid_r_gpu, grid_width);
    checkCudaErrors( cudaStreamSynchronize(_cuda_stream) );
    gettimeofday(&et4,NULL);
    _exec_logger->_cputime->t3_1_5_ligvoxgpu_kernel_fill_surf += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));

    //electro
    gettimeofday(&et3,NULL);

    if(_parameter->lig_elec_serial_flag == 0) {
        lig_vox_elec<<<nBlocks_ng3, nThreads, 0, _cuda_stream>>>(ng1, na, grid_width, _Charge_gpu, atom_coord_rotated_gpu, grid_i_gpu);
    } else {
        lig_vox_elec_serial<<<nBlocks_ng3, nThreads, 0, _cuda_stream>>>(ng1, na, grid_width, _Charge_gpu, atom_coord_rotated_gpu, grid_i_gpu);
    }

    checkCudaErrors( cudaStreamSynchronize(_cuda_stream) );
    gettimeofday(&et4,NULL);
    _exec_logger->_cputime->t3_1_6_ligvoxgpu_kernel_elec += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));

    //set Voxel grid[ng3] into center of FFT grid[nf3]
    gettimeofday(&et3,NULL);
    ligand_voxel_set<<<nBlocks_ng3, nThreads, 0, _cuda_stream>>>(ng1,CUFFTin_gpu,grid_r_gpu,grid_i_gpu);
    checkCudaErrors( cudaStreamSynchronize(_cuda_stream) );
    gettimeofday(&et4,NULL);
    _exec_logger->_cputime->t3_1_7_ligvoxgpu_kernel_set_array += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));

}


//============================================================================//
void FFTProcessTable::ligand_data_transfer_gpu(float *grid_coord)
//============================================================================//
{
    const int ng1 = _Num_fft / 2;
    const int na = _ligand->num_atoms();
    const int num_gpu = _parallel->num_gpu();
    const int nproc2 = _parallel->nproc2();
    const float   rcore2 = 1.5;           // ZDOCK parameter
    const float   rsurf2 = 1.0;           // ZDOCK parameter
    struct timeval et1, et2;

    float radius_core2[na];
    float radius_surf2[na];

    for(int i = 0; i < na; i++) {
        radius_core2[i] = _ligand->_Radius[i] * _ligand->_Radius[i] * rcore2;
        radius_surf2[i] = _ligand->_Radius[i] * _ligand->_Radius[i] * rsurf2;
    }

    gettimeofday(&et1,NULL);
    int myid2 = omp_get_thread_num();
    cudaSetDevice(myid2 % num_gpu);
    checkCudaErrors( cudaMemcpyAsync(radius_core2_gpu, radius_core2, sizeof(float)*na, cudaMemcpyHostToDevice, _cuda_stream) );
    checkCudaErrors( cudaMemcpyAsync(radius_surf2_gpu, radius_surf2, sizeof(float)*na, cudaMemcpyHostToDevice, _cuda_stream) );
    checkCudaErrors( cudaMemcpyAsync(_Charge_gpu, _ligand->_Charge, sizeof(float)*na, cudaMemcpyHostToDevice, _cuda_stream) );
    checkCudaErrors( cudaMemcpyAsync(grid_coord_gpu, grid_coord, sizeof(float)*ng1, cudaMemcpyHostToDevice, _cuda_stream) );
    checkCudaErrors( cudaMemcpyAsync(atom_coord_orig_gpu, _ligand->_Coordinate, sizeof(float)*na*3, cudaMemcpyHostToDevice, _cuda_stream) );
    checkCudaErrors( cudaMemcpyAsync(mole_center_coord_gpu, _ligand->_Center, sizeof(float)*3, cudaMemcpyHostToDevice, _cuda_stream) );

    gettimeofday(&et2,NULL);
    _exec_logger->_cputime->t6_data_transfer_lig += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

    return;
}

#endif /* CUFFT */



//============================================================================//
void FFTProcessTable::fft_memory_free()
//============================================================================//
{
    const size_t nproc2    = _parallel->nproc2();
    const int num_gpu = _parallel->num_gpu();
    const size_t nf3 = _Num_fft * _Num_fft * _Num_fft;

#ifndef CUFFT
    fftwf_destroy_plan(plan_fftw_forward);
    fftwf_destroy_plan(plan_fftw_inverse);

    _exec_logger->record_free(sizeof(float)*nf3*2);

#else

    //const int num_sort = _parameter->_Num_sort;
    const int nThreads = NUM_THREADS;
    const int nBlocks_nf3 = (nf3 + (nThreads-1)) / nThreads;

    int myid2 = omp_get_thread_num();
    cudaSetDevice(myid2 % num_gpu);

    cufftDestroy(cufft_plan);

    checkCudaErrors( cudaStreamDestroy(_cuda_stream));

    checkCudaErrors( cudaFree(CUFFTin_gpu));
    checkCudaErrors( cudaFree(CUFFTout_gpu));
    checkCudaErrors( cudaFree(_FFT_rec_r_gpu));
    checkCudaErrors( cudaFree(_FFT_rec_i_gpu));


    checkCudaErrors( cudaFree(grid_r_gpu));
    checkCudaErrors( cudaFree(grid_i_gpu));
    checkCudaErrors( cudaFree(grid_coord_gpu));

    checkCudaErrors( cudaFree(radius_core2_gpu));
    checkCudaErrors( cudaFree(radius_surf2_gpu));
    checkCudaErrors( cudaFree(_Charge_gpu));

    checkCudaErrors( cudaFree(xd_gpu));
    checkCudaErrors( cudaFree(yd_gpu));

    checkCudaErrors( cudaFree(zd_gpu));

    checkCudaErrors( cudaFree(atom_coord_rotated_gpu));
    checkCudaErrors( cudaFree(atom_coord_orig_gpu));
    checkCudaErrors( cudaFree(mole_center_coord_gpu));
    checkCudaErrors( cudaFree(ligand_rotation_angle_gpu));

    checkCudaErrors( cudaFree(top_score_gpu));
    checkCudaErrors( cudaFree(top_index_gpu));

    delete [] top_score_host;
    delete [] top_index_host;


    _exec_logger->record_free( sizeof(float)*nBlocks_nf3 + sizeof(int)*nBlocks_nf3 );

#endif

    return;
}
