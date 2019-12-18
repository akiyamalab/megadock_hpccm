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

#include "fft_process.h"

#include <string.h>
#include <unistd.h>

#ifdef CUFFT

#include "cuda_kernel.cu"

#endif

//============================================================================//
bool ascending_struct(const SortScore &s1,const SortScore &s2)  // for sort
//============================================================================//
{
    return s1.score > s2.score;
}

//============================================================================//
template<class P> void FFTProcess<P>::sort_index(float *fwork,int *iwork)
//============================================================================//
{
    const int     no  = _parameter->_Num_output;
    const int     nt  = (_Top.size() < no) ? _Top.size() : no;

    partial_sort(_Top.begin(),_Top.begin()+nt,_Top.end(),
                 ascending_struct);

    return;
}

//============================================================================//
template<class P> void FFTProcess<P>::top_score_clean()
//============================================================================//
{
    const int     num_sort  = _parameter->_Num_sort;
    const int     num_angle  = _parameter->_Num_rot_angles;
    const int     no  = _parameter->_Num_output;
    int   num_toprank;

    num_toprank = num_angle * num_sort;
    if( no > num_toprank ) num_toprank = no;

    for( int j = 0 ; j < num_toprank ; j++ ) {
        _Top[j].score = 0.0;
    }

    return;
}

// explicit instantiations
template class FFTProcess<ParameterPDB>;
template class FFTProcess<ParameterTable>;
