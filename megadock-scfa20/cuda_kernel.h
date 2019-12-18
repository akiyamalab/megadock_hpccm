/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  cuda_kernel.cu
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef CudaKernel_h
#define CudaKernel_h 1

__global__ void lig_vox_fill(int ng1
                             ,int na
                             ,float delta
                             ,float *radius2
                             ,float *xd
                             ,float *yd
                             ,float *zd
                             ,float *grid_coord
                             ,float *atom_coord_rotated
                             ,float *grid_r
    						 ,float grid_width);



__global__ void lig_rotation(int na, float *theta, float *atom_coord_orig, float *mole_center_coord, float *atom_coord_rotated);

__global__ void ligvoxgpu_copy_htod(const int na, const float *const theta, const int ng1, const float *const atom_coord_orig, const float *const mole_center_coord, float *atom_coord_rotated, float *xd, float *yd, float *zd, const float *const grid_coord);

__global__ void lig_calc_dis_atomgrid(int na, int ng1, float *xd, float *yd, float *zd, float *grid_coord, float *atom_coord_rotated);


__global__ void lig_vox_init_grid(int ng3,float *grid_r,float *grid_i);


__global__ void lig_vox_init_fft(int nf3,cufftComplex *lig_in);


__global__ void lig_vox_init(int ng3,int nf3,float *grid_r,float *grid_i,cufftComplex *lig_in);

__global__ void ligand_voxel_set(int ng1
                                 ,cufftComplex *lig_in
                                 ,float *grid_r
                                 ,float *grid_i);



__global__ void lig_vox_surface_cut_CtoT(int ng1, float delta, float *grid_r);


__global__ void lig_vox_elec(int ng1,int na,float grid_width,float *_Charge,float *atom_coord_rotated,float *grid_i);


__global__ void lig_vox_elec_serial(int ng1,int na,float grid_width,float *_Charge,float *atom_coord_rotated,float *grid_i);



__device__ void lig_vox_surface_cut_TtoO(int ng3, float delta, float *grid_r);


__global__ void convolution_gpu(int nf3, float *rec_r, float *rec_i, cufftComplex *lig_out, cufftComplex *lig_in);


__global__ void max_pos_single(int nf3, cufftComplex *out, float *score, int *pos);


__global__ void max_pos_multi_set(int nf3, cufftComplex *out, float *temp_score, int *temp_index);


//, std::vector<cufftComplex> *temp_result , thrust::vector<cufftComplex> *temp_result
//thrust::device_ptr<cufftComplex> *temp_result cufftComplex *temp_result,thrust::device_ptr<cufftComplex> temp_result
__global__ void max_pos_multi(int nf3, cufftComplex *out, float *score, int *pos,const int num_sort,const int offset);




#endif
