/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : DockingTable
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#include "docking_table.h"

//============================================================================//
void DockingTable::initialize()
//============================================================================//
{
    const float grid_width = _parameter->grid_width;

    _Num_grid = max(_receptor->num_grid(),_ligand->num_grid());
    _parameter->_Num_grid = _Num_grid;
    _parameter->_Num_fft = _Num_grid * 2;

    _fft_process  = new FFTProcessTable(_exec_logger,_parallel,_parameter,_receptor,_ligand);

    maxsize_voxel();

    return;
}

//============================================================================//
void DockingTable::maxsize_voxel()
//============================================================================//
{
    const int num_fft = _Num_grid * 2;
    const size_t natom = max(_receptor->num_atoms(),_ligand->num_atoms());

    _parameter->_Num_atom_max = natom;

    const size_t nag = natom * _Num_grid;
    const size_t ng3 = _Num_grid * _Num_grid * _Num_grid;
    const size_t nf3 = num_fft * num_fft * num_fft;

    const int num_sort  = _parameter->_Num_sort;
    const int num_angle  = _parameter->_Num_rot_angles;
    const size_t nproc2  = _parallel->nproc2();
    long double fmem;
    long double imem;
    const int MB = 1048576; // 1024*1024

    fmem  = nproc2*(ng3*3+natom*3+nag*3) / MB;
    fmem += nproc2*2*nf3 / MB;
    fmem += 2*nf3 / MB;
    fmem += nproc2*4*nf3 / MB;
    fmem += num_angle*num_sort / MB;

    imem  = nproc2*(ng3*2+natom*4) / MB;
    imem += 4*num_angle*num_sort / MB;

    long double mb = fmem*4+imem*2;
    mb += 10.0;       // margine

    _fft_process->alloc_array(num_fft);
    alloc_array(natom,nag,ng3);

    return;
}

//============================================================================//
void DockingTable::alloc_array(const int &natom,const int &nag,const size_t &ng3)
//============================================================================//
{
    const size_t nproc2    = _parallel->nproc2();

    _Grid_coord = new float[_Num_grid];
    if( !_Grid_coord ) {
        cerr << "[ERROR] Out of memory. Number of voxels for 1 axis = ("
            << _Num_grid
            << ") for (_Grid_coord) in docking.cpp!!\n";
        exit(1);
    }

    _Mol_coord = new float[natom*3];
    if( !_Mol_coord ) {
        cerr << "[ERROR] Out of memory. Number of atoms times 3 axis = ("
            << natom*3
            << ") for (_Mol_coord) in docking.cpp!!\n";
        exit(1);
    }

    _Memfw = ng3*3+natom*3+nag*3;
    _Fwork = new float[_Memfw];
    if( !_Fwork ) {
        cerr << "[ERROR] Out of memory. Number of float arrays = ("
            << _Memfw*nproc2
            << ") for (_Fwork) in docking.cpp!!\n";
        exit(1);
    }

    _Memiw = ng3*2+natom*4;
    _Iwork = new int[_Memiw*nproc2];
    if( !_Iwork ) {
        cerr << "[ERROR] Out of memory. Number of int arrays = ("
            << _Memiw*nproc2
            << ") for (_Iwork) in docking.cpp!!\n";
        exit(1);
    }

    _exec_logger->record_malloc( sizeof(float)*natom*3 ); //_Mol_coord
    _exec_logger->record_malloc( sizeof(float)*_Memfw*nproc2 + sizeof(int)*_Memiw*nproc2 ); //_F/Iwork

    return;
}

//============================================================================//
void DockingTable::rec_init()
//============================================================================//
{
    //const int ng1 = 4;
    const int ng1 = _Num_grid;
    const int ng2 = ng1 * ng1;
    const int ng3 = ng2 * ng1;
    float search_length   = _parameter->grid_width * _Num_grid;
    struct timeval et1, et2;
    gettimeofday(&et1,NULL);

    for( int i = 0 ; i < _Num_grid ; i++ ) {
        _Grid_coord[i] = search_length * ( -0.5 + ( 0.5 + i )/_Num_grid );
        //printf(" grid coord[%d], %f\n",i,_Grid_coord[i]);
    }

    _receptor->shift_center(_Mol_coord);
    create_voxel(_receptor);


    gettimeofday(&et2,NULL);
    //printf(" [Rec-1] %10.5f\n",(et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6)));
    gettimeofday(&et1,NULL);

    _fft_process->receptor_fft(_Fwork,&_Fwork[ng3]);

    gettimeofday(&et2,NULL);
    //printf(" [Rec-2] %10.5f\n",(et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6)));
    gettimeofday(&et1,NULL);

    return;
}

//============================================================================//
void DockingTable::create_voxel(Protein<ParameterTable> *rprotein)
//============================================================================//
{
    struct timeval et1, et2;
    if( _parameter->_Score_func == 3 ) {   // rPSC + CHARMM19 electrostatic + Receptor ACE
        const float beta = -2800.0;                 // MEGADOCK parameter
        gettimeofday(&et1,NULL);
        rprotein->rpscace(_parameter->_ACE_ratio, _Num_grid, _Grid_coord, _Mol_coord,_Iwork,
                _Fwork, _parameter->_rPSC_param_rec_core, _parameter->_rPSC_param_lig_core
                , _parameter->_Old_voxel_flag);
        gettimeofday(&et2,NULL);
        //printf(" vox rpsc %10.5f\n", (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6)));
        gettimeofday(&et1,NULL);
        rprotein->electro(beta,_parameter->_Elec_ratio,_Num_grid,_Grid_coord,_Mol_coord,
                _Iwork,_Fwork, _parameter->_Old_voxel_flag);
        gettimeofday(&et2,NULL);
        //printf(" vox elec %10.5f\n", (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6)));
    }

    return;
}

//============================================================================//
void DockingTable::dockz()
//============================================================================//
{
    struct timeval et1, et2;
    struct timeval et3, et4;
    float theta[3];
    int xyz,myid2=0;
    const size_t nproc2  = _parallel->nproc2();
    const int num_gpu = _parallel->num_gpu();

    const int ng1 = _Num_grid;
    const int ng2 = ng1 * ng1;
    const int ng3 = ng2 * ng1;

    //const int nc = max(_parameter->_Num_rot_angles / 10, 1);

    _fft_process->top_score_clean();

    int loop_count;
    float each_thread_time;
    struct timeval time_thread[2];

    //printf(" [%2d Threads]\n",nproc2);

#ifdef CUFFT
    if(num_gpu > 0) {
        gettimeofday(&et1,NULL);
        _fft_process->ligand_data_transfer_gpu(_Grid_coord); //transfer ligand infomation(atom radius, charge, original coord, grid coord)
        gettimeofday(&et2,NULL);
        _exec_logger->_cputime->t3_1_1_ligvoxgpu_copy_htod += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
    }
#endif //CUFFT

    myid2 = omp_get_thread_num();
    for( int ang = 0 ; ang < _parameter->_Num_rot_angles ; ang++ ) {
        gettimeofday(&time_thread[0],NULL);
        loop_count++;

        _fft_process->rotation_index(ang); // _Current_rot_angle_num[myid2] = ang;

        for(xyz = 0 ; xyz < 3 ; xyz++ ) {
            theta[xyz] = _parameter->_Zangle[3*ang+xyz];
        }

#ifdef CUFFT
        _fft_process->cuda_fft(_Fwork, _Fwork,_Grid_coord,_Mol_coord,theta,myid2);
#else

        gettimeofday(&et1,NULL); // ligand_voxelization --------------------

        ligand_rotationz(theta);
        create_voxel(_ligand);

        gettimeofday(&et2,NULL);
        _exec_logger->_cputime->t3_1_ligand_voxelization += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

        gettimeofday(&et1,NULL);
        _fft_process->ligand_preparation(_Fwork,_Fwork);
        gettimeofday(&et2,NULL);
        _exec_logger->_cputime->t3_1_ligand_voxelization += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

        gettimeofday(&et1,NULL); // ligand_fft -----------------------------
        _fft_process->fft3d(-1.0);
        gettimeofday(&et2,NULL);
        _exec_logger->_cputime->t3_2_fftprocess_ligand_fft += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

        gettimeofday(&et1,NULL); // convolution ----------------------------
        _fft_process->convolution();
        gettimeofday(&et2,NULL);
        //comout for separation of convolution and IDFT
        _exec_logger->_cputime->t3_3_fftprocess_convolution += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

        gettimeofday(&et1,NULL); // fft_inverse ----------------------------
        _fft_process->fft3d(1.0);
        gettimeofday(&et2,NULL);
        //comout for separation of convolution and IDFT
        _exec_logger->_cputime->t3_4_fftprocess_fft_inverse += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

        gettimeofday(&et1,NULL); // score_sort -----------------------------
        _fft_process->score_sort();
        gettimeofday(&et2,NULL);
        _exec_logger->_cputime->t3_5_fftprocess_score_sort += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
#endif



        gettimeofday(&time_thread[1],NULL);
        each_thread_time += (time_thread[1].tv_sec-time_thread[0].tv_sec + (float)((time_thread[1].tv_usec-time_thread[0].tv_usec)*1e-6));
        //printf("%f\n",(time_thread[1].tv_sec-time_thread[0].tv_sec + (float)((time_thread[1].tv_usec-time_thread[0].tv_usec)*1e-6)));
    }

    gettimeofday(&et1,NULL);
    _fft_process->sort_index(_Fwork,_Iwork);
    gettimeofday(&et2,NULL);
    _exec_logger->_cputime->t3_6_fftprocess_sort_index += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

    const float average_angle = (float)_parameter->_Num_rot_angles;
    float stddev_angle = 0.0;
    stddev_angle = pow((float)(loop_count-average_angle),2);

    stddev_angle = (float)(sqrt(stddev_angle));

    if(_parameter->tem_flag2==1) {
        printf("\n Load balance (Ave : %5.1f angles , Stddev : %5.1f angles)\n\n",average_angle,stddev_angle);
        printf(" Thread %2d : %5d angles, (%4.0f%%), %10.2f sec\n",myid2,loop_count,(float)(loop_count*100.0/average_angle),each_thread_time);
        printf("%d,",loop_count);
        printf("%f,",stddev_angle);
        printf("\n");
    }
    return;
}

//============================================================================//
void DockingTable::dock_memory_free()
//============================================================================//
{
    const size_t nproc2    = _parallel->nproc2();
    const int natom = _parameter->_Num_atom_max;

    delete [] _Mol_coord;

    _exec_logger->record_free( sizeof(float)*natom*3 ); //_Mol_coord

    _fft_process->fft_memory_free();

    return;
}


//============================================================================//
void DockingTable::ligand_rotationz(float *theta)
//============================================================================//
{
    const int nl = _ligand->num_atoms();
    float     x, y, z;

    //theta[0]: rotation around Z axis, theta[1]: X axis, theta[2]: Z axis

    const float r11 = cos(theta[0])*cos(theta[2])  -  sin(theta[0])*cos(theta[1])*sin(theta[2]);
    const float r21 = sin(theta[0])*cos(theta[2])  +  cos(theta[0])*cos(theta[1])*sin(theta[2]);
    const float r31 = sin(theta[1])*sin(theta[2]);
    const float r12 = -cos(theta[0])*sin(theta[2])  -  sin(theta[0])*cos(theta[1])*cos(theta[2]);
    const float r22 = -sin(theta[0])*sin(theta[2])  +  cos(theta[0])*cos(theta[1])*cos(theta[2]);
    const float r32 = sin(theta[1])*cos(theta[2]);
    const float r13 = sin(theta[0])*sin(theta[1]);
    const float r23 = -cos(theta[0])*sin(theta[1]);
    const float r33 = cos(theta[1]);

    for( int i = 0 ; i < nl ; i++ ) {
        x = _ligand->coordinate(i,0) - _ligand->_Center[0];
        y = _ligand->coordinate(i,1) - _ligand->_Center[1];
        z = _ligand->coordinate(i,2) - _ligand->_Center[2];
        _Mol_coord[3*i  ] = r11 * x + r12 * y + r13 * z;
        _Mol_coord[3*i+1] = r21 * x + r22 * y + r23 * z;
        _Mol_coord[3*i+2] = r31 * x + r32 * y + r33 * z;
    }

    return;
}

//============================================================================//
void DockingTable::output()
//============================================================================//
{
    char      buff[256];
    float     coord[3], coordl[3], theta[3];

    const int nf1 = _Num_grid * 2;
    const float   grid_width     = _parameter->grid_width;
    const int num_output  = _parameter->_Num_output;

    _exec_logger->_RLOut_file = _parameter->_RLOut_file;

    FILE *Output = fopen(_parameter->_RLOut_file.c_str(), "w");
    if( !Output ) {
        cerr << "[ERROR] Output file [" << _parameter->_RLOut_file << "] not open!!" << endl;
        exit(1);
    }

    fprintf(Output,"%d\t%.2f\n",nf1,grid_width); // Line 1 : Number of FFT, Grid width

    for( int xyz = 0 ; xyz < 3 ; xyz++ ) {
        theta[xyz] = _receptor->angle(xyz);
    }

    fprintf(Output,"\t%f\t%f\t%f\n",theta[0],theta[1],theta[2]); // Line 2 : angle of Receptor

    for( int xyz = 0 ; xyz < 3 ; xyz++ ) {
        coord[xyz] = _receptor->center(xyz);
        coordl[xyz] = _ligand->center(xyz);
    }

    fprintf(Output,"%s\t%f\t%f\t%f\n",_receptor->input_file().c_str(),coord[0],coord[1],coord[2]); // Line 3 : Receptor PDB file, center coordinate

    fprintf(Output,"%s\t%f\t%f\t%f\n",_ligand->input_file().c_str(),coordl[0],coordl[1],coordl[2]); // Line 4 : Ligand PDB file, center coordinate

    for( int rank = 0 ; rank < num_output ; rank++ ) {
        const int n = _fft_process->top_index(rank,0);
        for( int k = 0 ; k < 3 ; k++ ) {
            theta[k] = _parameter->_Zangle[3*n+k];
            if( theta[k] >= PI ) {
                theta[k] -= PI*2.0;
            }
        }
        fprintf(Output,"%f\t%f\t%f\t%d\t%d\t%d\t%.2f\n",theta[0],theta[1],theta[2], // Line 5- : angle, coord, score
                _fft_process->top_index(rank,1),_fft_process->top_index(rank,2),
                _fft_process->top_index(rank,3),_fft_process->top_score(rank));
    }

    fclose(Output);

    return;
}

//============================================================================//
void DockingTable::output_detail()
//============================================================================//
{

    // ##############################
    // #Output detailed result      #
    // ##############################

    float score_rpsc, score_elec, score_ace, score_total;
    float score_rpsc_plus, score_rpsc_minus;
    float score_elec_plus, score_elec_minus;
    float score_ace_plus, score_ace_minus;
    int count_rpsc_plus, count_rpsc_minus;
    int count_elec_plus, count_elec_minus;
    int count_ace_plus, count_ace_minus;
    float theta[3], coordr[3], coordl[3], diff_dis[3];
    float *rec_r, *rec_i, *rec_r_rpsc, *lig_r, *lig_i;
    float *lig_atom_coord_orig, *lig_atom_coord_temp;
    int la, lb, ma, mb, na, nb;
    int cor_coord[3], tem_coord[3];
    float rmsd_lca,rmsd_lbb;
    int natom_lca = 0,natom_lbb = 0;
    int near_native;
    float nn_thr = 5.0;
    char buff[256];

    const int ng1 = _Num_grid;
    const int ng2 = ng1 * ng1;
    const int ng3 = ng2 * ng1;
    const int nf1 = ng1 * 2;
    const int natom = _ligand->num_atoms();
    const int na3 = natom * 3;
    const float rmsd2_thr =(float) pow(10.0,2.0) * natom;
    //const int nf2 = nf1 * nf1;
    //const int nf3 = nf2 * nf1;

    const float grid_width = _parameter->grid_width;
    const int num_output = _parameter->_Num_output;
    //const int orig_ACE_ratio = _parameter->_ACE_ratio;

    rec_r = new float[ng3];
    rec_r_rpsc = new float[ng3];
    rec_i = new float[ng3];
    lig_r = new float[ng3];
    lig_i = new float[ng3];
    lig_atom_coord_orig = new float[na3];
    lig_atom_coord_temp = new float[na3];


    FILE *Output_detail = fopen(_parameter->_RLOut_file_detail.c_str(), "w");
    if( !Output_detail ) {
        cerr << "[ERROR] Output_detail file [" << _parameter->_RLOut_file_detail << "] not open!!" << endl;
        exit(1);
    }
    FILE *Output_csv = fopen(_parameter->_RLOut_file_csv.c_str(), "w");
    if( !Output_csv ) {
        cerr << "[ERROR] Output_csv file [" << _parameter->_RLOut_file_csv << "] not open!!" << endl;
        exit(1);
    }

    fprintf(Output_detail,"%d\t%.2f\n",nf1,grid_width); // Line 1 : Number of FFT, Grid width

    for( int xyz = 0 ; xyz < 3 ; xyz++ ) {
        theta[xyz] = _receptor->angle(xyz);
    }

    fprintf(Output_detail,"\t%f\t%f\t%f\n",theta[0],theta[1],theta[2]); // Line 2 : angle of Receptor

    for( int xyz = 0 ; xyz < 3 ; xyz++ ) {
        coordr[xyz] = _receptor->center(xyz);
        coordl[xyz] = _ligand->center(xyz);
    }

    fprintf(Output_detail,"%s\t%f\t%f\t%f\n",_receptor->input_file().c_str(),coordr[0],coordr[1],coordr[2]); // Line 3 : Receptor PDB file, center coordinate

    fprintf(Output_detail,"%s\t%f\t%f\t%f\n",_ligand->input_file().c_str(),coordl[0],coordl[1],coordl[2]); // Line 4 : Ligand PDB file, center coordinate


    fprintf(Output_detail,"\nElec score ratio: %6.2f\nACE  score ratio: %6.2f\n",_parameter->_Elec_ratio,_parameter->_ACE_ratio);


    // Line 5- : show "correct" position and positions around correct position

    for(int i=0; i<3; i++) { // calculate correct position
        diff_dis[i] = coordr[i] - coordl[i];
        cor_coord[i] =round(diff_dis[i] / grid_width);
        if(cor_coord[i] < 0) cor_coord[i] += nf1;
    }
    fprintf(Output_detail,"\ncorrect position: angle = (%5.2f,%5.2f,%5.2f), trans = (%3d,%3d,%3d) \n\n",theta[0],theta[1],theta[2],cor_coord[0],cor_coord[1],cor_coord[2]);

    fprintf(Output_detail," *** score around correct position *** \n\n");

    string column_list( "Rank   Angle               Trans         total       rPSC     rPSC(+)  rPSC(-)  ELEC       RecACE     LCARMSD NearNat" );
    //sprintf(buff," \t \t \t \t \t \t \t \t \ttotal\trPSC\tELEC\tRecACE\tRMSD\tNearNat");
    fprintf(Output_detail,"%s\n",column_list.c_str());

    // Prepare voxel scores -------------------------------------------------
    _receptor->shift_center(_Mol_coord);
    create_voxel(_receptor); // Receptor --------
    for (int i=0; i<ng3; i++) {
        rec_r[i] = _Fwork[i];
        rec_i[i] = _Fwork[i+ng3];
    }

    _parameter->_ACE_ratio = 0.0; // Receptor(except ACE) --------
    create_voxel(_receptor);

    for (int i=0; i<ng3; i++) rec_r_rpsc[i] = _Fwork[i];

    for(int j = 0 ; j < 3 ; j++ ) theta[j] = 0.0; // Ligand --------
    ligand_rotationz(theta);

    //for RMSD calc
    for (int i=0; i<natom; i++) {
        for (int xyz=0; xyz<3; xyz++) {
            const int j = i*3+xyz;
            lig_atom_coord_orig[j] = _ligand->coordinate(i,xyz);
        }
        if (_ligand->_Atomtype[i] == 1) { //_Atomtype 1:CA, 2:C,N, 0:others
            natom_lca++;
        }
        if (_ligand->_Atomtype[i] > 0) {
            natom_lbb++;
        }
    }

    create_voxel(_ligand);
    for (int i=0; i<ng3; i++) {
        lig_r[i] = _Fwork[i];
        lig_i[i] = _Fwork[i+ng3];
    }

    for(int rank=0; rank<7; rank++) {
        score_rpsc = 0.0;
        score_elec = 0.0;
        score_ace = 0.0;
        score_total = 0.0;

        score_rpsc_plus = 0.0;
        score_rpsc_minus = 0.0;

        if(rank==1) cor_coord[0]--;
        if(rank==2) cor_coord[0]=cor_coord[0]+2;
        if(rank==3) {
            cor_coord[0]--;
            cor_coord[1]--;
        }
        if(rank==4) cor_coord[1]=cor_coord[1]+2;
        if(rank==5) {
            cor_coord[1]--;
            cor_coord[2]--;
        }
        if(rank==6) cor_coord[2]=cor_coord[2]+2;
        //printf("\n %d, %d, %d \n",cor_coord[0],cor_coord[1],cor_coord[2]);

        for( int xyz = 0 ; xyz < 3 ; xyz++ ) tem_coord[xyz] = (3*ng1-cor_coord[xyz]) % (2*ng1);

        if (tem_coord[0]<ng1) {
            la=ng1-tem_coord[0];
            lb=ng1;
        } else {
            la=0;
            lb=2*ng1-tem_coord[0];
        }
        if (tem_coord[1]<ng1) {
            ma=ng1-tem_coord[1];
            mb=ng1;
        } else {
            ma=0;
            mb=2*ng1-tem_coord[1];
        }
        if (tem_coord[2]<ng1) {
            na=ng1-tem_coord[2];
            nb=ng1;
        } else {
            na=0;
            nb=2*ng1-tem_coord[2];
        }
        for (int l=la; l < lb; l++) {
            const int rxc = ng2 * (l+tem_coord[0]-ng1);
            const int lxc = ng2 * l;
            for (int m=ma; m < mb; m++) {
                const int ryc = rxc + ng1 * (m+tem_coord[1]-ng1);
                const int lyc = lxc + ng1 * m;
                for (int n=na; n < nb; n++) {
                    const int rzc=ryc+(tem_coord[2]+n-ng1);
                    const int lzc=lyc+n;

                    const float rpsc = rec_r_rpsc[rzc]*_Fwork[lzc];
                    if (rpsc < 0) {
                        score_rpsc_minus += rpsc;
                        count_rpsc_minus++;
                    } else if(rpsc > 0) {
                        score_rpsc_plus += rpsc;
                        count_rpsc_plus++;
                    }


                    score_rpsc += rec_r_rpsc[rzc]*lig_r[lzc];
                    score_elec += rec_i[rzc]*lig_i[lzc];
                    score_ace += (rec_r[rzc] - rec_r_rpsc[rzc])*lig_r[lzc];
                }
            }
        }
        score_total = score_rpsc + score_elec + score_ace;

        //RMSD calc
        rmsd_lca = 0.0;
        rmsd_lbb = 0.0;
        float r=0.0;
        for (int i=0; i<natom; i++) {
            for (int xyz=0; xyz<3; xyz++) {
                const int j = i*3+xyz;
                if (_ligand->_Atomtype[i] > 0) { //_Atomtype 1:CA, 2:C,N, 0:others
                    lig_atom_coord_temp[j] = _ligand->coordinate(i,xyz) + grid_width*(tem_coord[xyz]-ng1) + diff_dis[xyz];
                    r+=pow(lig_atom_coord_temp[j] - lig_atom_coord_orig[j],2.0);
                    //printf("r %f r/n %f rd %f\n",r,r/j,sqrt(r/j));

                    if (_ligand->_Atomtype[i] == 1) rmsd_lca += pow(lig_atom_coord_temp[j] - lig_atom_coord_orig[j],2.0);
                    rmsd_lbb += pow(lig_atom_coord_temp[j] - lig_atom_coord_orig[j],2.0);
                    //printf("lca %f %f lbb %f %f\n",rmsd_lca,sqrt(rmsd_lca/natom_lca),rmsd_lbb,sqrt(rmsd_lbb/natom_lbb));
                }
            }
        }
        rmsd_lca = sqrt(rmsd_lca/natom_lca);
        rmsd_lbb = sqrt(rmsd_lbb/natom_lbb);
        //printf("lca %f lbb %f all %f\n",rmsd_lca,rmsd_lbb, sqrt(r/natom));
        if(rmsd_lca < nn_thr) {
            near_native = 1;
        } else {
            near_native = 0;
        }
        //RMSD calc
        fprintf(Output_detail,"%5d: (%5.2f,%5.2f,%5.2f) (%3d,%3d,%3d) %11.2f %8.0f %8.0f %8.0f %10.2f %10.2f %7.3f %7d\n"
                ,rank+1,theta[0],theta[1],theta[2],cor_coord[0],cor_coord[1],cor_coord[2]
                ,score_total,score_rpsc,score_rpsc_plus,score_rpsc_minus,score_elec,score_ace,rmsd_lca,near_native);
    }

    // Line 19- : docking results
    fprintf(Output_detail,"\n *** docking results *** \n\n");
    fprintf(Output_detail,"%s\n",column_list.c_str());
    fprintf(Output_csv,"Rank,Anglex,Angley,Anglez,Transi,Transj,Transk,total,rPSC,rPSC_gain,rPSC_penalty,ELEC ,RecACE,LCARMSD,NearNat\n");

    for( int rank = 0 ; rank < num_output ; rank++ ) {
        const int angle_num = _fft_process->top_index(rank,0); // prepare angle & coord
        for( int xyz = 0 ; xyz < 3 ; xyz++ ) {
            theta[xyz] = _parameter->_Zangle[3*angle_num+xyz];
            cor_coord[xyz] = _fft_process->top_index(rank,xyz + 1);
            tem_coord[xyz] = (3*ng1-cor_coord[xyz]) % (2*ng1);
        }

        ligand_rotationz(theta);
        create_voxel(_ligand);

        if (tem_coord[0]<ng1) {
            la=ng1-tem_coord[0];
            lb=ng1;
        } else {
            la=0;
            lb=2*ng1-tem_coord[0];
        }
        if (tem_coord[1]<ng1) {
            ma=ng1-tem_coord[1];
            mb=ng1;
        } else {
            ma=0;
            mb=2*ng1-tem_coord[1];
        }
        if (tem_coord[2]<ng1) {
            na=ng1-tem_coord[2];
            nb=ng1;
        } else {
            na=0;
            nb=2*ng1-tem_coord[2];
        }
        //const int tem_num_grid3 = (lb-la)*(mb-ma)*(nb-na);

        score_rpsc = 0.0;
        score_elec = 0.0;
        score_ace = 0.0;
        score_total = 0.0;

        score_rpsc_plus = 0.0;
        score_rpsc_minus = 0.0;
        score_elec_plus = 0.0;
        score_elec_minus = 0.0;
        score_ace_plus = 0.0;
        score_ace_minus = 0.0;

        count_rpsc_plus = 0;
        count_rpsc_minus = 0;
        count_elec_plus = 0;
        count_elec_minus = 0;
        count_ace_plus = 0;
        count_ace_minus = 0;

        for (int l=la; l < lb; l++) {
            const int rxc = ng2 * (l+tem_coord[0]-ng1);
            const int lxc = ng2 * l;
            for (int m=ma; m < mb; m++) {
                const int ryc = rxc + ng1 * (m+tem_coord[1]-ng1);
                const int lyc = lxc + ng1 * m;
                for (int n=na, rzc=ryc+(tem_coord[2]+na-ng1), lzc=lyc+na, lzc_i=lyc+na+ng3; n < nb; n++, rzc++, lzc++, lzc_i++) { // rzc & lzc increment for speed
                    const float rpsc = rec_r_rpsc[rzc]*_Fwork[lzc];
                    const float elec = rec_i[rzc]*_Fwork[lzc_i];
                    const float ace = (rec_r[rzc] - rec_r_rpsc[rzc])*_Fwork[lzc];
                    if (rpsc < 0) {
                        score_rpsc_minus += rpsc;
                        count_rpsc_minus++;
                    } else if(rpsc > 0) {
                        score_rpsc_plus += rpsc;
                        count_rpsc_plus++;
                    }
                    if (elec < 0) {
                        score_elec_minus += elec;
                        count_elec_minus++;
                    } else if(elec > 0) {
                        score_elec_plus += elec;
                        count_elec_plus++;
                    }
                    if (ace < 0) {
                        score_ace_minus += ace;
                        count_ace_minus++;
                    } else if(ace > 0) {
                        score_ace_plus += ace;
                        count_ace_plus++;
                    }

                }
            }
        }

        score_rpsc =score_rpsc_plus + score_rpsc_minus;
        score_elec =score_elec_plus + score_elec_minus;
        score_ace = score_ace_plus + score_ace_minus;

        score_total = score_rpsc + score_elec + score_ace;

        for( int xyz = 0 ; xyz < 3 ; xyz++ ) { // PI wo koeru theta ha nai kiga suru
            if( theta[xyz] >= PI ) {
                theta[xyz] -= PI*2.0;
            }
        }

        //*
        //RMSD calc
        rmsd_lca = 0.0;
        rmsd_lbb = 0.0;
        for (int i=0; i<natom; i++) {
            for (int xyz=0; xyz<3; xyz++) {
                const int j = i*3+xyz;
                if (_ligand->_Atomtype[i] > 0) { //_Atomtype 1:CA, 2:C,N, 0:others
                    lig_atom_coord_temp[j] = _ligand->coordinate(i,xyz) + grid_width*(tem_coord[xyz]-ng1) + diff_dis[xyz];

                    if (_ligand->_Atomtype[i] == 1) rmsd_lca += pow(lig_atom_coord_temp[j] - lig_atom_coord_orig[j],2.0);
                    rmsd_lbb += pow(lig_atom_coord_temp[j] - lig_atom_coord_orig[j],2.0);
                }
            }
        }
        rmsd_lca = sqrt(rmsd_lca/natom_lca);
        rmsd_lbb = sqrt(rmsd_lbb/natom_lbb);
        if(rmsd_lca < nn_thr) {
            near_native = 1;
        } else {
            near_native = 0;
        }
        //*/

        fprintf(Output_detail,"%5d: (%5.2f,%5.2f,%5.2f) (%3d,%3d,%3d) %11.2f %8.0f %8.0f %8.0f %10.2f %10.2f %7.3f %7d\n"
                ,rank+1,theta[0],theta[1],theta[2],cor_coord[0],cor_coord[1],cor_coord[2]
                ,score_total,score_rpsc,score_rpsc_plus,score_rpsc_minus,score_elec,score_ace,rmsd_lca,near_native);
        fprintf(Output_csv,"%d, %f,%f,%f,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%d\n"
                ,rank+1,theta[0],theta[1],theta[2],cor_coord[0],cor_coord[1],cor_coord[2]
                ,score_total,score_rpsc,score_rpsc_plus,score_rpsc_minus,score_elec,score_ace,rmsd_lca,near_native);

    }

    fclose(Output_detail);
    fclose(Output_csv);

    delete [] rec_r;
    delete [] rec_r_rpsc;
    delete [] rec_i;
    delete [] lig_r;
    delete [] lig_i;
    delete [] lig_atom_coord_temp;
    delete [] lig_atom_coord_orig;

    return;
}


//============================================================================//
void DockingTable::output_calc_time_log()
//============================================================================//
{
    char buff[4096];
    const int ng1 = _Num_grid;
    const int nf1 = ng1 * 2;
    const int num_per_angle = _parameter->_Num_sort;
    const int natom = _ligand->num_atoms();


    time_t timer;
    struct tm *local;
    timer = time(NULL);
    local = localtime(&timer);
    int newfile_flag = 0;

    string  rfile = _receptor->input_file();
    string  lfile = _ligand->input_file();
    int ipr, ipl;

    while(1) {
        ipr   = rfile.rfind("/");

        if( ipr == (int) string::npos ) {
            break;
        } else {
            rfile = rfile.substr(ipr+1);
        }
    }

    ipr   = rfile.rfind(".");
    rfile = rfile.substr(0,ipr);

    while(1) {
        ipl   = lfile.rfind("/");

        if( ipl == (int) string::npos ) {
            break;
        } else {
            lfile = lfile.substr(ipl+1);
        }
    }

    ipl   = lfile.rfind(".");
    lfile = lfile.substr(0,ipr);

    const float t_total = _exec_logger->_cputime->t1_initialize
        +_exec_logger->_cputime->t2_receptor_process
        +_exec_logger->_cputime->t3_docking_total
        +_exec_logger->_cputime->t4_docking_output_detail
        +_exec_logger->_cputime->t5_docking_output;

    string calc_id = _parameter->calc_id;

    string file_path;

    file_path="./megadock_calc_time_log.csv";

    ifstream ifs(file_path.c_str(), ios::in);
    if (!ifs) {
        newfile_flag = 1;
    }
    ifs.close();

    ofstream Output_detail(file_path.c_str(), std::ios::app); // open Output file (append mode)
    if( !Output_detail ) {
        printf("cannot open calc-log file\n");
    } else {
        sprintf(buff,"%s,%d,%d,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f"
                , calc_id.c_str()
                , _parallel->nproc2(), _parallel->num_gpu()
                , rfile.c_str(), lfile.c_str()
                , ng1, nf1, natom
                , num_per_angle, _parameter->_Rotation_angle_set, _parameter->_Old_voxel_flag
                , _parameter->fft_base_set, _parameter->lig_elec_serial_flag
                , t_total
                ,_exec_logger->_cputime->t1_initialize
                ,_exec_logger->_cputime->t2_receptor_process
                ,_exec_logger->_cputime->t3_docking_total
                ,_exec_logger->_cputime->t4_docking_output_detail
                ,_exec_logger->_cputime->t5_docking_output
                ,_exec_logger->_cputime->t3_1_ligand_voxelization
                ,_exec_logger->_cputime->t3_2_fftprocess_ligand_fft
                ,_exec_logger->_cputime->t3_3_fftprocess_convolution
                ,_exec_logger->_cputime->t3_4_fftprocess_fft_inverse
                ,_exec_logger->_cputime->t3_5_fftprocess_score_sort
                ,_exec_logger->_cputime->t3_6_fftprocess_sort_index
                ,_exec_logger->_cputime->t3_1_1_ligvoxgpu_copy_htod
                ,_exec_logger->_cputime->t3_1_2_ligvoxgpu_kernel_init
                ,_exec_logger->_cputime->t3_1_3_ligvoxgpu_kernel_fill_core
                ,_exec_logger->_cputime->t3_1_4_ligvoxgpu_kernel_cut_surf
                ,_exec_logger->_cputime->t3_1_5_ligvoxgpu_kernel_fill_surf
                ,_exec_logger->_cputime->t3_1_6_ligvoxgpu_kernel_elec
                ,_exec_logger->_cputime->t3_1_7_ligvoxgpu_kernel_set_array
                ,_exec_logger->_cputime->t7_offload_transfer
                );
        Output_detail << buff << endl;
    }

    Output_detail.close(); // close Output file

    return;
}

