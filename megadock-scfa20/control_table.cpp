/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : ControlTable
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#include "control_table.h"

//============================================================================//
void ControlTable::initialize(bool verbose)
//============================================================================//
{
    struct timeval et1, et2;
    gettimeofday(&et1,NULL);

    // Number of processors limitation
    const int thread_limit = _parameter->_Num_thread_limit;
    const int gpu_limit = _parameter->_Num_GPU_limit;

    if(_parallel->nproc2() > thread_limit) {
        _parallel->nproc2(thread_limit);
    }

    if(_parallel->num_gpu() > gpu_limit || _parallel->num_gpu() > _parallel->nproc2()) {
        _parallel->num_gpu( min(gpu_limit, (int)_parallel->nproc2()) );
    }
    if (verbose)
        printf("# Using %3d CPU cores, %d GPUs\n\n", _parallel->nproc2(), _parallel->num_gpu());

    gettimeofday(&et2,NULL);
    _exec_logger->_cputime->t1_initialize += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
}

//============================================================================//
#ifdef MPI_DP
void ControlTable::prepare()
#else
void ControlTable::prepare(string rec_file, string lig_file, string out_file)
#endif
//============================================================================//
{
    int ngrid;
    vector<int> ngrid_table;

    struct timeval et1, et2;
    gettimeofday(&et1,NULL);

#ifdef MPI_DP
    string rec_file = _parameter->_RecPDB_file;
    string lig_file = _parameter->_LigPDB_file;
#else
    if (out_file == "") {
        _parameter->output_file_name(rec_file, lig_file);
    } else {
        string detail_ext = ".detail", csv_ext = ".csv";
        _parameter->_RLOut_file = out_file;
        _parameter->_RLOut_file_detail = out_file+detail_ext;
        _parameter->_RLOut_file_csv = out_file+csv_ext;
        if(_parameter->_RLOut_file.length() > 4) {
            if(_parameter->_RLOut_file.substr(_parameter->_RLOut_file.length()-4)==".out") {
                _parameter->_RLOut_file_detail = _parameter->_RLOut_file.substr(0,_parameter->_RLOut_file.length()-4)+detail_ext;;
                _parameter->_RLOut_file_csv = _parameter->_RLOut_file.substr(0,_parameter->_RLOut_file.length()-4)+csv_ext;;
            }
        }
    }
#endif
    // Receptor<ParameterTable>
    _receptor = new Receptor<ParameterTable>(rec_file);
    _exec_logger->rec_filename = rec_file;
    _receptor->initialize(_parameter);
    _exec_logger->record_malloc( sizeof(float)*_receptor->num_atoms()*3 ); //Atom coordinate

    // Ligand<ParameterTable>
    _ligand = new Ligand<ParameterTable>(lig_file);
    _exec_logger->lig_filename = lig_file;
    _ligand->initialize(_parameter);
    _exec_logger->record_malloc( sizeof(float)*_ligand->num_atoms()*3 ); //Atom coordinate

    _exec_logger->_Num_fft_flag = _parameter->_Num_fft_flag;
    if( !_parameter->_Num_fft_flag ) {
        switch (_parameter->fft_base_set) {
        case 13:
            gridtable_13base_normal(ngrid,ngrid_table);
            break;
        case 7:
            gridtable_07base_normal(ngrid,ngrid_table);
            break;
        case 11:
            gridtable_11base_normal(ngrid,ngrid_table);
            break;
        case 0:
            gridtable_fftw_custom(ngrid,ngrid_table);
            break;
        case 1:
            gridtable_cufft_custom(ngrid,ngrid_table);
            break;
        }
        autogridr(ngrid,ngrid_table);
        autogridl(ngrid,ngrid_table);
    } else {
        checkgridr();
        checkgridl();
    }

    // DockingTable
    _docking = new DockingTable(_exec_logger,_parallel,_parameter,_receptor,_ligand);
    _docking->initialize();

    gettimeofday(&et2,NULL);
    _exec_logger->_cputime->t1_initialize += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
}

//============================================================================//
void ControlTable::autogridr(const int &ngrid,vector<int> &ngrid_table)
//============================================================================//
{
    int       num_grid = 1;
    float     size, size_rec = 0.0;

    for( int i = 0 ; i < 3 ; i++ ) {
        size = _receptor->edge(i,1) - _receptor->edge(i,0);
        
        //printf(" %f, %f\n",_receptor->edge(i,1),_receptor->edge(i,0));

        if( size > size_rec ) {
            size_rec = size;
        }
    }

    _exec_logger->rec_max_size = size_rec;

    size_rec += 2.0 * _parameter->_Grid_space_rec;
    _exec_logger->rec_voxel_size = size_rec;

    num_grid = 1 + int(size_rec / _parameter->grid_width);

    for( int i = 0 ; i < ngrid ; i++ ) {
        if( ngrid_table[i] >= num_grid ) {
            num_grid = ngrid_table[i];
            break;
        }
    }

    _receptor->num_grid(num_grid);
    _exec_logger->rec_num_grid = num_grid;

    return;
}

//============================================================================//
void ControlTable::autogridl(const int &ngrid,vector<int> &ngrid_table)
//============================================================================//
{
    int       num_grid = 1;
    float     size_lig = 0.0;
    float     x1, y1, z1, x2, y2, z2, d2;
    const int na  = _ligand->num_atoms();

    for( int i = 0 ; i < na-1 ; i++ ) {
        x1 = _ligand->coordinate(i,0);
        y1 = _ligand->coordinate(i,1);
        z1 = _ligand->coordinate(i,2);

        for( int j = i+1 ; j < na ; j++ ) {
            x2 = _ligand->coordinate(j,0);
            y2 = _ligand->coordinate(j,1);
            z2 = _ligand->coordinate(j,2);

            d2 = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1);

            if( d2 > size_lig ) {
                size_lig = d2;
            }
        }
    }

    size_lig = sqrt(size_lig);
    _exec_logger->lig_max_size = size_lig;

    size_lig += 2.0 * _parameter->_Grid_space_lig;
    _exec_logger->lig_voxel_size = size_lig;
    
    _parameter->ligand_max_edge = size_lig;

    num_grid = 1 + int(size_lig / _parameter->grid_width);

    for( int i = 0 ; i < ngrid ; i++ ) {
        if( ngrid_table[i] >= num_grid ) {
            num_grid = ngrid_table[i];
            break;
        }
    }

    _ligand->num_grid(num_grid);
    _exec_logger->lig_num_grid = num_grid;

    return;
}

//============================================================================//
void ControlTable::checkgridr()
//============================================================================//
{
    float     size, size_rec  = 0.0;
    const int num_grid    = _parameter->_Num_grid;
    const float   search_length   = _parameter->grid_width * num_grid;

    for( int i = 0 ; i < 3 ; i++ ) {
        size = _receptor->edge(i,1) - _receptor->edge(i,0);

        if( size > size_rec ) {
            size_rec = size;
        }
    }

    _exec_logger->rec_max_size = size_rec;

    size_rec += 2.0*_parameter->_Grid_space_rec;
    _exec_logger->rec_voxel_size = size_rec;

    if( size_rec > search_length ) {
        cerr << "[ERROR] Receptor data is too big!!\n";
        exit(1);
    }

    _receptor->num_grid(num_grid);
    _exec_logger->rec_num_grid = num_grid;


    return;
}

//============================================================================//
void ControlTable::checkgridl()
//============================================================================//
{
    float     size_lig = 0.0;
    float     x1, y1, z1, x2, y2, z2, d2;
    const int na  = _ligand->num_atoms();
    const int num_grid    = _parameter->_Num_grid;
    const float   search_length   = _parameter->grid_width * num_grid;

    for( int i = 0 ; i < na-1 ; i++ ) {
        x1 = _ligand->coordinate(i,0);
        y1 = _ligand->coordinate(i,1);
        z1 = _ligand->coordinate(i,2);

        for( int j = i+1 ; j < na ; j++ ) {
            x2 = _ligand->coordinate(j,0);
            y2 = _ligand->coordinate(j,1);
            z2 = _ligand->coordinate(j,2);

            d2 = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1);

            if( d2 > size_lig ) {
                size_lig = d2;
            }
        }
    }

    size_lig = sqrt(size_lig);
    _exec_logger->lig_max_size = size_lig;

    size_lig += 2.0*_parameter->_Grid_space_lig;
    _exec_logger->lig_voxel_size = size_lig;

    if( size_lig > search_length ) {
        cerr << "[ERROR] Ligand data is too big!!\n";
        exit(1);
    }

    _ligand->num_grid(num_grid);
    _exec_logger->lig_num_grid = num_grid;

    _exec_logger->grid_width = _parameter->grid_width;

    return;
}

//============================================================================//
void ControlTable::execute()
//============================================================================//
{
    struct timeval et1, et2;

    gettimeofday(&et1,NULL); // Receptor process (voxelization, forward FFT of Receptor)
    _docking->rec_init();
    gettimeofday(&et2,NULL);
    _exec_logger->_cputime->t2_receptor_process += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

    gettimeofday(&et1,NULL); // docking
    _docking->dockz();
    gettimeofday(&et2,NULL);
    _exec_logger->_cputime->t3_docking_total += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

    if(_parameter->detail_output_flag == 1) { // detailed result output
        gettimeofday(&et1,NULL);
        _docking->output_detail();
        gettimeofday(&et2,NULL);
        _exec_logger->_cputime->t4_docking_output_detail += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
    }

    if(_parameter->calc_time_log_output_flag >= 1) { // calculation info
        _docking->output_calc_time_log();
    }

    gettimeofday(&et1,NULL); // normal result output
    _docking->output();
    gettimeofday(&et2,NULL);
    _exec_logger->_cputime->t5_docking_output += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

    _docking->dock_memory_free();

    const int ng1 = _parameter->_Num_grid;
    const int ng3 = ng1*ng1*ng1;
    const int nf1 = ng1*2;
    const int nf3 = nf1*nf1*nf1;
    const int nproc2 = _parallel->nproc2();
    const int natom = _parameter->_Num_atom_max;
    const int nag = natom * ng1;
    const size_t _Memfw = ng3*3+natom*3+nag*3;
    const size_t _Memiw = ng3*2+natom*4;

    //delete docking include delete fft_process, _FFT_rec_r/i[nf3], _FFTWin/out[nf3*nproc2]
    _exec_logger->record_free( sizeof(float)*nf3*2 + sizeof(fftwf_complex)*nf3*2*nproc2);
#ifdef CUFFT
    _exec_logger->record_free( sizeof(cufftComplex)*nf3*2 ); //_in/outBuf
#endif
    _exec_logger->record_free( sizeof(float)*_Memfw*nproc2 + sizeof(int)*_Memiw*nproc2 ); //_F/Iwork
    delete _docking;
    _exec_logger->record_free( sizeof(float)*_ligand->num_atoms()*3 );
    delete _ligand;
    _exec_logger->record_free( sizeof(float)*_receptor->num_atoms()*3 );
    delete _receptor;
    _exec_logger->record_free( sizeof(float)*_parameter->_Num_rot_angles*3 + sizeof(unordered_map<string,float>)*(_parameter->_Charmmr.size() + _parameter->_Charmmc.size() + _parameter->_ACE.size()) );
    //delete _parameter;


    return;
}
