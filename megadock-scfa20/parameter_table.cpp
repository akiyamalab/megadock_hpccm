/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : ParameterTable
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#include "parameter_table.h"

#include <string>
using namespace std;

//============================================================================//
void ParameterTable::process_args(int argc, char *argv[])
//============================================================================//
{
#ifdef MPI_DP
    for( int i = 0 ; i < 3 ; i++ ) { //Receptor PDB, Ligand PDB, Outfile
        _IO_flag[i] = 0;
    }
#endif

    optind = 1;
    int ch;

    std::string pdb_ext_r_b("_r_b.pdb"); // for input support
    std::string pdb_ext_l_b("_l_b.pdb");
    std::string pdb_ext_r_u("_r_u.pdb");
    std::string pdb_ext_l_u("_l_u.pdb");
    std::string detail_ext(".detail");
    std::string csv_ext(".csv");

    while( (ch = getopt(argc, argv,
                        /*
                            Line 1:input files
                            Line 2:output
                            Line 3:parameters
                            Line 4:calculation options
                        */
#ifdef MPI_DP
                        "\
    R:L:B:U:\
    o:ON:t:E:\
    a:b:e:d:F:v:\
    Dhy:z:r:f:i:jkl:S:G:T:"
#else
                        "\
    I:\
    ON:t:E:\
    a:b:e:d:F:v:\
    Dhy:z:r:f:i:jkl:S:G:T:"
#endif
                        )) != -1 ) {
        switch( ch ) {
            //-----------------------------------------------------------------
            // input PDB file -------------------------------------------------

#ifdef MPI_DP
        case 'R':
            _RecPDB_file  = optarg;
            _IO_flag[0]   = 1;
            break;
        case 'L':
            _LigPDB_file  = optarg;
            _IO_flag[1]   = 1;
            break;
        case 'B': // for input support
            _RecPDB_file  = optarg+pdb_ext_r_b;
            _IO_flag[0]   = 1;
            _LigPDB_file  = optarg+pdb_ext_l_b;
            _IO_flag[1]   = 1;
            break;
        case 'U': // for input support
            _RecPDB_file  = optarg+pdb_ext_r_u;
            _IO_flag[0]   = 1;
            _LigPDB_file  = optarg+pdb_ext_l_u;
            _IO_flag[1]   = 1;
            break;

            //-----------------------------------------------------------------
            // output options -------------------------------------------------

        case 'o':
            _RLOut_file = optarg;
            _RLOut_file_detail = optarg+detail_ext;
            _RLOut_file_csv = optarg+csv_ext;
            if(_RLOut_file.length() > 4) {
                if(_RLOut_file.substr(_RLOut_file.length()-4)==".out") {
                    _RLOut_file_detail = _RLOut_file.substr(0,_RLOut_file.length()-4)+detail_ext;;
                    _RLOut_file_csv = _RLOut_file.substr(0,_RLOut_file.length()-4)+csv_ext;;
                }
            }
            _IO_flag[2]   = 1;
            break;
#else
        case 'I':
            _Table_file   = optarg;
            break;
#endif

            //-----------------------------------------------------------------
            // output options -------------------------------------------------

        case 'O':
            detail_output_flag = 1;
            break;
        case 'N':
            _Num_output   = atoi(optarg);
            _Num_output_flag  = 1;
            cout << "# Number of output = " << _Num_output << endl;
            break;
        case 't':
            _Num_sort = atoi(optarg);
            cout << "# Set number of scores per one angle = " << _Num_sort << endl;
            break;
        case 'E':
            calc_time_log_output_flag = atoi(optarg);
            break;
        case 'i':
            calc_id = optarg;
            break;

            //-----------------------------------------------------------------
            // setting parameters ---------------------------------------------

        case 'e':
            _Elec_ratio   = atof(optarg);
            if (f1_flag == 1) {
                printf("Do not use -f 1 and -e option at the same time.\n");
                exit(1);
            }
            f1_flag = -1; // do not use rPSC only
            cout << "# Set electric term ratio = "<< _Elec_ratio << endl;
            break;
        case 'd':
            _ACE_ratio = atof(optarg);
            if (f1_flag == 1 || f2_flag == 1) {
                printf("Do not use -f 1 or 2 and -d option at the same time.\n");
                exit(1);
            }
            f1_flag = -1; // do not use rPSC only
            f2_flag = -1; // do not use rPSC and Elec function
            cout << "# Set ACE term ratio = " << _ACE_ratio << endl;
            break;
        case 'a':
            _rPSC_param_rec_core = atof(optarg); // Receptor core
            break;
        case 'b':
            _rPSC_param_lig_core = atof(optarg); // Ligand core
            break;
        case 'F':
            _Num_fft  = atoi(optarg);
            _Num_grid = _Num_fft / 2;
            _Num_fft_flag = 1;
            cout << "# Number of FFT N  = " << _Num_fft << endl;
            break;
        case 'v':
            grid_width    = atof(optarg);
            cout << "# Set voxel size   = " << grid_width << endl;
            break;

            //-----------------------------------------------------------------
            // setting calculation mode ---------------------------------------
        case 'T':
            _Num_thread_limit = atoi(optarg);
            if(_Num_thread_limit<1) {
                printf("Please set a positive number (Number of cores).\n");
                exit(1);
            }
            omp_set_num_threads(_Num_thread_limit);
            break;
        case 'G':
            _Num_GPU_limit = atoi(optarg);
            if(_Num_GPU_limit<0) {
                printf("Please set a non-negative number (Number of GPUs you want to use).\n");
                exit(1);
            }
            break;
        case 'f':
            _Score_func   = atoi(optarg);
            cout << "#Set score Function = "   << _Score_func << endl;
            assert( 1 <= _Score_func && _Score_func <= 3 );
            if(_Score_func == 1) {
                _Elec_ratio   = 0.0;
                _ACE_ratio    = 0.0;
                if (f1_flag == -1) {
                    printf("Do not use -f 1 and -e/-d option at the same time.\n");
                    exit(1);
                }
                f1_flag       = 1;
            } else if(_Score_func == 2) {
                _ACE_ratio    = 0.0;
                if (f2_flag == -1) {
                    printf("Do not use -f 2 and -d option at the same time.\n");
                    exit(1);
                }
                f2_flag       = 1;
            }
            _Score_func = 3;
            break;
        case 'r':
            _Rotation_angle_set = atoi(optarg);
            if(_Rotation_angle_set == 54000) {
                cout << "# Set 54,000 rotational angles (6 degree)"<< endl;
            } else if(_Rotation_angle_set == 1) {
                cout << "# Set 1 rotational angle (test mode)"<< endl;
            } else if(_Rotation_angle_set == 3) {
                cout << "# Set 3 rotational angle (test mode)"<< endl;
            } else if(_Rotation_angle_set == 24) {
                cout << "# Set 24 rotational angles (test mode)"<< endl;
            } else if(_Rotation_angle_set == 360) {
                cout << "# Set 360 rotational angles (test mode)"<< endl;
            }
            break;
        case 'D':
            _Rotation_angle_set = 54000;
            cout << "# Set 54,000 rotational angles (6 degree)"<< endl;
            break;

            //-----------------------------------------------------------------
            // help -----------------------------------------------------------
        case 'j':
            tem_flag1 = 1;
            break;
        case 'k':
            tem_flag2 = 1;
            break;
        case 'S':
            lig_elec_serial_flag = atoi(optarg);
            break;
        case 'l':
            fft_base_set = atoi(optarg);
            break;

        case 'h':
            break;
        }
    }
#ifdef MPI_DP
    pdb_step();
#endif
}

#ifdef MPI_DP
//============================================================================//
void ParameterTable::pdb_step()
//============================================================================//
{
    if( !_IO_flag[0] ) {
        cerr << "[ERROR] Receptor PDB file is not specified!!" << endl;
        usage();
        exit(1);
    }

    if( !_IO_flag[1] ) {
        cerr << "[ERROR] Ligand PDB file is not specified!!" << endl;
        usage();
        exit(1);
    }

    if( !_IO_flag[2] ) {
        string  rfile = _RecPDB_file;
        string  lfile = _LigPDB_file;
        string  ofile;
        int     ipr;
        int     ipl;

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
        rfile = rfile + "-";

        while(1) {
            ipl   = lfile.rfind("/");

            if( ipl == (int) string::npos ) {
                break;
            } else {
                lfile = lfile.substr(ipl+1);
            }
        }

        ipl   = lfile.rfind(".");
        lfile = lfile.substr(0,ipl);

        ofile = rfile + lfile + ".out";
        _RLOut_file = ofile;
        ofile = rfile + lfile + ".detail";
        _RLOut_file_detail = ofile;
        ofile = rfile + lfile + ".csv";
        _RLOut_file_csv = ofile;
    }

    //cout << "#Receptor = " << _RecPDB_file << endl;
    //cout << "#Ligand   = " << _LigPDB_file << endl;
    cout << "# Output file = " << _RLOut_file << endl;

    return;
}
#endif

//============================================================================//
void ParameterTable::initialize(ParameterTable *pparameter)
//============================================================================//
{
#ifndef MPI_DP
    _Table_file = pparameter->_Table_file;
#endif
    _RLOut_file = pparameter->_RLOut_file;
    _RLOut_file_detail = pparameter->_RLOut_file_detail; 
    _RLOut_file_csv = pparameter->_RLOut_file_csv;
    calc_id = pparameter->calc_id;
    detail_output_flag = pparameter->detail_output_flag;
    calc_time_log_output_flag = pparameter->calc_time_log_output_flag;

    _Num_grid = pparameter->_Num_grid;
    _Num_fft = pparameter->_Num_fft;
    _Num_fft_flag = pparameter->_Num_fft_flag;
    _Num_atom_max = pparameter->_Num_atom_max;
    _Num_output = pparameter->_Num_output;
    _Num_output_flag = pparameter->_Num_output_flag;
    _Num_thread_limit = pparameter->_Num_thread_limit;
    _Num_GPU_limit = pparameter->_Num_GPU_limit;

    _Score_func = pparameter->_Score_func;
    _Num_sort = pparameter->_Num_sort;
    _Elec_ratio = pparameter->_Elec_ratio;
    _ACE_ratio = pparameter->_ACE_ratio;
    grid_width = pparameter->grid_width;
    ligand_max_edge = pparameter->ligand_max_edge;
    _Rotation_angle_set = pparameter->_Rotation_angle_set;
    fft_base_set = pparameter->fft_base_set;
    lig_elec_serial_flag = pparameter->lig_elec_serial_flag;
    fft_library_type = pparameter->fft_library_type;

    tem_flag1 = pparameter->tem_flag1;
    tem_flag2 = pparameter->tem_flag2;
    tem_flag3 = pparameter->tem_flag3;
    tem_flag4 = pparameter->tem_flag4;
    f1_flag = pparameter->f1_flag;
    f2_flag = pparameter->f2_flag;
    
    _Old_voxel_flag = pparameter->_Old_voxel_flag;
    _Grid_space_rec = pparameter->_Grid_space_rec;
    _Grid_space_lig = pparameter->_Grid_space_lig;
    
    _rPSC_param_rec_core = pparameter->_rPSC_param_rec_core;
    _rPSC_param_lig_core = pparameter->_rPSC_param_lig_core;

    _Num_rot_angles = pparameter->_Num_rot_angles;
    _Charmmr = pparameter->_Charmmr;
    _Charmmc = pparameter->_Charmmc;
    _ACE = pparameter->_ACE;


    _Zangle = new float[_Num_rot_angles*3];
    for( int i = 0 ; i < _Num_rot_angles*3 ; i++ ) {
        _Zangle[i] = pparameter->_Zangle[i];
    }
}

//============================================================================//
void ParameterTable::output_file_name(const string rec_file, const string lig_file)
//============================================================================//
{
    string  rfile = rec_file;
    string  lfile = lig_file;
    string  ofile;
    int     ipr;
    int     ipl;

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
    rfile = rfile + "-";

    while(1) {
        ipl   = lfile.rfind("/");

        if( ipl == (int) string::npos ) {
            break;
        } else {
            lfile = lfile.substr(ipl+1);
        }
    }

    ipl   = lfile.rfind(".");
    lfile = lfile.substr(0,ipl);

    ofile = rfile + lfile + ".out";
    _RLOut_file = ofile;
    ofile = rfile + lfile + ".detail";
    _RLOut_file_detail = ofile;
    ofile = rfile + lfile + ".csv";
    _RLOut_file_csv = ofile;

    return;
}
