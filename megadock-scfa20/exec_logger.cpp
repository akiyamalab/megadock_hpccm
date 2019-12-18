/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : ExecLogger
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#include "exec_logger.h"

//============================================================================//
void ExecLogger::initialize()
//============================================================================//
{
    _cputime->initialize();
	
    mb = 0.0;

#ifdef CUFFT
    devmem_free  = 0;
    devmem_total = 0;
    devmem_use   = 0;
#endif

    _RLOut_file = "";
    _Num_fft_flag = 0;
    rec_filename = "";
    lig_filename = "";
    rec_max_size = 0.0;
    lig_max_size = 0.0;
    rec_voxel_size = 0.0;
    lig_voxel_size = 0.0;
    rec_num_grid = 0;
    lig_num_grid = 0;
    grid_width = 0.0;
    return;
}

//============================================================================//
void ExecLogger::output(const int myid2)
//============================================================================//
{
    printf("# Output file = %s\n", _RLOut_file.c_str());
    if (!_Num_fft_flag) {
        printf("\nReceptor = %s\n"
               "Receptor max size = %f\n"
               "Required voxel size = %f\n"
               "Number of grid = %d\n"
               "FFT N = %d\n"
               "\nLigand = %s\n"
               "Ligand max size = %f\n"
               "Required voxel size = %f\n"
               "Number of grid = %d\n"
               "FFT N = %d\n",
               rec_filename.c_str(),
               rec_max_size,
               rec_voxel_size,
               rec_num_grid,
               rec_num_grid * 2,
               lig_filename.c_str(),
               lig_max_size,
               lig_voxel_size,
               lig_num_grid,
               lig_num_grid * 2
              );

    } else {
        cout << "\nReceptor max size = " << rec_max_size << endl;
        cout << "Required voxel size = " << rec_voxel_size << endl;
        cout << "\n(Receptor)\n";
        cout << "Number of grid = " << rec_num_grid << endl;
        cout << "FFT N = " << rec_num_grid*2 << endl;
        cout << "Grid size = " << grid_width << endl;
        cout << "\nLigand max size = " << lig_max_size << endl;
        cout << "Required voxel size = " << lig_voxel_size << endl;
        cout << "\n(Ligand)\n";
        cout << "Number of grid = " << lig_num_grid << endl;
        cout << "FFT N = " << lig_num_grid*2 << endl;
        cout << "Grid size = " << grid_width << endl;
    }
    if ( mb < 1000 )
        printf("Memory requirement (/node)  = %.1Lf MB\n",mb); // approximate value
    else 
        printf("Memory requirement (/node)  = %.1Lf GB\n",mb/1024); // approximate value

#ifdef CUFFT
    printf("# GPU Memory : Use %3.1f MB (%4.1f%%), Free %3.1f MB (%4.1f%%), Total %3.1f MB\n",(float)devmem_use/1024.0/1024.0,(float)(100*devmem_use/devmem_total), (float)devmem_free/1024.0/1024.0, (float)(100*devmem_free/devmem_total), (float)devmem_total/1024.0/1024.0);
#endif

    printf(
            "\n---------- Start docking calculations\n"
            "\nLigand = %s\n"
            "Target receptors:\n"
            " %s\n\n", lig_filename.c_str(), rec_filename.c_str()
          );
    //if( !((ang+1)%nc) || _parameter->tem_flag1==1 ) {
    //    printf("   >Ligand rotation = %5d / %5d (%2d)\n",ang+1,_parameter->_Num_rot_angles,myid2);
    //}
    printf("   Thread ID = %2d\n\n",myid2);

    _cputime->output();
}
