/*
 * Copyright (C) 2019 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : Mpidp
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#include "mpidp.h"

#define LASTUPDATED "2014/4/30"
#define MASTER_PROCESS_ID 0
#define MASTER_THREAD_ID 0

#ifndef SYSTEMCALL
int application(int argc,char *argv[]);
#endif

//============================================================================//
int main(int argc,char *argv[])
//============================================================================//
{
    Mpidp		mpidp;
    ofstream	logout;				// log file stream
    double      stime, etime;
    int 		nproc, myid, resultlen;		// for MPI parameters
    char		hostname[MPI_MAX_PROCESSOR_NAME];
    
    struct timeval et1, et2;
    // Preparation using MPI
    MPI_Init(&argc,&argv);
    stime = MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(hostname,&resultlen);
    char *hostall = new char [nproc*MPI_MAX_PROCESSOR_NAME];
    MPI_Gather(hostname,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,
               hostall,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,
               0,MPI_COMM_WORLD);

    // for master process
    if( myid == MASTER_PROCESS_ID ){
        gettimeofday(&et1,NULL);
        string log_file = "./master.log";		// Default log file name
        for( int i = 1 ; i < argc ; i++ ) {
            if( !strncmp(argv[i],"-lg",3) ) {
                log_file = argv[++i];
            }
        }
        
        logout.open(log_file.c_str());
        if( !logout ) {
            cerr << "[ERROR] Log file [" << log_file << "] was not opened!!" << endl;
            MPI_Abort(MPI_COMM_WORLD,1);
            exit(1);
        }
        
        logout << " MEGADOCK ver. 4.0 Master Process"<< endl;
        logout << "      megadock@bi.cs.titech.ac.jp   last updated: "
               << LASTUPDATED << endl << endl;
        logout << "#RANK = " << nproc << endl;
        
        int nprocess = 1;				// # of processes in one core
        string *shost = new string[nproc];		// Hostname
        shost[0] = &hostall[0];
        for( int i = 1 ; i < nproc ; i++ ) {
            shost[i] = &hostall[i*MPI_MAX_PROCESSOR_NAME];
            if( shost[i] == shost[0] ) {
                nprocess ++;
            }
            else {
                break;
            }
        }
        logout << "#Node = " << nproc/nprocess
               << " (#RANK/Node = " << nprocess << ")" << endl;
        
        logout << "\n used nodes list(id) :";
        for( int i = 0 ; i < nproc ; i++ ) {
            if( i % 5 == 0 ) {
                logout << endl;
            }
            logout << "  " << &hostall[i*MPI_MAX_PROCESSOR_NAME] << "(" << i << ")";
        }
        logout << endl << endl;
        logout.flush();
        delete [] shost;
    }
    
    if( myid == MASTER_PROCESS_ID ) {			// for master
        //int ntry;				// Upper limit of the number of retries
        int eflag;				// return flag
        
        // Read command options and JOB list
        mpidp.read_table(argc,argv,logout);
        
        eflag = mpidp.master0(argc, argv, nproc);	// = 0 (MPI_Finalize) or 1 (MPI_Abort)
    }
    else {				// for workers
        mpidp.worker(myid,hostname,argc,argv);
    }
    
    etime = MPI_Wtime();

    printf("\nTotal time (process %5d)       = %8.2f sec.\n", myid, etime - stime);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    if( myid == MASTER_PROCESS_ID ) {
        gettimeofday(&et2,NULL);
        const float total_time = (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
        logout << "\nElapsed time  = " << total_time << " sec." << endl;
        printf("\nTotal time (entire process)       = %8.2f sec.\n", total_time);
    }
    
    return 0;
}


//============================================================================//
void Mpidp::read_table(int argc,char *argv[],ofstream &logout)
// Read command options and JOB list
//============================================================================//
{
    int		csize = 0;
    int		ndata = 0;
    string	table;
    
    _Title = "MasterProc";      // TITLE (initialization)
    _Param = "MPIDP";			// PARAM data (initialization)
    _Psize = _Param.size() + 1;
    _Smallest_chunk_size = 1;
    _Largest_chunk_size = LONG_MAX;
    _Out_option = 0;
    //_Ntry = 0;                  // Number retrying limit
    _Worker_life = 3;			// Worker life (default=3)
    
    // for MPIDP options
    for( int i = 1 ; i < argc ; i++ ) {
        if( !strncmp(argv[i],"-tb",3) ) {
            _Table_file = argv[++i];
            logout << "Table file    : -tb " << _Table_file << endl;
        }
        else if( !strncmp(argv[i],"-ch",3) ) {
            _Smallest_chunk_size = atoi(argv[++i]);
            if (_Smallest_chunk_size <= 0) {
                cerr << "[ERROR] Smallest chunk size must be a positive integer!" << endl;
                exit(1);
            }
            logout << "Smallest chunk size    : -ch " << _Smallest_chunk_size << endl;
        }
        else if( !strncmp(argv[i],"-lc",3) ) {
            _Largest_chunk_size = atoi(argv[++i]);
            if (_Largest_chunk_size <= 0) {
                cerr << "[ERROR] Largest chunk size must be a positive integer!" << endl;
                exit(1);
            }
            logout << "Largest chunk size    : -ic " << _Largest_chunk_size << endl;
        }
        else if( !strncmp(argv[i],"-ot",3) ) {
            _Out_option = atoi(argv[++i]);
            logout << "Output option : -ot " << _Out_option << endl;
        }
        else if( !strncmp(argv[i],"-wl",3) ) {
            _Worker_life = atoi(argv[++i]);
            logout << "Worker life   : -wl " << _Worker_life << endl;
        }
        else if( !strncmp(argv[i],"-pg",3) ) {
            logout << "Program name  : -pg " << argv[++i] << endl;
        }
        else if( !strncmp(argv[i],"-lg",3) ) {
            logout << "Log file      : -lg " << argv[++i] << endl;
        }
    }
    if (_Smallest_chunk_size > _Largest_chunk_size) {
        cerr << "[ERROR] Largest chunk size cannot be smaller than smallest chunk size!" << endl;
        exit(1);
    }
    
    // for other(application's) options
    int oflag = 0;
    for( int i = 1 ; i < argc ; i++ ) {
        if( !strncmp(argv[i],"-tb",3) ||
            !strncmp(argv[i],"-ch",3) ||
            !strncmp(argv[i],"-lc",3) ||
            !strncmp(argv[i],"-ot",3) ||
            !strncmp(argv[i],"-rt",3) ||
            !strncmp(argv[i],"-wl",3) ||
            !strncmp(argv[i],"-pg",3) ||
            !strncmp(argv[i],"-lg",3) ) {
            i++;
        }
        else {
            if( oflag == 0 ) {
                logout << "Other options :";
                oflag = 1;
            }
            logout << " " << argv[i];
        }
    }
    
    if( oflag == 1 ) {
        logout << endl;
    }
    logout << endl;
    
    // open JOB list file
    ifstream Input(_Table_file.c_str(),ios::in);
    if( !Input ) {
        cerr << "[ERROR] Table file [" << _Table_file
             << "] was not opened!!" << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
        exit(1);
    }

    _Num_pair = 0;

    // read JOB list file (calculate csize)
    while(1) {
        if( !getline(Input,table) ) break;
        
        table = erase_space(table,7);
        
        if( !strncmp(table.c_str(),"TITLE=",6) ||
            !strncmp(table.c_str(),"Title=",6) ||
            !strncmp(table.c_str(),"title=",6) ) {
            _Title = table.substr(6);
            logout << "TITLE=" << _Title << endl;
        }
        else if( !strncmp(table.c_str(),"PARAM=",6) ||
                 !strncmp(table.c_str(),"Param=",6) ||
                 !strncmp(table.c_str(),"param=",6) ) {
            _Param = table.substr(6);
            _Psize = _Param.size() + 1;
            logout << "PARAM=" << _Param << endl;
        }
        else {
            csize = (table.size() > csize) ? table.size() : csize;
            _Num_pair++;
        }
    }
    
    logout << endl;
    Input.close();
    _Csize = csize + 1;
    
    int line_index = 0;
    _Table_list.reserve((long) _Csize * _Num_pair);
    
    ifstream Input2(_Table_file.c_str(),ios::in);
    if( !Input2 ) {
        cerr << "[ERROR] Table file [" << _Table_file
             << "] was not opened!!" << endl;
        MPI_Abort(MPI_COMM_WORLD,1);
        exit(1);
    }

    // read JOB list file (read lines)
    while(1) {
        if( !getline(Input2,table) ) break;
        
        table = erase_space(table,7);
        
        if( strncmp(table.c_str(),"TITLE=",6) &&
            strncmp(table.c_str(),"Title=",6) &&
            strncmp(table.c_str(),"title=",6) &&
            strncmp(table.c_str(),"PARAM=",6) &&
            strncmp(table.c_str(),"Param=",6) &&
            strncmp(table.c_str(),"param=",6) )
        {
            copy(table.begin(), table.begin() + _Csize, back_inserter(_Table_list));
            _Table_list.back() = 0;
        }
    }
    Input2.close();

    for( int i = 0 ; i < _Csize; i++ ) {
        if( _Table_list[i] == '\t' ) {
            ndata ++;
        }
    }
    
    // According to _Name data
    _Ndata = ndata + 1;
    
    //ntry = _Ntry;
    
    return;
}

//============================================================================//
string Mpidp::erase_space(const string &s0,const int ip)
// Deletion of spaces
//============================================================================//
{
    int 	n;
    string	s1;
    
    s1 = s0;
    
    while(1) {
        n = s1.find(" ");
        
        if( n == std::string::npos ) {
            break;
        }
        else if( n < ip ) {
            s1.erase(n,1);
        }
        else {
            break;
        }
    }
    
    return s1;
}

//============================================================================//
int Mpidp::master0(int argc, char *argv[], const int &nproc)
// master process for NO retry
//============================================================================//
{
    int	wid;				// Worker id
    int		nproc2;
    
#ifdef _OPENMP
    #pragma omp parallel
    {
        nproc2 = omp_get_num_threads();
        if(omp_get_thread_num() == 0) {
            cout << "# Using OpenMP parallelization: " << nproc2 << " threads." << endl;
        }
    }
    //printf("#OpenMP version %d\n", _OPENMP);
#else
    nproc2 = 1;
#endif //#ifdef _OPENMP

    char	*param = new char[_Psize];	// for PARAM data
    
    strcpy(param,_Param.c_str());
    
    // The calculation condition is sent to workers. 
    for( int i = 1 ; i < nproc ; i++ ) {
        MPI_Send(&_Psize,1,MPI_INT,i,100,MPI_COMM_WORLD);
        MPI_Send(param,_Psize,MPI_CHAR,i,200,MPI_COMM_WORLD);
        MPI_Send(&_Csize,1,MPI_INT,i,300,MPI_COMM_WORLD);
        MPI_Send(&_Ndata,1,MPI_INT,i,400,MPI_COMM_WORLD);
        MPI_Send(&_Out_option,1,MPI_INT,i,420,MPI_COMM_WORLD);
    }


    /*
    const long rem = _Num_pair % nproc;
    const long quo = _Num_pair / nproc;
    _Tlist_size = quo * _Csize;
    int begin_index = 0;
    vector<MPI_Request> reqs(nproc - 1);
    vector<MPI_Status> stats(nproc - 1);
    for (int i = 0; i < nproc - 1; i++) {
        const long tlist_size = _Tlist_size + ((i < rem) ? _Csize : 0);
        MPI_Send(&tlist_size, 1, MPI_LONG, i + 1, 430, MPI_COMM_WORLD);
        MPI_Isend(_Table_list.data() + begin_index + _Tlist_size , tlist_size, MPI_CHAR, i + 1, 450, MPI_COMM_WORLD, &reqs[i]);
        begin_index += tlist_size;
    }


    MPI_Waitall(nproc - 1, &reqs[0], &stats[0]);

    _Table_list.resize(_Tlist_size);
    _Table_list.shrink_to_fit();
    */

    // Correction of a bug
    int arglen = max(_Csize,_Psize);
    for( int i = 0 ; i < argc ; i++ ) {
        if( strlen(argv[i]) >= arglen ) {
            arglen = strlen(argv[i]) + 1;
        }
    }
    // Bug fix
    int nwargv = argc + _Ndata*2;
    for( int i = 1 ; i < _Psize-2 ; i++ ) {
        if( param[i] == ' ' ) nwargv++;
    }

    _App = new Application(nproc2);
    _App->initialize();

#pragma omp parallel
    {
        int myid2 = omp_get_thread_num();
        if (myid2 == 0)
            master_thread(nproc);
        else
            worker_thread(nwargv, arglen, argc, argv, MASTER_PROCESS_ID, myid2);
    }
    delete [] param;
    delete _App;

    return 0;
}

//============================================================================//
void Mpidp::worker(int &myid,char *hostname,int argc,char *argv[])
// worker process for function call
//============================================================================//
{
    int		nproc2;
    
#ifdef _OPENMP
    #pragma omp parallel
    {
        nproc2 = omp_get_num_threads();
        if(omp_get_thread_num() == 0) {
            cout << "# Using OpenMP parallelization: " << nproc2 << " threads." << endl;
        }
    }
    //printf("#OpenMP version %d\n", _OPENMP);
#else
    nproc2 = 1;
#endif //#ifdef _OPENMP
    
    MPI_Recv(&_Psize,1,MPI_INT,MASTER_PROCESS_ID,100,MPI_COMM_WORLD,&_Status);
    char *param = new char[_Psize];
    MPI_Recv(param,_Psize,MPI_CHAR,MASTER_THREAD_ID,200,MPI_COMM_WORLD,&_Status);
    _Param = param;
    
    if( !strncmp(param,"MPIDP",5) ) {
        cerr << "[ERROR] [PARAM=] was not found in table file!!" << endl;
        exit(1);
    }
    
    MPI_Recv(&_Csize,1,MPI_INT,MASTER_PROCESS_ID,300,MPI_COMM_WORLD,&_Status);
    MPI_Recv(&_Ndata,1,MPI_INT,MASTER_PROCESS_ID,400,MPI_COMM_WORLD,&_Status);
    MPI_Recv(&_Out_option,1,MPI_INT,MASTER_PROCESS_ID,420,MPI_COMM_WORLD,&_Status);
    
    // Correction of a bug
    int arglen = max(_Csize,_Psize);
    for( int i = 0 ; i < argc ; i++ ) {
        if( strlen(argv[i]) >= arglen ) {
            arglen = strlen(argv[i]) + 1;
        }
    }
    // Bug fix
    int nwargv = argc + _Ndata*2;
    for( int i = 1 ; i < _Psize-2 ; i++ ) {
        if( param[i] == ' ' ) nwargv++;
    }

    /*
    MPI_Recv(&_Tlist_size, 1, MPI_LONG, MASTER_PROCESS_ID, 430, MPI_COMM_WORLD, &_Status);
    _Table_list.resize(_Tlist_size);
    MPI_Recv((void*)_Table_list.data(), _Tlist_size, MPI_CHAR, MASTER_PROCESS_ID, 450, MPI_COMM_WORLD, &_Status);
    */
    _App = new Application(nproc2);
    _App->initialize();

#pragma omp parallel
    {
        int myid2 = omp_get_thread_num();
        worker_thread(nwargv, arglen, argc, argv, myid, myid2);
    }
    
    delete [] param;
    delete _App;
    
    return;
}

//============================================================================//
void Mpidp::master_thread(const int nproc)
//============================================================================//
{
    int countdown = nproc - 1;
    bool end_flag = false;

    while (countdown > 0) {
        int child_id;
        MPI_Recv(&child_id, 1, MPI_INT, MPI_ANY_SOURCE, 500, MPI_COMM_WORLD, &_Status); // receive a request from a child process

#pragma omp critical (distribute0)
        {
            if (!has_next_task() && !get_end_flag()) {
                end_flag_on();
            }
            if (get_end_flag()) {
                long send_size = 1;
                MPI_Send(&send_size, 1, MPI_LONG, child_id, 600, MPI_COMM_WORLD); // send 1
                char end_char = '\0';
                MPI_Send(&end_char, 1, MPI_CHAR, child_id, 700, MPI_COMM_WORLD); // send 0 (end of file message)
                countdown--;
            } else {
                const long table_list_size = _Table_list.size() / _Csize; 
                const long chunk_size = min(_Largest_chunk_size, max(_Smallest_chunk_size, table_list_size / nproc));
                long send_size = min(chunk_size * _Csize, (long) _Table_list.size());
                MPI_Send(&send_size, 1, MPI_LONG, child_id, 600, MPI_COMM_WORLD); // send number of tasks to send
                MPI_Send(_Table_list.data() + _Table_list.size() - send_size, send_size, MPI_CHAR, child_id, 700, MPI_COMM_WORLD); // send tasks
                _Table_list.resize(_Table_list.size() - send_size);
                _Table_list.shrink_to_fit();
            }
        }
    }
}

//============================================================================//
void Mpidp::worker_thread(const int nwargv, const int arglen, const int argc, char *argv[], const int myid, const int myid2)
//============================================================================//
{
    struct timeval et3, et4;
    int		wargc;		// # of application command line parameters
    char **wargv;

#pragma omp critical (alloc)
    {
        wargv = new char*[nwargv];
        for( int i = 0 ; i < nwargv ; i++ ) {
            wargv[i] = new char[arglen];
        }
    }
    int argc2 = argument(argc,argv,wargv); // # of mpidp comand line parameters

    gettimeofday(&et3,NULL);

    while(1) {

#pragma omp critical (distribute0)
        {
            if (!has_next_task() && !get_end_flag()) {
                if (myid == 0) {
                    end_flag_on();
                } else {
                    request_tasks(myid);
                }
            }
            if (!get_end_flag()) {
                char ctable[_Csize];
                strncpy(ctable, _Table_list.data() + _Table_list.size() - _Csize, _Csize);

                wargc = for_worker(ctable,argc2,wargv);
                _Table_list.resize(_Table_list.size() - _Csize);
                _Table_list.shrink_to_fit();
            }
        }
        if (get_end_flag()) break;

        try {
            throw _App->application(wargc,wargv,myid2);	// application's main function
        }
        catch(int e) {
            if (e) {
                cerr << "_App->application was not successfully finished" << endl;
                MPI_Abort(MPI_COMM_WORLD,1);
            }
        }
        catch(char *e) {
            cerr << "[ERROR] [application] exception : " << e << endl;
            MPI_Abort(MPI_COMM_WORLD,1);
        }

    }
    delete [] wargv;

    gettimeofday(&et4,NULL);
    const float total_time = (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));
    printf("\nTotal time ([proc / thread] %5d / %2d)       = %8.2f sec.\n", myid, myid2, total_time);
}


//============================================================================//
int Mpidp::argument(int argc,char *argv[],char **wargv)
// Procedure of options for function call version
//============================================================================//
{
    int ic = 0;
    
    for( int i = 0 ; i < argc ; i++ ) {
        if( !strncmp(argv[i],"-tb",3) ||
            !strncmp(argv[i],"-ch",3) ||
            !strncmp(argv[i],"-lc",3) ||
            !strncmp(argv[i],"-ot",3) ||
            !strncmp(argv[i],"-rt",3) ||
            !strncmp(argv[i],"-wl",3) ||
            !strncmp(argv[i],"-lg",3) ) {
            i++;
        }
        else {
            strcpy(wargv[ic++],argv[i]);
        }
    }
    
    return ic;
}

//============================================================================//
int Mpidp::argument(int argc,char *argv[],string &main_argv)
// Procedure of options for system call version
//============================================================================//
{
    int iflag = 1;
    
    for( int i = 1 ; i < argc ; i++ ) {
        if( !strncmp(argv[i],"-pg",3) ) {
            main_argv = argv[++i];
            iflag = 0;
        }
        else if( !strncmp(argv[i],"-tb",3) ||
                 !strncmp(argv[i],"-ch",3) ||
                 !strncmp(argv[i],"-lc",3) ||
                 !strncmp(argv[i],"-ot",3) ||
                 !strncmp(argv[i],"-rt",3) ||
                 !strncmp(argv[i],"-wl",3) ||
                 !strncmp(argv[i],"-lg",3) ) {
            i++;
        }
        else {
            main_argv += ' ';
            main_argv += argv[i];
        }
    }
    
    return iflag;
}

//============================================================================//
int Mpidp::for_worker(char *ctable,int argc2,char **wargv)
// Preparation using function call version by workers
//============================================================================//
{
    string	position;
    string	tstock;
    string	sparam = _Param;
    char		tag[5], *elem;
    
    //strcpy(_Name,strtok(ctable,"\t"));
    position = "$1";
    tstock = strtok(ctable,"\t");
    sparam = replace_pattern(sparam,position,tstock);
    
    for( int i = 1 ; i < _Ndata ; i++ ) {
        sprintf(tag,"$%d",i+1);
        position = tag;
        tstock = strtok(NULL,"\t");
        
        // character is replaced
        sparam = replace_pattern(sparam,position,tstock);
    }
    
    char param[sparam.size()+1];
    //char *param = new char[sparam.size()+1];
    strcpy(param,sparam.c_str());
    strcpy(wargv[argc2++],strtok(param," "));
    
    while( (elem = strtok(NULL," ")) ) {
        strcpy(wargv[argc2++],elem);
    }
    
    return argc2;
}


//============================================================================//
string Mpidp::replace_pattern(const string &pattern,const string &position,
                              const string &option)
// the pattern of a character is replaced 
//============================================================================//
{
    string	command;
    int	pos_before = 0;
    int	pos = 0;
    int	len = position.size();
    
    while( (pos = pattern.find(position, pos)) != std::string::npos) {
        command.append(pattern, pos_before, pos - pos_before);
        command.append(option);
        pos += len;
        pos_before = pos;
    }
    
    command.append(pattern, pos_before, pattern.size() - pos_before);
    
    return command;
}

//============================================================================//
void Mpidp::request_tasks(const int &myid)
//============================================================================//
{
    MPI_Send(&myid, 1, MPI_INT, 0, 500, MPI_COMM_WORLD); // request num_of_tasks_per_request tasks
    int task_msg_size;
    MPI_Recv(&task_msg_size, 1, MPI_LONG, 0, 600, MPI_COMM_WORLD, &_Status); // receive task_msg_size
    _Table_list.resize(task_msg_size);
    MPI_Recv(_Table_list.data(), task_msg_size, MPI_CHAR, 0, 700, MPI_COMM_WORLD, &_Status); // receive tasks
    if (_Table_list[0] == 0) { // if received end of file message
        end_flag_on();
    }
}
