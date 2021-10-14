#ifdef MI_MALLOC
#include "mimalloc-new-delete.h"
#endif 

#include <mpi.h>
#include <omp.h>
#include <mkl.h>
#include "dmrg_su2_mpi.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <sys/types.h>
#include <time.h>
#include <cstring>



using namespace std;
int max_dcut=8000,myrank=0,psize=12,mkl=1,comm_rank,comm_size,thread,bdry=0,jobid,memory_flag=0;
double delta=1,qdelta;
double t1,t2,preparetime=0,svdtime=0,lanczostime=0;
dmrg_su2* chain_ptr;
tensor ***spin_op,***cgc_coef_left,***cgc_coef_rght,*cgc_coef_singlet,*identity;
double **spin_op_trace,****fac_operator_onsite_left,****fac_operator_onsite_rght,*****fac_operator_transformation_left,*****fac_operator_transformation_rght,****fac_operator_pairup_left,****fac_operator_pairup_rght,***fac_hamilt_vec,*****fac_permutation_left,*****fac_permutation_rght;


extern "C"{
  double ran_();
  void initran_(int*);  
}

int main(int argc,char** argv){
  int i,sec,read,read2,d,dr,len,exci,niter,lx,ly,nstep;
  bool pass;
  double t3, t4,enr1,enr2;
  ofstream fout;
  char name[100],id[10],dk[10];

  MPI::Init(argc,argv);
  comm_rank=MPI::COMM_WORLD.Get_rank();
  comm_size=MPI::COMM_WORLD.Get_size();
  cout<<"my comm_rank="<<comm_rank<<endl;
  if(comm_rank==0)
    cout<<"my comm_size="<<comm_size<<endl;

  t1=clock();
  i=(long)t1;
  initran_(&i);
  for(i=1;i<argc;i++){
    if(strcmp(argv[i],"-d")==0){
      max_dcut=atoi(argv[i+1]);
      if(comm_rank==0)
	cout<<"max_dcut="<<max_dcut<<endl;
    }
    else if(strcmp(argv[i],"-dr")==0){
      read=atoi(argv[i+1]);
      if(comm_rank==0)
	cout<<"read="<<read<<endl;
    }
    else if(strcmp(argv[i],"-dge")==0){
      read2=atoi(argv[i+1]);
      if(comm_rank==0)
	cout<<"read2="<<read2<<endl;
    }
    else if(strcmp(argv[i],"-boundary")==0){
      bdry=atoi(argv[i+1]);
      if(comm_rank==0)
	cout<<"boundary="<<bdry<<endl;
    }
    else if(strcmp(argv[i],"-lx")==0){
      lx=atoi(argv[i+1]);
      if(comm_rank==0)
	cout<<"lx="<<lx<<endl;
    }
    else if(strcmp(argv[i],"-ly")==0){
      ly=atoi(argv[i+1]);
      if(comm_rank==0)
	cout<<"ly="<<ly<<endl;
    }
    else if(strcmp(argv[i],"-sec")==0){
      sec=atoi(argv[i+1]);
      if(comm_rank==0)
	cout<<"sec="<<sec<<endl;
    }
    else if(strcmp(argv[i],"-exci")==0){
      exci=atoi(argv[i+1]);
      if(comm_rank==0)
	cout<<"exci="<<exci<<endl;
    }
    else if(strcmp(argv[i],"-niter")==0){
      niter=atoi(argv[i+1]);
      if(comm_rank==0)
	cout<<"niter="<<niter<<endl;
    }
    else if(strcmp(argv[i],"-jcoup")==0){
      delta=atof(argv[i+1]);
      if(comm_rank==0)
	cout<<"jcoup="<<delta<<endl;
    }
    else if(strcmp(argv[i],"-qcoup")==0){
      qdelta=atof(argv[i+1]);
      if(comm_rank==0)
	cout<<"qcoup="<<qdelta<<endl;
    }
    else if(strcmp(argv[i],"-jobid")==0){
      jobid=atoi(argv[i+1]);
      if(comm_rank==0)
	cout<<"jobid="<<jobid<<endl;
    }
    else if(strcmp(argv[i],"-memory_flag")==0){
      memory_flag=atoi(argv[i+1]);
      if(comm_rank==0)
	cout<<"memory_flag="<<memory_flag<<endl;
    }
    else if(strcmp(argv[i],"-psize")==0){
      psize=atoi(argv[i+1]);
      if(comm_rank==0)
	cout<<"psize="<<psize<<endl;
    }
    else if(strcmp(argv[i],"-mkl")==0){
      mkl=atoi(argv[i+1]);
      if(comm_rank==0)
	cout<<"mkl="<<mkl<<endl;
    }
  }

  sprintf(dk,"%d",comm_rank+3);
  strcpy(name,"mkdir ./data");
  strcat(name,dk);
  //strcat(name,"/");
  //sprintf(id,"%d",jobid);
  //strcat(name,id);
  system(name);

  omp_set_num_threads(psize);
  mkl_set_num_threads(mkl);
  // This is needed to allow MKL to run in parallel 
  omp_set_nested(1);
 // omp_set_max_active_levels(2);
  mkl_set_dynamic(0);

  dmrg_su2 chain(lx,ly,sec,read,read2,exci);
  chain_ptr=&chain;
  MPI::COMM_WORLD.Barrier();
  if(read==0&&exci==0){
    chain.do_idmrg();
  }
  else{
    pass=chain.read_mps(read,read2);
    MPI::COMM_WORLD.Barrier();
    if(comm_rank==0){
      cout<<comm_rank<<" pass="<<pass<<endl;
      chain.read_enr(read,read2);
      chain.read_ww(read,read2);
    }
    MPI::COMM_WORLD.Barrier();
    if(comm_rank==1){
      cout<<comm_rank<<" pass="<<pass<<endl;
      chain.read_enr(read,read2);
      chain.read_ww(read,read2);
    }
    MPI::COMM_WORLD.Barrier();
    t3=omp_get_wtime();
    if(pass) chain.prepare_sweep();
    else chain.do_idmrg();
    t4=omp_get_wtime();
    if(comm_rank==0)cout<<"presweep time="<<t4-t3<<endl;
  }
  enr1=chain.get_enr();
  for(i=0;i<niter;i++){
    t3=omp_get_wtime();
    chain.sweep();
    t4=omp_get_wtime();
    if(comm_rank==0){
      fout.open("time.dat",ios::app);
      if(fout.is_open()){
	fout<<jobid<<"\tsweep="<<i<<"\tmemory_flag="<<memory_flag<<"\ttime="<<t4-t3<<endl;
	fout.close();
      }
    }
    enr2=chain.get_enr();
    if(comm_rank==0){
      MPI::COMM_WORLD.Send(&enr2,1,MPI_DOUBLE,1,0);
      MPI::COMM_WORLD.Send(&enr1,1,MPI_DOUBLE,1,0);
    }
    else if(comm_rank==1){
      MPI::COMM_WORLD.Recv(&enr2,1,MPI_DOUBLE,0,0);
      MPI::COMM_WORLD.Recv(&enr1,1,MPI_DOUBLE,0,0);
    }
    if(fabs(enr1-enr2)<1.e-6)
      break;
    enr1=enr2;
  }

  sprintf(dk,"%d",comm_rank+3);
  strcpy(name,"rm -rf ./data");
  strcat(name,dk);
  //strcat(name,"/");
  //sprintf(id,"%d",jobid);
  //strcat(name,id);
  system(name);

  MPI::Finalize();
}
