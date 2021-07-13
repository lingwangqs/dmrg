#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include "dmrg_su2_mpi.hpp"

using namespace std;

extern "C"{
  void dstev_(char*,int*,double*,double*,double*,int*,double*,int*);
  double ran_();
}

extern dmrg_su2* chain_ptr;
extern int comm_rank;
//--------------------------------------------------------------------------------------
lanczos_su2::lanczos_su2(){
//--------------------------------------------------------------------------------------
}

//--------------------------------------------------------------------------------------
lanczos_su2::~lanczos_su2(){
//--------------------------------------------------------------------------------------
  delete []eig;
  delete []vec;
  delete []nnl;
  delete []aal;
  delete []ff;
}

//--------------------------------------------------------------------------------------
void lanczos_su2::initialize_lanczos(tensor_su2& vin,int ml){
//--------------------------------------------------------------------------------------
  int ne1,ne2;
  vin.get_nelement(ne1,ne2);
  if(ne1>ml)
    mlanc=ml;
  else
    mlanc=ne1+1;
  aal=new double[mlanc];
  nnl=new double[mlanc+1];
  ff=new tensor_su2[mlanc+1];
  eig=new double[mlanc];
  vec=new double[mlanc*mlanc];
  ff[0]=vin;
  ff[0].normalize_vector();
}

//--------------------------------------------------------------------------------------
void lanczos_su2::lanczos1(int il, int ir, tensor_su2& vout,int minlan){
//--------------------------------------------------------------------------------------
  int m,n,i,j,flag;
  double q1,q2;
  tensor_su2 tmp;
  for(i=0;i<mlanc;i++)
    aal[i]=0;
  for(i=0;i<mlanc+1;i++)
    nnl[i]=0;
  nnl[0]=1;
  diag_op(il,ir,ff[0],ff[1]);
  aal[0]=ff[0].inner_prod(ff[1]);
  tmp=ff[0];
  tmp*=aal[0];
  ff[1]-=tmp;
  nnl[1]=sqrt(ff[1].inner_prod(ff[1]));
  ff[1]/=nnl[1];
  if(nnl[1]<1.e-8){
    //cout<<"lanczos1 exit at m=1"<<endl;
    diatridiag(1);
    compute_eigenvector(il,ir,vout,1,n,flag);
    return;
  }
  n=4;
  for(m=2;m<mlanc;m++){
    //cout<<"lanczos1 m="<<m<<endl;
    diag_op(il,ir,ff[m-1],ff[m]);
    aal[m-1]=ff[m-1].inner_prod(ff[m]);
    tmp=ff[m-1];
    tmp*=aal[m-1];
    ff[m]-=tmp;
    tmp=ff[m-2];
    tmp*=nnl[m-1];
    ff[m]-=tmp;
    nnl[m]=sqrt(ff[m].inner_prod(ff[m]));
    ff[m]/=nnl[m];
    //cout<<"aal["<<m-1<<"]="<<aal[m-1]<<"\tnnl["<<m<<"]="<<nnl[m]<<endl;
    //orthogonalize ff[m] with all previous ff[j]
    for(j=0;j<m;j++){
      q1=ff[m].inner_prod(ff[j]);
      q2=1./sqrt(1.-q1*q1);
      //cout<<"q1="<<q1<<"\tq2="<<q2<<endl;
      tmp=ff[j];
      tmp*=q1;
      ff[m]-=tmp;
      ff[m]*=q2;
      //ff[m]=(ff[m]-ff[j]*q1)*q2;
    }
    if(nnl[m]<1.e-8){
      diatridiag(m);
      compute_eigenvector(il,ir,vout,m,n,flag);
      break;
    }
    else{
      if(mlanc>m&&m<n)continue;
      diatridiag(m);
      compute_eigenvector(il,ir,vout,m,n,flag);
      if(flag)break;
    }
  }
}

//--------------------------------------------------------------------------------------
void lanczos_su2::diatridiag(int n){
//--------------------------------------------------------------------------------------
  double *d,*e,*work;
  int i,j,info;
  ofstream fout;
  char jobz;
  jobz='V';
  d=new double[n];
  e=new double[n];
  work=new double[2*n];
  for(i=0;i<n;i++){
    d[i]=aal[i];
    e[i]=nnl[i+1];
  }
  dstev_(&jobz,&n,d,e,vec,&n,work,&info);
  for(i=0;i<n;i++){
    eig[i]=d[i];
    //if(comm_rank==0)cout<<"eig["<<i<<"]="<<eig[i]<<endl;
  }
  //fout.open("out.dat",ios::app);
  //fout<<n<<"\t"<<setprecision(8)<<eig[0]<<"\t"<<eig[1]<<"\t"<<eig[2]<<endl;
  //fout.close();
  delete []d;
  delete []e;
  delete []work;
}

//--------------------------------------------------------------------------------------
void lanczos_su2::compute_eigenvector(int il, int ir, tensor_su2& vout, int mlanc_curr, int& n, int& stp){
//--------------------------------------------------------------------------------------
  ofstream fout;
  double *eigval,*olap;
  int i,j,k;
  char name[100],rank[10];
  tensor_su2 tmp;
  eigval=new double[mlanc];
  olap=new double[mlanc];
  for(k=0;k<1;k++){//only compute the gs vector
    for(j=0;j<mlanc_curr;j++)
      if(j==0){
	vout=ff[j];
	vout*=vec[j+k*mlanc_curr];
      }
      else{
	tmp=ff[j];
	tmp*=vec[j+k*mlanc_curr];
	vout+=tmp;
      }
    check_eigenvector(il,ir,vout,eigval[k],olap[k]);
    if(olap[k]>1-1.e-10||mlanc_curr==mlanc-1){
      stp=1;
      strcpy(name,"out-");
      sprintf(rank,"%d",comm_rank);
      strcat(name,rank);
      strcat(name,".dat");
      fout.open(name,ios::app);
      fout<<il<<"\t"<<ir<<"\t"<<setprecision(12)<<"\t"<<mlanc_curr<<"\t"<<-eigval[k]<<"\t"<<olap[k]<<endl;
      fout.close();
    }
    else{
      if(olap[k]<1-5.e-9)n=mlanc_curr+5;
      stp=0;
    }
  }
  delete []eigval;
  delete []olap;
}

//--------------------------------------------------------------------------------------
void lanczos_su2::check_eigenvector(int il, int ir, tensor_su2& vout, double& eigval, double& olap){
//--------------------------------------------------------------------------------------
  tensor_su2 vec;
  double prod,overlap,nor;
  nor=sqrt(vout.inner_prod(vout));
  vout/=nor;
  diag_op(il,ir,vout,vec);
  eigval=vec.inner_prod(vout);
  prod=sqrt(vec.inner_prod(vec));
  overlap=eigval/prod;
  olap=fabs(overlap);
}

//------------------------------------------------------------------------------
void lanczos_su2::diag_op(int i0,int i1, tensor_su2& v1, tensor_su2& v2){
//------------------------------------------------------------------------------
  chain_ptr->hamiltonian_vector_multiplication_idmrg(i0,i1,v1,v2);
}
