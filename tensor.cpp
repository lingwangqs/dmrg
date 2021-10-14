#include <omp.h>
#include "tensor.hpp"
#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision
#include <stdlib.h>
#include <string.h>
#include <math.h>

 
using namespace std;

extern "C"{
  void dgemm_(char*,char*,int*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);
  void dsyev_(char*,char*,int*,double*,int*,double*,double*,int*,int*);
  void dsyevx_(char*,char*,char*,int*,double*,int*,double*,double*,int*,int*,double*,int*,double*,double*,int*,double*,int*,int*,int*,int*);
  void dsyrk_(char*,char*,int*,int*,double*,double*,int*,double*,double*,int*);
  void dposv_(char*,int*,int*,double*,int*,double*,int*,int*);
  double ran_();
}
void obtain_symmetric_matrix_eigenvector_2(double*,double*,int,int);
void obtain_symmetric_matrix_eigenvector(double*,double*,int);
void get_tensor_index(int&,int,int*,int*);
void get_bond_index(int,int,int*,int*);
double factorial(int);

extern int max_dcut,myrank,psize;
double tolerance=1.e-8;
//--------------------------------------------------------------------------------------
tensor::tensor(){
//--------------------------------------------------------------------------------------
  telement=NULL;
  bonddim=NULL;
  nelement=0;
  nbond=0;
}

//--------------------------------------------------------------------------------------
tensor::tensor(int nb,int *bdim){
//--------------------------------------------------------------------------------------
  int i,j,k,l;
  nbond=nb;
  bonddim=new int[nbond];
  nelement=1;
  for(i=0;i<nbond;i++){
    bonddim[i]=bdim[i];
    nelement*=bonddim[i];
  }
  telement=new double[nelement];
  for(i=0;i<nelement;i++)
    telement[i]=ran_()*2.-1.;
}

//--------------------------------------------------------------------------------------
tensor::tensor(int nb, int *bdim, double *tele){
//--------------------------------------------------------------------------------------
  int i,j,k,label,*dd;
  nbond=nb;
  bonddim=new int[nbond];
  nelement=1;
  for(i=0;i<nbond;i++){
    bonddim[i]=bdim[i];
    nelement*=bonddim[i];
  }
  telement=new double[nelement];
  //#pragma omp parallel for default(shared) private(i) 
  memcpy(telement,tele,sizeof(double) * nelement) ;
}

//--------------------------------------------------------------------------------------
void tensor::copy(tensor& t1){
//--------------------------------------------------------------------------------------
  int nb;
  nb=t1.get_nbond();
  copy(nb,t1.get_bonddim_ptr(),t1.getptr());
}

//--------------------------------------------------------------------------------------
void tensor::copy(int nb,int *bdim,double* tele){
//--------------------------------------------------------------------------------------
  int i;
  nbond=nb;
  if(bonddim!=NULL)delete []bonddim;
  bonddim=new int[nbond];
  nelement=1;
  for(i=0;i<nbond;i++){
    bonddim[i]=bdim[i];
    nelement*=bonddim[i];
  }
  if(telement!=NULL)delete []telement;
  telement=new double[nelement];
  //#pragma omp parallel for default(shared) private(i) 
  for(i=0;i<nelement;i++)
    telement[i]=tele[i];
}

//--------------------------------------------------------------------------------------
void tensor::shift_copy(int d1, int d2, int shift0, int shift2, int shift1, int shift3, tensor& t1 ){
  //--------------------------------------------------------------------------------------
  int i,j,k,l;
  for(i=0;i<d1;i++)
    for(j=0;j<d2;j++){
      k=shift0+i+(shift2+j)*bonddim[0];
      l=shift1+i+(shift3+j)*t1.get_bonddim(0);
      telement[k]=t1.get_telement(l);
    }
}

//--------------------------------------------------------------------------------------
void tensor::shift_copy(int ind, int dim, int shift, tensor& t1){
//--------------------------------------------------------------------------------------
  int i,j,i0,nele,*bdim,*aa;
  clean();
  nbond=t1.get_nbond();
  nele=t1.get_nelement();
  bonddim=new int[nbond];
  bdim=new int[nbond];
  aa=new int[nbond];
  nelement=1;
  for(i=0;i<nbond;i++){
    bdim[i]=t1.get_bonddim(i);
    bonddim[i]=bdim[i];
    if(i==ind){
      bonddim[i]=dim;
      if(dim<bdim[i]+shift){
	cout<<"tensor::shift_copy, something is wrong shift_copy"<<endl;
	cout<<"dim="<<dim<<"\tshift="<<shift<<"\tbdim="<<bdim[i]<<endl;
	exit(0);
      }
    }
    nelement*=bonddim[i];
  }
  telement=new double[nelement];
  for(i=0;i<nelement;i++)
    telement[i]=0;
  for(i=0;i<nele;i++){
    i0=i;
    for(j=0;j<nbond;j++){
      aa[j]=i0%bdim[j];
      i0/=bdim[j];
      if(j==ind)aa[j]+=shift;
    }
    i0=0;
    for(j=nbond-1;j>=0;j--){
      i0*=bonddim[j];
      i0+=aa[j];
    }
    telement[i0]=t1.get_telement(i);
  }
  delete []bdim;
  delete []aa;
}

//--------------------------------------------------------------------------------------
void tensor::random_init(){
//--------------------------------------------------------------------------------------
  int i;
  for(i=0;i<nelement;i++)
    telement[i]=(ran_()*2.-1.);
}

//--------------------------------------------------------------------------------------
void tensor::alloc_space(int nb,int *bdim){
//--------------------------------------------------------------------------------------
  int i;
  nbond=nb;
  if(bonddim!=NULL)delete []bonddim;
  bonddim=new int[nbond];
  nelement=1;
  for(i=0;i<nbond;i++){
    bonddim[i]=bdim[i];
    nelement*=bonddim[i];
  }
  if(telement!=NULL) delete []telement;
  telement=new double[nelement];
  for(i=0;i<nelement;i++)
    telement[i]=0;
}

//--------------------------------------------------------------------------------------
void tensor::clean(){
//--------------------------------------------------------------------------------------
  if(telement!=NULL)  delete []telement;
  if(bonddim!=NULL) delete []bonddim;
  telement=NULL;
  bonddim=NULL;
  nelement=0;
  nbond=0;
}

//--------------------------------------------------------------------------------------
tensor::~tensor(){
//--------------------------------------------------------------------------------------
  if(telement!=NULL)  delete []telement;
  if(bonddim!=NULL) delete []bonddim;
}

//--------------------------------------------------------------------------------------
void tensor::obtaintelement(double* arr,int istart,int len){
//--------------------------------------------------------------------------------------
  int i;
  for(i=0;i<len;i++)arr[i]=telement[istart+i];
}

//--------------------------------------------------------------------------------------
void tensor::obtaintelement(double* arr){
//--------------------------------------------------------------------------------------
  int i;
  for(i=0;i<nelement;i++)
    arr[i]=telement[i];
}

//--------------------------------------------------------------------------------------
double tensor::get_norm(){
//--------------------------------------------------------------------------------------
  int i;
  double nor;
  nor=0;
  for(i=0;i<nelement;i++)
    nor+=telement[i]*telement[i];
  return nor;
}

//--------------------------------------------------------------------------------------
tensor& tensor::operator = (const tensor& t1){
//--------------------------------------------------------------------------------------
  int i,j;
  clean();
  if(t1.get_nbond()==0)return *this;
  nbond=t1.get_nbond();
  nelement=t1.get_nelement();
  bonddim=new int[nbond];
  telement=new double[nelement];
  for(i=0;i<nbond;i++)
    bonddim[i]=t1.get_bonddim(i);
  //#pragma omp parallel for default(shared) private(i) 
  memcpy(telement, t1.telement, sizeof(double) * nelement) ;
  //for(i=0;i<nelement;i++) {
  //  telement[i]=t1.get_telement(i);
  //}
  return *this;
}

//--------------------------------------------------------------------------------------
bool tensor::operator != (const tensor& t1){
//--------------------------------------------------------------------------------------
  int i;

  if(t1.get_nbond()!=nbond)return true;
  if(t1.get_nelement()!=nelement)return true;
  for(i=0;i<nelement;i++)
    //if(fabs(fabs(t1.get_telement(i))-fabs(telement[i]))>tolerance)return true;
    if(fabs(t1.telement[i]-telement[i])>tolerance)return true;
  return false;
}

//--------------------------------------------------------------------------------------
bool tensor::operator == (const tensor& t1){
//--------------------------------------------------------------------------------------
  if(*this!=t1)return false;
  else return true;
}

//--------------------------------------------------------------------------------------
tensor& tensor::operator = (double zero){
//--------------------------------------------------------------------------------------
  int i;
  for(i=0;i<nelement;i++) {
    telement[i]=zero;
  }
  return *this;
}

//--------------------------------------------------------------------------------------
tensor& tensor::operator += (const tensor& t1 ){
//--------------------------------------------------------------------------------------
  int i;
  tensor t2;
  if(nelement!=t1.get_nelement()){
    cout<<"+= nelement not consistent"<<endl;
    print();
    t2=t1;
    t2.print();
    exit(0);
  }
  if(t1.get_nbond()==0) {
    return *this;
  }
  for(i=0;i<nelement;i++) {
    telement[i] += t1.telement[i];
  }
  return *this;
}

//--------------------------------------------------------------------------------------
tensor& tensor::operator -= (const tensor& t1 ){
//--------------------------------------------------------------------------------------
  int i;
  if(nelement!=t1.get_nelement()){
    cout<<"+= nelement not consistent"<<endl;
    exit(0);
  }
  if(t1.get_nbond()==0)return *this;
  
  for(i=0;i<nelement;i++)
    telement[i] -= t1.telement[i];
  return *this;
}

//--------------------------------------------------------------------------------------
tensor tensor::operator * (double aa){
//--------------------------------------------------------------------------------------
  int i;
  tensor t2;
  if(nelement==0)return t2;
  t2.nelement=nelement;
  t2.telement=new double[nelement];
  t2.nbond=nbond;
  t2.bonddim=new int[nbond];
  for(i=0;i<nbond;i++)
    t2.bonddim[i]=bonddim[i];
  for(i=0;i<nelement;i++)
    t2.telement[i]=telement[i]*aa;
  return t2;
}

//--------------------------------------------------------------------------------------
tensor tensor::operator / (double aa){
//--------------------------------------------------------------------------------------
  int i;
  tensor t2;
  if(nelement==0)return t2;
  t2.nelement=nelement;
  t2.telement=new double[nelement];
  t2.nbond=nbond;
  t2.bonddim=new int[nbond];
  aa = 1.0 / aa ;
  for(i=0;i<nbond;i++)
    t2.bonddim[i]=bonddim[i];
  for(i=0;i<nelement;i++)
    t2.telement[i]=telement[i] * aa;
  return t2;
}

//--------------------------------------------------------------------------------------
tensor tensor::operator + (const tensor& t1){
//--------------------------------------------------------------------------------------
  int i;
  tensor t2;
  if(nelement==0)return t2;
  t2.nelement=nelement;
  t2.telement=new double[nelement];
  t2.nbond=nbond;
  t2.bonddim=new int[nbond];
  for(i=0;i<nbond;i++)
    t2.bonddim[i]=bonddim[i];
  for(i=0;i<nelement;i++)
    t2.telement[i]=telement[i]+t1.telement[i];
  return t2;
}

//--------------------------------------------------------------------------------------
tensor tensor::operator - (const tensor& t1){
//--------------------------------------------------------------------------------------
  int i;
  tensor t2;
  if(nelement==0)return t2;
  t2.nelement=nelement;
  t2.telement=new double[nelement];
  t2.nbond=nbond;
  t2.bonddim=new int[nbond];
  for(i=0;i<nbond;i++)
    t2.bonddim[i]=bonddim[i];
  for(i=0;i<nelement;i++)
    t2.telement[i]=telement[i]-t1.telement[i];
  return t2;
}

//--------------------------------------------------------------------------------------
tensor& tensor::operator *= (double c){
//--------------------------------------------------------------------------------------
  int i;
  //#pragma omp parallel for default(shared) private(i) 
  //for(i=0;i<nelement;i++)telement[i]*=c;
  for(i=0;i<nelement;i++) {
    telement[i]*=c;
  }
  return *this;
}

//--------------------------------------------------------------------------------------
tensor& tensor::operator /= (double c){
//--------------------------------------------------------------------------------------
  int i ;
  //#pragma omp parallel for default(shared) private(i) 
  c = 1.0 / c ;
  for(i=0;i<nelement;i++) {
    telement[i]*=c;
  }
  return *this;
}

//--------------------------------------------------------------------------------------
void tensor::print(){
//--------------------------------------------------------------------------------------
  int i;
  for(i=0;i<nbond;i++)
    cout<<i<<"\t"<<bonddim[i]<<endl;
  cout<<"nelement="<<nelement<<endl;
  cout<<"nbond="<<nbond<<endl;
  ///*
  cout<<"==============="<<endl;
  for(i=0;i<nelement;i++)
    if(fabs(telement[i])>tolerance)
      cout<<i<<"\t"<<setprecision(20)<<telement[i]<<endl;
  cout<<"==============="<<endl;
  //*/

}

 
//--------------------------------------------------------------------------------------
void tensor::exchangeindex(int ind1,int ind2){
//--------------------------------------------------------------------------------------
  int i,j,*bdim,**aa;
  double *tele;
  if(ind1==ind2)return;
  aa=new int*[psize];
  bdim=new int[nbond];
  for(i=0;i<psize;i++)
    aa[i]=new int[nbond];
  for(i=0;i<nbond;i++)
    bdim[i]=bonddim[i];
  i=bdim[ind1];
  bdim[ind1]=bdim[ind2];
  bdim[ind2]=i;
  tele=new double[nelement];
  //#pragma omp parallel for default(shared) private(i,j,myrank) schedule(static,1)
  for(i=0;i<nelement;i++){
    //myrank=omp_get_thread_num();
    myrank=0;
    get_bond_index(i,nbond,bonddim,aa[myrank]);
    j=aa[myrank][ind1];
    aa[myrank][ind1]=aa[myrank][ind2];
    aa[myrank][ind2]=j;
    get_tensor_index(j,nbond,bdim,aa[myrank]);
    tele[j]=telement[i];
  }
  delete []telement;
  telement=tele;
  delete []bonddim;
  bonddim=bdim;
  for(i=0;i<psize;i++)
    delete []aa[i];
  delete []aa;
}  



//--------------------------------------------------------------------------------------
void tensor::shift(int i0,int i1){
//--------------------------------------------------------------------------------------
  int i,j,ishift,*bdim,*aa,*bb;
  double *tele;
  if(i0==i1)return;
  bdim=new int[nbond];
  if(i1>i0)ishift=i1-i0;
  else ishift=nbond-(i0-i1);
  for(i=0;i<nbond;i++)
    bdim[(i+ishift)%nbond]=bonddim[i];
  aa=new int[nbond];
  bb=new int[nbond];
  tele=new double[nelement];
  for(i=0;i<nelement;i++){
    get_bond_index(i,nbond,bonddim,aa);
    for(j=0;j<nbond;j++)
      bb[(j+ishift)%nbond]=aa[j];
    get_tensor_index(j,nbond,bdim,bb);
    tele[j]=telement[i];
  }
  delete []telement;
  telement=tele;
  delete []bonddim;
  bonddim=bdim;
  delete []aa;
  delete []bb;
}

//--------------------------------------------------------------------------------------
void tensor::mergeindex(int ind1,int ind2){
//--------------------------------------------------------------------------------------
  //allow to merge only when ind2=ind1+1
  int i,j;
  if(ind2!=ind1+1){
    cout<<"merge index with ind1="<<ind1<<" ind2="<<ind2<<endl;
    exit(0);
  }
  bonddim[ind1]*=bonddim[ind2];
  nbond--;
  for(i=ind2;i<nbond;i++)
    bonddim[i]=bonddim[i+1];
}

//--------------------------------------------------------------------------------------
void tensor::separateindex(int ind, int dim1, int dim2){
//--------------------------------------------------------------------------------------
  int *bdim,i,j;
  if(dim1*dim2!=bonddim[ind]){
    cout<<"separateindex: the bonddim are not consistent "<<bonddim[ind]<<" "<<dim1<<" "<<dim2<<endl;
    exit(0);
  }
  bdim=new int[nbond+1];
  for(i=0;i<ind;i++)
    bdim[i]=bonddim[i];
  for(i=ind+1;i<nbond;i++)
    bdim[i+1]=bonddim[i];
  nbond++;
  bdim[ind]=dim1;
  bdim[ind+1]=dim2;
  delete []bonddim;
  bonddim=bdim;
}

//--------------------------------------------------------------------------------------
tensor& tensor::contract(tensor& t1, int i1, tensor& t2, int i2){
//--------------------------------------------------------------------------------------
  if(t1.bonddim[i1]!=t2.bonddim[i2]){
    cout<<"tensor::contract the bond dimension is not consistent"<<endl;
    cout<<t1.bonddim[i1]<<"\t"<<t2.bonddim[i2]<<endl;
    exit(0);
  }
  tensor tp1,tp2;
  char transa,transb;
  int lda,ldb,nele,nbond1,*bdim;
  int m,n,k,i,j;
  double alpha=1,beta=0,*tele,*tele1,*tele2;
  
  tp1=t1;
  tp2=t2;
  if(i1==0){
    transa='T';
    lda=t1.bonddim[i1];
  }
  else if(i1==t1.nbond-1){
    transa='N';
    lda=t1.nelement/t1.bonddim[i1];
  }
  else{
    tp1.shift(i1,t1.nbond-1);
    transa='N';
    lda=t1.nelement/t1.bonddim[i1];
  }
  if(i2==0){
    transb='N';
    ldb=t2.bonddim[i2];
  }
  else if(i2==t2.nbond-1){
    transb='T';
    ldb=t2.nelement/t2.bonddim[i2];
  }
  else{
    tp2.shift(i2,0);
    transb='N';
    ldb=t2.bonddim[i2];
  }
  m=t1.nelement/t1.bonddim[i1];
  n=t2.nelement/t2.bonddim[i2];
  k=t1.bonddim[i1];
  nele=m*n;
  tele=new double[nele];
  dgemm_(&transa,&transb,&m,&n,&k,&alpha,tp1.getptr(),&lda,tp2.getptr(),&ldb,&beta,tele,&m);
  nbond1=t1.nbond+t2.nbond-2;
  bdim=new int[nbond1];
  j=t1.nbond-1;
  for(i=0;i<j;i++)
    bdim[i]=t1.bonddim[(i1+1+i)%t1.nbond];    
  k=t2.nbond-1;
  for(i=0;i<k;i++)
    bdim[j+i]=t2.bonddim[(i2+1+i)%t2.nbond];
  if(telement!=NULL)delete []telement;
  telement=tele;
  if(bonddim!=NULL)delete []bonddim;
  bonddim=bdim;
  nbond=nbond1;
  nelement=nele;
  return *this;
}

//--------------------------------------------------------------------------------------
tensor& tensor::contract_dmrg_overlap_initial(tensor& t1, tensor& t2, int flag){
//--------------------------------------------------------------------------------------
  //flag=0 contract index 0,1 of t1,t2, flag=1 contract index 1,2 of t1,t2
  if(t1.nbond!=3||t2.nbond!=3){
    cout<<"tensor::contract_dmrg_overlap_initial can not perform"<<endl;
    t1.print();
    t2.print();
    exit(0);
  }
  else if(t1.bonddim[0]!=t2.bonddim[0]||t1.bonddim[1]!=t2.bonddim[1]){
    cout<<"tensor::contract_dmrg_overlap_initial the bond dimensions are not consistent"<<endl;
    t1.print();
    t2.print();
    exit(0);
  }
  char transa,transb;
  int nele,nbond1,*bdim;
  int m,n,k,i,j;
  double alpha=1,beta=0,*tele;
  k=t1.bonddim[0]*t1.bonddim[1];
  m=t1.nelement/k;
  n=t2.nelement/k;
  nele=m*n;
  tele=new double[nele];
  transa='T';
  transb='N';
  dgemm_(&transa,&transb,&m,&n,&k,&alpha,t1.getptr(),&k,t2.getptr(),&k,&beta,tele,&m);
  nbond1=2;
  bdim=new int[nbond1];
  bdim[0]=t1.get_bonddim(2);
  bdim[1]=t2.get_bonddim(2);
  if(telement!=NULL)delete []telement;
  telement=tele;
  if(bonddim!=NULL)delete []bonddim;
  bonddim=bdim;
  nbond=nbond1;
  nelement=nele;
  return *this;
}

//--------------------------------------------------------------------------------------
void tensor::contract_dmrg_operator_initial(tensor& t1, tensor& t2, tensor& t3, int flag){
//--------------------------------------------------------------------------------------
  if(flag==0&&(t1.bonddim[1]!=t3.bonddim[0]||t2.bonddim[1]!=t3.bonddim[2])){
    cout<<"tensor::contract_dmrg_operator_initial bonddim not consistent"<<endl;
    exit(0);
  }
  else if(flag==1&&(t1.bonddim[0]!=t3.bonddim[0]||t2.bonddim[0]!=t3.bonddim[2])){
    cout<<"tensor::contract_dmrg_operator_initial bonddim not consistent"<<endl;
    exit(0);
  }
  tensor tmp1,tmp2;
  if(t3.bonddim[0]==1&&t3.bonddim[1]==1&&t3.bonddim[2]==1){
    this->contract_dmrg_overlap_initial(t1,t2,flag);
    (*this)*=t3.getptr()[0];
    this->separateindex(0,bonddim[0],1);
  }
  else if(flag==0){
    tmp1.contract(t1,1,t3,0);
    tmp1.exchangeindex(1,2);
    tmp1.mergeindex(2,3);
    tmp2=t2;
    tmp2.mergeindex(0,1);
    this->contract(tmp1,2,tmp2,0);
  }
  else if(flag==1){
    tmp1.contract(t1,0,t3,0);
    tmp1.shift(3,0);
    tmp1.mergeindex(0,1);
    tmp2=t2;
    tmp2.mergeindex(0,1);
    this->contract(tmp1,0,tmp2,0);
  }
}

//--------------------------------------------------------------------------------------
void tensor::contract_dmrg_operator_transformation(tensor& t1, tensor& t2, tensor& t3, int flag){
//--------------------------------------------------------------------------------------
  if(flag==0&&(t1.bonddim[0]!=t3.bonddim[0]||t2.bonddim[0]!=t3.bonddim[2])){
    cout<<"tensor::contract_dmrg_operator_transformation bonddim not consistent"<<endl;
    exit(0);
  }
  else if(flag==1&&(t1.bonddim[1]!=t3.bonddim[0]||t2.bonddim[1]!=t3.bonddim[2])){
    cout<<"tensor::contract_dmrg_operator_transformation bonddim not consistent"<<endl;
    exit(0);
  }
  tensor tmp1,tmp2;
  char transa,transb;
  int nele,nele1,nbond1,*bdim;
  int m,n,k;
  double alpha=1,beta=0,*tele,*tele1;
  //resume from here tonight
  if(t3.bonddim[1]==1&&(flag==0&&t1.bonddim[1]==1||flag==1&&t1.bonddim[0]==1)){
    k=t3.bonddim[0];
    m=t1.nelement/k;
    n=t3.nelement/k;
    nele1=m*n;
    tele1=new double[nele1];
    transa='T';
    transb='N';
    dgemm_(&transa,&transb,&m,&n,&k,&alpha,t1.getptr(),&k,t3.getptr(),&k,&beta,tele1,&m);
    k=t3.bonddim[2];
    m=nele1/k;
    n=t2.nelement/k;
    nele=m*n;
    tele=new double[nele];
    transa='N';
    transb='N';
    dgemm_(&transa,&transb,&m,&n,&k,&alpha,tele1,&m,t2.getptr(),&k,&beta,tele,&m);
    nbond1=3;
    bdim=new int[nbond1];
    bdim[0]=t1.bonddim[2];
    bdim[1]=t3.bonddim[1];
    bdim[2]=t2.bonddim[2];
    if(telement!=NULL)delete []telement;
    telement=tele;
    if(bonddim!=NULL)delete []bonddim;
    bonddim=bdim;
    nbond=nbond1;
    nelement=nele;
    delete []tele1;
  }
  else if(flag==0){
    tmp1.contract(t1,0,t3,0);
    tmp1.shift(1,0);
    tmp1.mergeindex(2,3);
    tmp2=t2;
    tmp2.mergeindex(0,1);
    this->contract(tmp1,2,tmp2,0);
  }
  else if(flag==1){
    tmp1.contract(t1,1,t3,0);
    tmp1.exchangeindex(1,2);
    tmp1.mergeindex(2,3);
    tmp2=t2;
    tmp2.mergeindex(0,1);
    this->contract(tmp1,2,tmp2,0);
  }
}

//--------------------------------------------------------------------------------------
void tensor::contract_dmrg_operator_pairup(tensor& uu, tensor& vv, tensor& op1, tensor& op2, int flag){
//--------------------------------------------------------------------------------------
  //this routine realize operator pairup in dmrg procedure
  if(flag==0&&(uu.bonddim[0]!=op1.bonddim[0]||vv.bonddim[0]!=op1.bonddim[2])){
    cout<<"tensor::contract_dmrg_operator_pairup bonddim not consistent"<<endl;
    exit(0);
  }
  else if(flag==1&&(uu.bonddim[1]!=op1.bonddim[0]||vv.bonddim[1]!=op1.bonddim[2])){
    cout<<"tensor::contract_dmrg_operator_pairup bonddim not consistent"<<endl;
    exit(0);
  }
  tensor tmp1,tmp2,tmp3;
  if(op1.bonddim[1]==1&&op2.bonddim[0]==1&&op2.bonddim[1]==1&&op2.bonddim[2]==1){
    this->contract_dmrg_operator_transformation(uu,vv,op1,flag);
    this->mergeindex(0,1);
    (*this)*=op2.getptr()[0];
  }
  else if(flag==0){
    tmp1.contract(uu,0,op1,0);
    tmp1.exchangeindex(1,2);
    tmp1.mergeindex(0,1);
    tmp2=op2;
    tmp2.mergeindex(0,1);
    tmp3.contract(tmp1,0,tmp2,0);
    tmp3.mergeindex(1,2);
    tmp2=vv;
    tmp2.mergeindex(0,1);
    this->contract(tmp3,1,tmp2,0);
  }
  else if(flag==1){
    tmp1.contract(uu,1,op1,0);
    tmp1.shift(1,0);
    tmp1.mergeindex(0,1);
    tmp2=op2;
    tmp2.mergeindex(0,1);
    tmp3.contract(tmp2,0,tmp1,0);
    tmp3.mergeindex(0,1);
    tmp2=vv;
    tmp2.mergeindex(0,1);
    this->contract(tmp3,0,tmp2,0);
  }
}

//--------------------------------------------------------------------------------------
tensor& tensor::contract_v2(tensor& t1, int i1, tensor& t2, int i2){
//--------------------------------------------------------------------------------------
  //memory save version, t1, t2 can be destroyed
  if(t1.bonddim[i1]!=t2.bonddim[i2]){
    cout<<"tensor::contract the bond dimension is not consistant"<<endl;
    exit(0);
  }
  char transa,transb;
  int lda,ldb,nele,nbond1,*bdim;
  int m,n,k,i,j,lcontract;
  double alpha=1,beta=0,*tele;

  lcontract=t1.bonddim[i1];
  nbond1=t1.nbond+t2.nbond-2;
  bdim=new int[nbond1];
  j=t1.nbond-1;
  for(i=0;i<j;i++)
    bdim[i]=t1.bonddim[(i1+1+i)%t1.nbond];    
  k=t2.nbond-1;
  for(i=0;i<k;i++)
    bdim[j+i]=t2.bonddim[(i2+1+i)%t2.nbond];

  if(i1==0){
    transa='T';
    lda=lcontract;
  }
  else if(i1==t1.nbond-1){
    transa='N';
    lda=t1.nelement/lcontract;
  }
  else{
    t1.shift(i1,t1.nbond-1);
    transa='N';
    lda=t1.nelement/lcontract;
  }
  if(i2==0){
    transb='N';
    ldb=lcontract;
  }
  else if(i2==t2.nbond-1){
    transb='T';
    ldb=t2.nelement/lcontract;
  }
  else{
    t2.shift(i2,0);
    transb='N';
    ldb=lcontract;
  }
  m=t1.nelement/lcontract;
  n=t2.nelement/lcontract;
  k=lcontract;
  nele=m*n;
  tele=new double[nele];
  dgemm_(&transa,&transb,&m,&n,&k,&alpha,t1.telement,&lda,t2.telement,&ldb,&beta,tele,&m);
  if(telement!=NULL)delete []telement;
  telement=tele;
  if(bonddim!=NULL)delete []bonddim;
  bonddim=bdim;
  nbond=nbond1;
  nelement=nele;
  return *this;
}

//--------------------------------------------------------------------------------------
tensor& tensor::contractindex(int i1, int i2){
//--------------------------------------------------------------------------------------
  int i,j,k,nele,ni1,ni2,nj1,nj2,nj3,m2,m,s0,s1,s2;
  double *tele,tmp;
  if(bonddim[i1]!=bonddim[i2]){
    cout<<"contractindex bonddim not consistant "<<bonddim[i1]<<"\t"<<bonddim[i2]<<endl;
    exit(0);
  }
  if(i1>i2){
    i=i1;
    i1=i2;
    i2=i;
  }
  m=bonddim[i1];
  m2=m*m;
  nele=nelement/m2;
  ni1=1;
  for(i=0;i<i1;i++)ni1*=bonddim[i];
  ni2=ni1;
  for(i=i1+1;i<i2;i++)ni2*=bonddim[i];
  nj1=ni1*m;
  nj2=ni2*m;
  nj3=ni2*m2;
  tele=new double[nele];
  //#pragma omp parallel for default(shared) private(i,j,k,s0,s1,s2,tmp)
  for(i=0;i<nele;i++){
    tele[i]=0;
    s0=i%ni1;
    s1=(i%ni2)/ni1;
    s2=i/ni2;
    tmp=0;
    for(j=0;j<m;j++){
      k=s0+j*ni1+s1*nj1+j*nj2+s2*nj3;
      tmp+=telement[k];
    }
    tele[i]=tmp;
  }
  for(i=i1;i<i2-1;i++)
    bonddim[i]=bonddim[i+1];
  for(i=i2-1;i<nbond-2;i++)
    bonddim[i]=bonddim[i+2];
  nbond=nbond-2;
  if(telement!=NULL)  delete []telement;
  telement=tele;
  nelement=nele;
  return *this;
}

//--------------------------------------------------------------------------------------
void tensor::rescale(double nor){
//--------------------------------------------------------------------------------------
  int i;
  if(fabs(nor)>tolerance)
    //#pragma omp parallel for default(shared) private(i) 
    for(i=0;i<nelement;i++)
      telement[i]/=nor;
}

//--------------------------------------------------------------------------------------
double tensor::rescale(){
//--------------------------------------------------------------------------------------
  int i,id;
  double nor,tmp;
  nor=0;
  for(i=0;i<nelement;i++){
    tmp=telement[i];
    if(fabs(tmp)>fabs(nor))nor=tmp;
  }
  if(fabs(nor)<=tolerance){
    for(i=0;i<nelement;i++)
      telement[i]=0;
    nor=0;
  }
  else{
    for(i=0;i<nelement;i++)
      telement[i]/=nor;
  }
  return nor;
}

//--------------------------------------------------------------------------------------
void obtain_symmetric_matrix_eigenvector(double *mtr, double *w, int mdim){
//--------------------------------------------------------------------------------------
  //this function destroy the symmetric matrix mtr and output its eigenvectors
  char jobz,uplo;
  int n,lda,lwork,i,info;
  double *work,alpha=1,beta=0;
  jobz='V';
  uplo='U';
  n=mdim;
  lda=mdim;
  lwork=-1;
  work=new double[2];
  dsyev_(&jobz,&uplo,&n,mtr,&lda,w,work,&lwork,&info);
  lwork=(int)work[0];
  delete []work;
  work=new double[lwork];
  dsyev_(&jobz,&uplo,&n,mtr,&lda,w,work,&lwork,&info);
  if(info!=0){
    cout<<"dsyev info="<<info<<endl;
    //exit(0);
  }
  delete []work;
  if(info!=0){
    for(i=0;i<mdim;i++){
      cout<<i<<"\t"<<w[i]<<endl;
    }
    exit(0);
  }
}

//--------------------------------------------------------------------------------------
void obtain_symmetric_matrix_eigenvector_2(double *mtr, double *w, int mdim,int dcut){
//--------------------------------------------------------------------------------------
  //this function destroy the symmetric matrix mtr and output its eigenvectors
  int lwork,i,info,il,iu,m,*iwork,*ifail;
  double *work,alpha=1,beta=0,vl,vu,abstol,*z;
  char jobz,uplo,range;
  jobz='V';
  uplo='U';
  range='I';
  abstol=0;
  il=mdim-dcut+1;
  iu=mdim;
  iwork=new int[5*mdim];
  ifail=new int[mdim];
  z=new double[mdim*dcut];
  lwork=-1;
  work=new double[2];
  dsyevx_(&jobz,&range,&uplo,&mdim,mtr,&mdim,&vl,&vu,&il,&iu,&abstol,&m,w,z,&mdim,work,&lwork,iwork,ifail,&info);
  lwork=(int)(work[0])+100*mdim;
  delete []work;
  work=new double[lwork];
  dsyevx_(&jobz,&range,&uplo,&mdim,mtr,&mdim,&vl,&vu,&il,&iu,&abstol,&m,w,z,&mdim,work,&lwork,iwork,ifail,&info);
  if(info!=0){
    cout<<"dsyev info="<<info<<endl;
    //exit(0);
  }
  for(i=0;i<mdim*m;i++)mtr[i]=z[i];
  delete []work;
  delete []iwork;
  delete []ifail;
  delete []z;
  if(info!=0){
    for(i=0;i<dcut;i++){
      cout<<i<<"\t"<<w[i]<<endl;
    }
    exit(0);
  }
}

//--------------------------------------------------------------------------------------
double compare(tensor& t1, tensor& t2){
//--------------------------------------------------------------------------------------
  if(t1.get_nelement()!=t2.get_nelement())return 1;
  int i;
  double diff,nor;
  double tmp;
  diff=0;
  nor=0;
  for(i=0;i<t1.get_nelement();i++){
    tmp=t1.get_telement(i)-t2.get_telement(i);
    diff+=tmp*tmp;
    tmp=t1.get_telement(i)+t2.get_telement(i);
    nor+=tmp*tmp;
  }
  //cout<<"diff="<<diff<<"\tnorm="<<nor<<endl;
  return sqrt(diff/nor);
}

//--------------------------------------------------------------------------------------
void tensor::calculate_difference(tensor& t1, double& dif, double& nor){
//--------------------------------------------------------------------------------------
  int i;
  double tmp;

  for(i=0;i<nelement;i++){
    tmp=telement[i]-t1.get_telement(i);
    dif+=tmp*tmp;
    nor+=telement[i]*telement[i];
  }
}

//--------------------------------------------------------------------------------------
double tensor::inner_prod(tensor& t1){
//--------------------------------------------------------------------------------------
  if(t1.get_nelement()!=nelement){
    cout<<"inner_prod parameter wrong"<<endl;
    exit(0);
  }
  int i;
  double prod;
  prod=0;
  for(i=0;i<nelement;i++)
    prod+=telement[i]*t1.get_telement(i);
  return prod;
}

//--------------------------------------------------------------------------------------
double tensor::take_trace(){
//--------------------------------------------------------------------------------------
  if(nbond!=2||bonddim[0]!=bonddim[1]){
    cout<<"take_trace bonddim wrong"<<endl;
    exit(0);
  }
  int i,m;
  double prod;
  m=bonddim[0];
  prod=0;
  for(i=0;i<m;i++)
    prod+=telement[i+i*m];
  return prod;
}

//--------------------------------------------------------------------------------------
void tensor::direct_sum(int ind, tensor& t1, tensor& t2){
//--------------------------------------------------------------------------------------
  int i,j,m1,m2,nele1,nele2,nele,*aa,*bdim;
  double *tele;
  if(t1.nbond!=t2.nbond){
    cout<<"combine2 not consistent nbond "<<t1.nbond<<"\t"<<t2.nbond<<endl;
    exit(0);
  }
  for(i=0;i<t1.nbond;i++)
    if(i!=ind&&t1.bonddim[i]!=t2.bonddim[i]){
      t1.print();
      t2.print();
      cout<<"combine2 not consistent bonddim "<<t1.bonddim[i]<<"\t"<<t2.bonddim[i]<<endl;
      exit(0);
    }
  m1=t1.get_bonddim(ind);
  m2=t2.get_bonddim(ind);
  nbond=t1.get_nbond();
  aa=new int[nbond];
  bdim=new int[nbond];
  for(i=0;i<nbond;i++)
    bdim[i]=t1.get_bonddim(i);
  bdim[ind]=m1+m2;
  nele1=t1.get_nelement();
  nele2=t2.get_nelement();
  nele=nele1+nele2;
  tele=new double[nele];
  for(i=0;i<nele1;i++){
    get_bond_index(i,nbond,t1.get_bonddim_ptr(),aa);
    get_tensor_index(j,nbond,bdim,aa);
    tele[j]=t1.get_telement(i);
  }
  for(i=0;i<nele2;i++){
    get_bond_index(i,nbond,t2.get_bonddim_ptr(),aa);
    aa[ind]+=m1;
    get_tensor_index(j,nbond,bdim,aa);
    tele[j]=t2.get_telement(i);
  }
  if(telement!=NULL)
    delete []telement;
  telement=tele;
  nelement=nele;
  if(bonddim!=NULL)
    delete []bonddim;
  bonddim=bdim;
  delete []aa;
}

//--------------------------------------------------------------------------------------
void tensor::direct_sum(int ind, tensor& t2){
//--------------------------------------------------------------------------------------
  int i,j,nele2,nele,m1,m2,**aa,*bdim;
  double *tele;
  if(nbond==0){
    *this=t2;
    return;
  }
  if(t2.nbond!=nbond){
    cout<<"combine2 not consistant nbond "<<t2.nbond<<"\t"<<nbond<<endl;
    exit(0);
  }
  for(i=0;i<nbond;i++)
    if(i!=ind&&t2.bonddim[i]!=bonddim[i]){
      cout<<"combine2 not consistant bonddim "<<t2.bonddim[i]<<"\t"<<bonddim[i]<<endl;
      print();
      t2.print();
      exit(0);
    }
  m1=get_bonddim(ind);
  m2=t2.get_bonddim(ind);
  aa=new int*[psize];
  for(i=0;i<psize;i++)
    aa[i]=new int[nbond];
  bdim=new int[nbond];
  for(i=0;i<nbond;i++)
    bdim[i]=bonddim[i];
  bdim[ind]=m1+m2;
  nele2=t2.get_nelement();
  nele=nelement+nele2;
  tele=new double[nele];
  //#pragma omp parallel for default(shared) private(i,j,myrank) schedule(static,1)
  for(i=0;i<nelement;i++){
    //myrank=omp_get_thread_num();
    myrank=0;
    get_bond_index(i,nbond,bonddim,aa[myrank]);
    get_tensor_index(j,nbond,bdim,aa[myrank]);
    tele[j]=telement[i];
  }
  //#pragma omp parallel for default(shared) private(i,j,myrank) schedule(static,1)
  for(i=0;i<nele2;i++){
    //myrank=omp_get_thread_num();
    myrank=0;
    get_bond_index(i,nbond,t2.get_bonddim_ptr(),aa[myrank]);
    aa[myrank][ind]+=m1;
    get_tensor_index(j,nbond,bdim,aa[myrank]);
    tele[j]=t2.get_telement(i);
  }
  delete []telement;
  telement=tele;
  nelement=nele;
  delete []bonddim;
  bonddim=bdim;
  for(i=0;i<psize;i++)
    delete []aa[i];
  delete []aa;
}

//--------------------------------------------------------------------------------------
void tensor::direct_subtract(int ind, int m1, tensor& t1){
//--------------------------------------------------------------------------------------
  int i,j,*bdim1,*bdim2,*aa,nele1,nele2,m2,m,nb;
  double *tele1,*tele2;
  m=bonddim[ind];
  m2=m-m1;
  if(m2==0){
    t1=*this;
    this->clean();
    return;
  }
  if(m1==0)return;
  bdim1=new int[nbond];
  bdim2=new int[nbond];
  aa=new int[nbond];
  for(i=0;i<nbond;i++){
    bdim1[i]=bonddim[i];
    bdim2[i]=bonddim[i];
  }
  bdim1[ind]=m1;
  bdim2[ind]=m2;
  nele1=1;
  nele2=1;
  for(i=0;i<nbond;i++){
    nele1*=bdim1[i];
    nele2*=bdim2[i];
  }
  tele1=new double[nele1];
  tele2=new double[nele2];
  for(i=0;i<nele1;i++){
    get_bond_index(i,nbond,bdim1,aa);
    get_tensor_index(j,nbond,bonddim,aa);
    tele1[i]=telement[j];
  }
  for(i=0;i<nele2;i++){
    get_bond_index(i,nbond,bdim2,aa);
    aa[ind]+=m1;
    get_tensor_index(j,nbond,bonddim,aa);
    tele2[i]=telement[j];
  }
  nb=nbond;
  clean();
  copy(nb,bdim2,tele2);
  t1.copy(nb,bdim1,tele1);
  delete []aa;
  delete []bdim1;
  delete []bdim2;
  delete []tele1;
  delete []tele2;
}

//--------------------------------------------------------------------------------------
void tensor::tensor_product(tensor& t1, tensor& t2){
//--------------------------------------------------------------------------------------
  clean();
  int i,j,k,nele1,nele2,nbond1,nbond2;
  nele1=t1.get_nelement();
  nele2=t2.get_nelement();
  nbond1=t1.get_nbond();
  nbond2=t2.get_nbond();
  if(nele1==0||nele2==0)return;
  nelement=nele1*nele2;
  nbond=nbond1+nbond2;
  bonddim=new int[nbond];
  telement=new double[nelement];
  for(i=0;i<nelement;i++){
    j=i%nele1;
    k=i/nele1;
    telement[i]=t1.get_telement(j)*t2.get_telement(k);
  }
  for(i=0;i<nbond1;i++)
    bonddim[i]=t1.get_bonddim(i);
  for(i=nbond1;i<nbond;i++)
    bonddim[i]=t2.get_bonddim(i-nbond1);
}

//--------------------------------------------------------------------------------------
void get_bond_index(int i0, int nbond, int* bdim, int* aa){
//--------------------------------------------------------------------------------------
  int i;
  for(i=0;i<nbond;i++){
    aa[i]=i0%bdim[i];
    i0/=bdim[i];
  }
}

//--------------------------------------------------------------------------------------
void get_tensor_index(int& i0, int nbond, int* bdim, int* aa){
//--------------------------------------------------------------------------------------
  int i;
  i0=0;
  for(i=nbond-1;i>=0;i--){
    i0*=bdim[i];
    i0+=aa[i];
  }
}

//--------------------------------------------------------------------------------------
void tensor::get_telement(double* tele){
//--------------------------------------------------------------------------------------
  int i;
  for(i=0;i<nelement;i++)
    tele[i]=telement[i];
}

//--------------------------------------------------------------------------------------
void tensor::set_telement(double* tele){
//--------------------------------------------------------------------------------------
  int i;
  for(i=0;i<nelement;i++)
    telement[i]=tele[i];
}

//--------------------------------------------------------------------------------------
bool tensor::check_tensor_null(){
//--------------------------------------------------------------------------------------
  if(nbond==0)return true;
  else return false;
}

//--------------------------------------------------------------------------------------
bool tensor::is_null(){
//--------------------------------------------------------------------------------------
  if(nbond==0) return true;
  else return false;
}

//--------------------------------------------------------------------------------------
bool tensor::is_identity(){
//--------------------------------------------------------------------------------------
  if(nbond!=2)return false;
  if(bonddim[0]!=bonddim[1])return false;
  int i,j,m;
  m=bonddim[0];
  for(i=0;i<m;i++)
    for(j=0;j<m;j++){
      if(i!=j&&fabs(telement[i+j*m])>tolerance){
	cout<<"tensor::is_identity i="<<i<<" j="<<j<<" "<<setprecision(20)<<telement[i+j*m]<<endl;
	return false;
      }
      else if(i==j&&fabs(telement[i+j*m]-1)>tolerance){
	cout<<"tensor::is_identity i="<<i<<" j="<<j<<" "<<setprecision(20)<<telement[i+j*m]<<endl;
	return false;
      }
    }
  for(i=0;i<m*m;i++)telement[i]=0;
  for(i=0;i<m;i++)telement[i+i*m]=1;
  return true;
}

//--------------------------------------------------------------------------------------
bool tensor::is_minus_identity(){
//--------------------------------------------------------------------------------------
  if(nbond!=2)return false;
  if(bonddim[0]!=bonddim[1])return false;
  int i,j,m;
  m=bonddim[0];
  for(i=0;i<m;i++)
    for(j=0;j<m;j++){
      if(i!=j&&fabs(telement[i+j*m])>tolerance){
	cout<<"tensor::is_identity i="<<i<<" j="<<j<<" "<<setprecision(20)<<telement[i+j*m]<<endl;
	return false;
      }
      else if(i==j&&fabs(telement[i+j*m]+1)>tolerance){
	cout<<"tensor::is_identity i="<<i<<" j="<<j<<" "<<setprecision(20)<<telement[i+j*m]<<endl;
	return false;
      }
    }
  for(i=0;i<m*m;i++)telement[i]=0;
  for(i=0;i<m;i++)telement[i+i*m]=-1;
  return true;
}

//--------------------------------------------------------------------------------------
bool tensor::is_zero(){
//--------------------------------------------------------------------------------------
  int i;
  for(i=0;i<nelement;i++)
    if(fabs(telement[i])>tolerance)return false;
  for(i=0;i<nelement;i++)
    telement[i]=0;
  return true;
}

//--------------------------------------------------------------------------------------
bool tensor::is_one(){
//--------------------------------------------------------------------------------------
  int i;
  for(i=0;i<nelement;i++)
    if(fabs(telement[i]-1)>tolerance)return false;
  for(i=0;i<nelement;i++)
    telement[i]=1;
  return true;
}

//--------------------------------------------------------------------------------------
bool tensor::is_minus_one(){
//--------------------------------------------------------------------------------------
  int i;
  for(i=0;i<nelement;i++)
    if(fabs(telement[i]+1)>tolerance)return false;
  for(i=0;i<nelement;i++)
    telement[i]=-1;
  return true;
}

//--------------------------------------------------------------------------------------
bool tensor::is_proportional_to(tensor& t1,double& nor1){
//--------------------------------------------------------------------------------------
  int i,flag;
  double rat;
  if(nbond!=t1.get_nbond())return false;
  for(i=0;i<nbond;i++)
    if(bonddim[i]!=t1.get_bonddim(i))return false;
  if(this->is_zero()){
    for(i=0;i<nelement;i++)
      telement[i]=t1.get_telement(i);
    nor1=0;
    return true;
  }
  rat=0;
  flag=0;
  for(i=0;i<nelement;i++)
    if(fabs(telement[i])<=tolerance&&fabs(t1.get_telement(i))<=tolerance)
      continue;
    else if(fabs(telement[i])<=tolerance&&fabs(t1.get_telement(i))>tolerance)
      return false;
    else if(fabs(telement[i])>tolerance&&fabs(t1.get_telement(i))<=tolerance)
      return false;
    else if(flag==1&&fabs(telement[i]/t1.get_telement(i)-rat)>tolerance)
      return false;
    else if(flag==0){
      rat=telement[i]/t1.get_telement(i);
      flag=1;
    }
  for(i=0;i<nelement;i++)
    telement[i]=t1.get_telement(i);
  nor1=rat;
  return true;
}

//--------------------------------------------------------------------------------------
void tensor::shift_set_identity(int nb,int shift,int* bdim){
//--------------------------------------------------------------------------------------
  int i,bsize;
  clean();
  nbond=nb;
  bonddim=new int[nbond];
  nelement=1;
  for(i=0;i<nbond;i++){
    bonddim[i]=bdim[i];
    nelement*=bonddim[i];
  }
  telement=new double[nelement];
  bsize=nelement/bonddim[nbond-1];
  if(bonddim[nbond-1]>=bsize+shift){
    for(i=0;i<nelement;i++)
      telement[i]=0;
    for(i=0;i<bsize;i++)
      telement[i+(i+shift)*bsize]=1;
  }
  else{
    cout<<"tensor::shift_set_identity, can not do shift_set_identity"<<endl;
    exit(0);
  }
}

//--------------------------------------------------------------------------------------
void tensor::set_telement(int i0, int i1, int i2, double tele){
//--------------------------------------------------------------------------------------
  int i;
  if(nbond<3){
    cout<<"can not set singlet telement, nbond<3"<<endl;
    exit(0);
  }
  i=i0+i1*bonddim[0]+i2*bonddim[0]*bonddim[1];
  telement[i]=tele;
}

//--------------------------------------------------------------------------------------
void tensor::make_identity(int jj){
//--------------------------------------------------------------------------------------
  int i;
  clean();
  nbond=2;
  bonddim=new int[2];
  bonddim[0]=jj+1;
  bonddim[1]=jj+1;
  nelement=bonddim[0]*bonddim[1];
  telement=new double[nelement];
  for(i=0;i<nelement;i++)
    telement[i]=0;
  for(i=0;i<jj+1;i++)
    telement[i+i*(jj+1)]=1;
}

//--------------------------------------------------------------------------------------
void tensor::make_singlet(int jj){
//--------------------------------------------------------------------------------------
  int i,j,m1,m2,sgn;
  clean();
  nbond=2;
  bonddim=new int[nbond];
  bonddim[0]=jj+1;
  bonddim[1]=jj+1;
  nelement=bonddim[0]*bonddim[1];
  telement=new double[nelement];
  for(i=0;i<nelement;i++)
    telement[i]=0;
  for(i=0;i<jj+1;i++){
    m1=2*i-jj;
    for(j=0;j<jj+1;j++){
      m2=2*j-jj;
      if(m1+m2==0){
	if(((jj-m1)/2)%2==0)sgn=1;
	else sgn=-1;
	telement[i+j*(jj+1)]=sgn/sqrt((double)(jj+1));
      }
    }
  }
}

//--------------------------------------------------------------------------------------
void tensor::make_cgc(int j1, int j2, int j3){
//--------------------------------------------------------------------------------------
  //note that j1 j2 j3 are 2 times of the real angular moments
  int i,j,k,i0,m1,m2,m3,k1,k2,k3,k4,k5,k0,kmax;
  double fac1,fac2,fac3,f1,f2,f3,sumk,sgn;
  clean();
  nbond=3;
  bonddim=new int[nbond];
  nelement=(j1+1)*(j2+1)*(j3+1);
  bonddim[0]=j1+1;
  bonddim[1]=j2+1;
  bonddim[2]=j3+1;

  telement=new double[nelement];
  for(i=0;i<nelement;i++)
    telement[i]=0;

  fac1=sqrt((double)(j3+1)*factorial((j3+j1-j2)/2)*factorial((j3-j1+j2)/2)*factorial((j1+j2-j3)/2)/factorial((j1+j2+j3)/2+1));
  //cout<<"fac1="<<fac1<<endl;
  for(i=0;i<j1+1;i++){
    m1=2*i-j1;
    f1=factorial((j1+m1)/2)*factorial((j1-m1)/2);
    for(j=0;j<j2+1;j++){
      m2=2*j-j2;
      f2=factorial((j2+m2)/2)*factorial((j2-m2)/2);
      for(k=0;k<j3+1;k++){
	m3=2*k-j3;
	f3=factorial((j3+m3)/2)*factorial((j3-m3)/2);
	if(m1+m2-m3==0){
	  fac2=sqrt(f1*f2*f3);
	  //cout<<"fac2="<<fac2<<endl;
	  k1=(j1+j2-j3)/2;
	  k2=(j1-m1)/2;
	  k3=(j2+m2)/2;
	  k4=(j3-j2+m1)/2;
	  k5=(j3-j1-m2)/2;
	  kmax=(j1+j2+j3)/2;
	  sumk=0;
	  for(k0=0;k0<=kmax;k0++){
	    if(k1-k0>=0&&k2-k0>=0&&k3-k0>=0&&k4+k0>=0&&k5+k0>=0){
	      if(k0%2==0)sgn=1;
	      else sgn=-1;
	      fac3=sgn/(factorial(k0)*factorial(k1-k0)*factorial(k2-k0)*factorial(k3-k0)*factorial(k4+k0)*factorial(k5+k0));
	      //cout<<"fac3="<<fac3<<endl;
	      sumk+=fac3;
	    }
	  }
	  i0=i+j*bonddim[0]+k*bonddim[0]*bonddim[1];
	  telement[i0]=fac1*fac2*sumk;
	  //cout<<"telement["<<i0<<"]="<<telement[i0]<<endl;
	}
      }
    }
  }
}

//--------------------------------------------------------------------------------------
double factorial(int n){
//--------------------------------------------------------------------------------------
  int i;
  double fac;
  fac=1;
  if(n==0||n==1) return 1;
  for(i=2;i<=n;i++)
    fac*=(double)i;
  return fac;
}

//--------------------------------------------------------------------------------------
void tensor::multiply_cgc(int ind1,int ind2, int j1, int j2, int j3){
//--------------------------------------------------------------------------------------
  if(ind2!=ind1+1){
    cout<<"multiply_cgc indices positions wrong"<<endl;
    exit(0);
  }
  tensor tmp1,tmp2;
  mergeindex(ind1,ind2);
  tmp1.make_cgc(j1,j2,j3);
  tmp1.mergeindex(0,1);
  tmp2.contract(tmp1,0,*this,ind1);
  tmp2.shift(0,ind1);
  *this=tmp2;
}

//--------------------------------------------------------------------------------------
void tensor::multiply_singular_value(int leg, double* w){
//--------------------------------------------------------------------------------------
  int i,*aa;
  if(nbond==0)    return;
  else if(leg>=nbond){
    cout<<"multiply_singular_value wrong leg="<<leg<<"\tnbond="<<nbond<<endl;
    exit(0);
  }
  aa=new int[nbond];
  for(i=0;i<nelement;i++){
    get_bond_index(i,nbond,bonddim,aa);
    telement[i]*=w[aa[leg]];
  }
  delete []aa;
}

//--------------------------------------------------------------------------------------
void tensor::devide_singular_value(int leg, double* w){
//--------------------------------------------------------------------------------------
  int i,*aa;
  if(nbond==0)    return;
  else if(leg>=nbond){
    cout<<"devide_singular_value wrong leg="<<leg<<"\tnbond="<<nbond<<endl;
    exit(0);
  }
  aa=new int[nbond];
  for(i=0;i<nelement;i++){
    get_bond_index(i,nbond,bonddim,aa);
    telement[i]/=w[aa[leg]];
  }
  delete []aa;
}

//--------------------------------------------------------------------------------------
void tensor::svd(tensor &tu, double p1, tensor &tv, double p2, double* diag,int& dcut, int n0){
//--------------------------------------------------------------------------------------
  double *u,*up,*upp,*v,*vp,*vpp,*a,*ap,*app;
  double *w,tol=1.e-24,tol1=1.e-8,num,tmp,alpha=1,beta=0;
  int m,mp,mpp,n,np,npp,i,j,*bdim,max_dc;
  m=1;
  for(i=0;i<n0;i++)m*=bonddim[i];
  n=1;
  for(i=n0;i<nbond;i++)n*=bonddim[i];
  bdim=new int[nbond];
  if(m>n){
    u=new double[m*n];
    v=new double[n*n];
    w=new double[n];
    a=new double[n*n];
    dgemm_("T","N",&n,&n,&m,&alpha,telement,&m,telement,&m,&beta,v,&n);
    memcpy((char*)a,(char*)v,sizeof(double)*n*n);
    obtain_symmetric_matrix_eigenvector(v,w,n);
    np=0;
    for(i=n-1;i>=0;i--)
      if(w[i]/w[n-1]<tol1){
	np=i+1;
	break;
      }
    //if(np>0)
    //cout<<"np="<<np<<endl;
    //for(i=0;i<np;i++)
    //cout<<i<<"\t"<<w[i]/w[n-1]<<endl;
    if(np>0){
      ap=new double[np*n];
      vp=new double[np*np];
      //calculate vt*a*v=vp
      dgemm_("T","N",&np,&n,&n,&alpha,v,&n,a,&n,&beta,ap,&np);
      dgemm_("N","N",&np,&np,&n,&alpha,ap,&np,v,&n,&beta,vp,&np);      
      memcpy((char*)ap,(char*)vp,sizeof(double)*np*np);
      obtain_symmetric_matrix_eigenvector(vp,w,np);
      //for(i=0;i<np;i++)
      //cout<<i<<"\t"<<w[i]/w[n-1]<<endl;
      //calculate v*vp=v
      j=n-np;
      dgemm_("N","N",&np,&np,&np,&alpha,v,&n,vp,&np,&beta,a,&n);
      dgemm_("N","N",&j,&np,&np,&alpha,&(v[np]),&n,vp,&np,&beta,&(a[np]),&n);
      memcpy((char*)v,(char*)a,sizeof(double)*np*n);
      npp=0;
      for(i=np-1;i>=0;i--)
	if(w[i]/w[np-1]<tol1){
	  npp=i+1;
	  break;
	}
      //if(npp>0)
      //cout<<"npp="<<npp<<endl;
      if(npp>0){
	app=new double[npp*np];
	vpp=new double[npp*npp];
	//calcualte vpt*ap*vp=vpp
	dgemm_("T","N",&npp,&np,&np,&alpha,vp,&np,ap,&np,&beta,app,&npp);
	dgemm_("N","N",&npp,&npp,&np,&alpha,app,&npp,vp,&np,&beta,vpp,&npp);
	obtain_symmetric_matrix_eigenvector(vpp,w,npp);
	//for(i=0;i<npp;i++)
	//cout<<i<<"\t"<<w[i]/w[n-1]<<endl;
	//calculate v*vpp
	j=n-npp;
	dgemm_("N","N",&npp,&npp,&npp,&alpha,v,&n,vpp,&npp,&beta,a,&n);
	dgemm_("N","N",&j,&npp,&npp,&alpha,&(v[npp]),&n,vpp,&npp,&beta,&(a[npp]),&n);
	memcpy((char*)v,(char*)a,sizeof(double)*npp*n);
	delete []vpp;
	delete []app;
      }
      delete []vp;
      delete []ap;
    }
    delete []a;
    dgemm_("N","N",&m,&n,&n,&alpha,telement,&m,v,&n,&beta,u,&m);
    num=w[n-1];
    dcut=0;
    for(i=n-1;i>=0;i--){
      if(w[i]/num>tol){
	w[i]=sqrt(w[i]);
	dcut++;
      }
      else w[i]=0;
    }
    /*
    cout<<"###############"<<endl;
    for(i=0;i<dcut;i++)
      cout<<i<<"\t"<<w[i+(n-dcut)]<<endl;
    cout<<"###############"<<endl;
    */
    if(dcut>max_dcut)dcut=max_dcut;
    for(i=n-1;i>=n-dcut;i--)
      if(w[i]>0){
	tmp=pow(w[i],p1-1.);
	for(j=0;j<m;j++)
	  u[j+i*m]*=tmp;
	tmp=pow(w[i],p2);
	for(j=0;j<n;j++)
	  v[j+i*n]*=tmp;
      }
      else{
	for(j=0;j<m;j++)
	  u[j+i*m]=0;
	for(j=0;j<n;j++)
	  v[j+i*n]=0;
      }
  }
  else if(n>=m){
    u=new double[m*m];
    v=new double[m*n];
    w=new double[m];
    a=new double[m*m];
    dgemm_("N","T",&m,&m,&n,&alpha,telement,&m,telement,&m,&beta,u,&m);
    memcpy((char*)a,(char*)u,sizeof(double)*m*m);
    obtain_symmetric_matrix_eigenvector(u,w,m);
    mp=0;
    for(i=m-1;i>=0;i--)
      if(w[i]/w[m-1]<tol1){
	mp=i+1;
	break;
      }
    //if(mp>0)
    //cout<<"mp="<<mp<<endl;
    //for(i=0;i<mp;i++)
    //cout<<i<<"\t"<<w[i]/w[m-1]<<endl;
    if(mp>0){
      ap=new double[mp*m];
      up=new double[mp*mp];
      //calculate ut*a*u=up
      dgemm_("T","N",&mp,&m,&m,&alpha,u,&m,a,&m,&beta,ap,&mp);
      dgemm_("N","N",&mp,&mp,&m,&alpha,ap,&mp,u,&m,&beta,up,&mp);      
      memcpy((char*)ap,(char*)up,sizeof(double)*mp*mp);
      obtain_symmetric_matrix_eigenvector(up,w,mp);
      //for(i=0;i<mp;i++)
      //cout<<i<<"\t"<<w[i]/w[m-1]<<endl;
      //calculate u*up=u
      j=m-mp;
      dgemm_("N","N",&mp,&mp,&mp,&alpha,u,&m,up,&mp,&beta,a,&m);
      dgemm_("N","N",&j,&mp,&mp,&alpha,&(u[mp]),&m,up,&mp,&beta,&(a[mp]),&m);
      memcpy((char*)u,(char*)a,sizeof(double)*mp*m);
      mpp=0;
      for(i=mp-1;i>=0;i--)
	if(w[i]/w[mp-1]<tol1){
	  mpp=i+1;
	  break;
	}
      //if(mpp>0)
      //cout<<"mpp="<<mpp<<endl;
      if(mpp>0){
	app=new double[mpp*mp];
	upp=new double[mpp*mpp];
	//calcualte upt*ap*up=upp
	dgemm_("T","N",&mpp,&mp,&mp,&alpha,up,&mp,ap,&mp,&beta,app,&mpp);
	dgemm_("N","N",&mpp,&mpp,&mp,&alpha,app,&mpp,up,&mp,&beta,upp,&mpp);
	obtain_symmetric_matrix_eigenvector(upp,w,mpp);
	//for(i=0;i<mpp;i++)
	//cout<<i<<"\t"<<w[i]/w[m-1]<<endl;
	//calculate u*upp
	j=m-mpp;
	dgemm_("N","N",&mpp,&mpp,&mpp,&alpha,u,&m,upp,&mpp,&beta,a,&m);
	dgemm_("N","N",&j,&mpp,&mpp,&alpha,&(u[mpp]),&m,upp,&mpp,&beta,&(a[mpp]),&m);
	memcpy((char*)u,(char*)a,sizeof(double)*mpp*m);
	delete []upp;
	delete []app;
      }
      delete []up;
      delete []ap;
    }
    delete []a;
    dgemm_("T","N",&n,&m,&m,&alpha,telement,&m,u,&m,&beta,v,&n);
    num=w[m-1];
    dcut=0;
    for(i=m-1;i>=0;i--){
      if(w[i]/num>tol){
	w[i]=sqrt(w[i]);
	dcut++;
      }
      else w[i]=0;
    }
    /*
    cout<<"###############"<<endl;
    for(i=0;i<dcut;i++)
      cout<<i<<"\t"<<w[i+(m-dcut)]<<endl;
    cout<<"###############"<<endl;
    */
    if(dcut>max_dcut)dcut=max_dcut;
    for(i=m-1;i>=m-dcut;i--)
      if(w[i]>0){
	tmp=pow(w[i],p2-1.);
	for(j=0;j<n;j++)
	  v[j+i*n]*=tmp;
	tmp=pow(w[i],p1);
	for(j=0;j<m;j++)
	  u[j+i*m]*=tmp;
      }
      else{
	for(j=0;j<n;j++)
	  v[j+i*n]=0;
	for(j=0;j<m;j++)
	  u[j+i*m]=0;
      }
  }
  if(m<n)j=m;
  else j=n;
  for(i=0;i<dcut;i++){
    //diag[i]=w[i+(j-dcut)]/w[j-1];
    diag[i]=w[i+(j-dcut)];
  }
  for(i=0;i<n0;i++)bdim[i]=bonddim[i];
  bdim[n0]=dcut;
  tu.copy(n0+1,bdim,&(u[m*(j-dcut)]));
  //tu.print();
  for(i=n0;i<nbond;i++)bdim[i-n0]=bonddim[i];
  bdim[nbond-n0]=dcut;
  tv.copy(nbond-n0+1,bdim,&(v[n*(j-dcut)]));
  //tv.print();
  delete []bdim;
  delete []u;
  delete []v;
  delete []w;
}
