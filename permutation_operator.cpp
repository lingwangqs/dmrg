#include <omp.h>
#include <mpi.h>
#include "dmrg_su2_mpi.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include <time.h>
using namespace std;

void get_tensor_index(int&,int,int*,int*);
void get_bond_index(int,int,int*,int*);
extern tensor *cgc_coef_singlet;
extern int comm_rank;

//--------------------------------------------------------------------------------------
void tensor_su2::make_spinhalf_Qterm(tensor_su2& iden_op1, tensor_su2& plq1, tensor_su2& iden_op2, tensor_su2& plq2){
//--------------------------------------------------------------------------------------
//construct plqterm
  su2bond bb[3];
  int mom[1],dim[1];
  tensor_su2 tmp,tmp1,tmp2;
  mom[0]=2;
  dim[0]=1;
  bb[0].set_su2bond(1,1,mom,dim);

  tmp.make_spinor_start(1);
  plq1.operator_tensor_product_identity(tmp,bb[0]);
  plq1.exchangeindex(3,4);
  plq1.exchangeindex(2,3);
  plq1.fuse(1,2);
  plq1.get_su2bond(2,bb[1]);
  plq1.get_su2bond(3,bb[2]);
  iden_op1.fuse(bb[1],bb[2]);
  plq1.fuse(2,3);
  plq1.make_standard_cgc();
  plq1*=sqrt(3)/2;
  tmp.conjugate(1);
  plq2.operator_tensor_product_identity(tmp,bb[0]);
  plq2.exchangeindex(2,3);
  plq2.fuse(1,2);
  plq2.exchangeindex(1,3);
  plq2.get_su2bond(2,bb[1]);
  plq2.get_su2bond(3,bb[2]);
  iden_op2.fuse(bb[1],bb[2]);
  plq2.fuse(2,3);
  plq2.make_standard_cgc();
  plq2*=-sqrt(3)/2;
}

//--------------------------------------------------------------------------------------
void tensor_su2::make_spinhalf_permutation_rightmove(tensor_su2& iden_op){
//--------------------------------------------------------------------------------------
  //direction -1 out going, direction 1 in going
  //order: down, right, up, left
  su2bond *bb;
  int i,*angm,*bdim,*cdim,bdim2[2];
  double *tele,iden[4];
  tensor tmp1,tmp2,tmp3,tmp;
  bool check;

  clean();
  nbond=4;
  locspin=0;
  bb=new su2bond[nbond];
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  tele=new double[1];
  tele[0]=1;
  bdim[0]=1;
  bdim[1]=1;
  bdim[2]=1;
  bdim[3]=1;
  angm[0]=1;
  bb[0].set_su2bond(1,-1,angm,bdim);
  angm[0]=0;
  angm[1]=2;
  bb[1].set_su2bond(2,-1,angm,bdim);
  angm[0]=1;
  bb[2].set_su2bond(1,1,angm,bdim);
  angm[0]=0;
  angm[1]=2;
  bb[3].set_su2bond(2,1,angm,bdim);
  iden_op.fuse(bb[2],bb[3]);
  
  cgc.set_su2struct(nbond,locspin,bb);
  nten=cgc.get_nten();
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }

  iden[0]=1;
  iden[1]=0;
  iden[2]=0;
  iden[3]=1;
  bdim2[0]=2;
  bdim2[1]=2;
  tmp1.copy(2,bdim2,iden);
  tmp2.make_singlet(1);
  tmp3.direct_product(tmp2,tmp1);
  tmp.direct_product(tmp3,tmp2);
  tmp.exchangeindex(3,5);
  tmp.mergeindex(4,5);
  tmp.mergeindex(1,2);
  tmp*=-2;//here multiply -(2s+1)
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    if(check){
      tmp1.make_cgc(1,1,angm[1]);
      tmp2.make_cgc(1,1,angm[3]);
      tmp1.mergeindex(0,1);
      tmp2.mergeindex(0,1);
      tmp3.contract(tmp,1,tmp1,0);
      tcgc[i].contract(tmp3,1,tmp2,0);
      tarr[i].copy(nbond,bdim,tele);
    }
  }
  this->fuse(2,3);
  this->make_standard_cgc();
  delete []angm;
  delete []bdim;
  delete []cdim;
  delete []tele;
  delete []bb;
}

//--------------------------------------------------------------------------------------
void tensor_su2::make_spinhalf_permutation_leftmove(){
//--------------------------------------------------------------------------------------
  //direction -1 out going, direction 1 in going
  //order: down, right, up, left
  su2bond *bb;
  int i,*angm,*bdim,*cdim,bdim2[2];
  double *tele,iden[4];
  tensor tmp1,tmp2,tmp3,tmp;
  bool check;

  clean();
  nbond=4;
  locspin=0;
  bb=new su2bond[nbond];
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  tele=new double[1];
  tele[0]=1;
  bdim[0]=1;
  bdim[1]=1;
  bdim[2]=1;
  bdim[3]=1;
  angm[0]=1;
  bb[0].set_su2bond(1,-1,angm,bdim);
  angm[0]=0;
  angm[1]=2;
  bb[1].set_su2bond(2,-1,angm,bdim);
  angm[0]=1;
  bb[2].set_su2bond(1,1,angm,bdim);
  angm[0]=0;
  angm[1]=2;
  bb[3].set_su2bond(2,1,angm,bdim);

  cgc.set_su2struct(nbond,locspin,bb);
  nten=cgc.get_nten();
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }

  iden[0]=1;
  iden[1]=0;
  iden[2]=0;
  iden[3]=1;
  bdim2[0]=2;
  bdim2[1]=2;
  tmp1.copy(2,bdim2,iden);
  tmp3.direct_product(tmp1,tmp1);
  tmp.direct_product(tmp3,tmp1);
  tmp.exchangeindex(1,2);
  tmp.mergeindex(3,4);
  tmp.mergeindex(1,2);
  tmp.exchangeindex(0,1);
  tmp.shift(3,0);
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    if(check){
      tmp1.make_cgc(1,1,angm[1]);
      tmp2.make_cgc(1,1,angm[3]);
      tmp1.mergeindex(0,1);
      tmp2.mergeindex(0,1);
      tmp3.contract(tmp,1,tmp1,0);
      tcgc[i].contract(tmp3,1,tmp2,0);
      tarr[i].copy(nbond,bdim,tele);
    }
  }
  this->fuse(2,3);
  this->make_standard_cgc();
  delete []angm;
  delete []bdim;
  delete []cdim;
  delete []tele;
  delete []bb;
}

//--------------------------------------------------------------------------------------
void tensor_su2::make_spinhalf_permutation_rightmove_4legs(tensor_su2& iden_op){
//--------------------------------------------------------------------------------------
  //direction -1 out going, direction 1 in going
  //order: down, right, up, left
  su2bond *bb;
  int i,*angm,*bdim,*cdim,bdim2[2];
  double *tele,iden[4];
  tensor tmp1,tmp2,tmp3,tmp;
  bool check;

  clean();
  nbond=4;
  locspin=0;
  bb=new su2bond[nbond];
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  tele=new double[1];
  tele[0]=1;
  bdim[0]=1;
  bdim[1]=1;
  bdim[2]=1;
  bdim[3]=1;
  angm[0]=1;
  bb[0].set_su2bond(1,-1,angm,bdim);
  angm[0]=0;
  angm[1]=2;
  bb[1].set_su2bond(2,-1,angm,bdim);
  angm[0]=1;
  bb[2].set_su2bond(1,1,angm,bdim);
  angm[0]=0;
  angm[1]=2;
  bb[3].set_su2bond(2,1,angm,bdim);
  iden_op.fuse(bb[2],bb[3]);
  
  cgc.set_su2struct(nbond,locspin,bb);
  nten=cgc.get_nten();
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }

  iden[0]=1;
  iden[1]=0;
  iden[2]=0;
  iden[3]=1;
  bdim2[0]=2;
  bdim2[1]=2;
  tmp1.copy(2,bdim2,iden);
  tmp2.make_singlet(1);
  tmp3.direct_product(tmp2,tmp1);
  tmp.direct_product(tmp3,tmp2);
  tmp.exchangeindex(3,5);
  tmp.mergeindex(4,5);
  tmp.mergeindex(1,2);
  tmp*=-2;//here multiply -(2s+1)
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    if(check){
      tmp1.make_cgc(1,1,angm[1]);
      tmp2.make_cgc(1,1,angm[3]);
      tmp1.mergeindex(0,1);
      tmp2.mergeindex(0,1);
      tmp3.contract(tmp,1,tmp1,0);
      tcgc[i].contract(tmp3,1,tmp2,0);
      tarr[i].copy(nbond,bdim,tele);
    }
  }
  //this->fuse(2,3);
  //this->make_standard_cgc();
  delete []angm;
  delete []bdim;
  delete []cdim;
  delete []tele;
  delete []bb;
}

//--------------------------------------------------------------------------------------
void tensor_su2::make_spinhalf_permutation_leftmove_4legs(){
//--------------------------------------------------------------------------------------
  //direction -1 out going, direction 1 in going
  //order: down, right, up, left
  su2bond *bb;
  int i,*angm,*bdim,*cdim,bdim2[2];
  double *tele,iden[4];
  tensor tmp1,tmp2,tmp3,tmp;
  bool check;

  clean();
  nbond=4;
  locspin=0;
  bb=new su2bond[nbond];
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  tele=new double[1];
  tele[0]=1;
  bdim[0]=1;
  bdim[1]=1;
  bdim[2]=1;
  bdim[3]=1;
  angm[0]=1;
  bb[0].set_su2bond(1,-1,angm,bdim);
  angm[0]=0;
  angm[1]=2;
  bb[1].set_su2bond(2,-1,angm,bdim);
  angm[0]=1;
  bb[2].set_su2bond(1,1,angm,bdim);
  angm[0]=0;
  angm[1]=2;
  bb[3].set_su2bond(2,1,angm,bdim);

  cgc.set_su2struct(nbond,locspin,bb);
  nten=cgc.get_nten();
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }

  iden[0]=1;
  iden[1]=0;
  iden[2]=0;
  iden[3]=1;
  bdim2[0]=2;
  bdim2[1]=2;
  tmp1.copy(2,bdim2,iden);
  tmp3.direct_product(tmp1,tmp1);
  tmp.direct_product(tmp3,tmp1);
  tmp.exchangeindex(1,2);
  tmp.mergeindex(3,4);
  tmp.mergeindex(1,2);
  tmp.exchangeindex(0,1);
  tmp.shift(3,0);
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    if(check){
      tmp1.make_cgc(1,1,angm[1]);
      tmp2.make_cgc(1,1,angm[3]);
      tmp1.mergeindex(0,1);
      tmp2.mergeindex(0,1);
      tmp3.contract(tmp,1,tmp1,0);
      tcgc[i].contract(tmp3,1,tmp2,0);
      tarr[i].copy(nbond,bdim,tele);
    }
  }
  //this->fuse(2,3);
  //this->make_standard_cgc();
  delete []angm;
  delete []bdim;
  delete []cdim;
  delete []tele;
  delete []bb;
}

//--------------------------------------------------------------------------------------
void tensor_su2::make_spinhalf_permutation_leftend(){
//--------------------------------------------------------------------------------------
  //direction -1 out going, direction 1 in going
  //order: down, right, up, left
  su2bond *bb;
  int i,*angm,*bdim,*cdim,bdim2[2];
  double *tele,iden[4];
  tensor tmp1,tmp2,tmp3,tmp;
  bool check;

  clean();
  nbond=3;
  locspin=0;
  bb=new su2bond[nbond];
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  tele=new double[1];
  tele[0]=1;
  bdim[0]=1;
  bdim[1]=1;
  bdim[2]=1;
  bdim[3]=1;
  
  angm[0]=1;
  bb[0].set_su2bond(1,-1,angm,bdim);
  angm[0]=0;
  angm[1]=2;
  bb[1].set_su2bond(2,-1,angm,bdim);
  angm[0]=1;
  bb[2].set_su2bond(1,1,angm,bdim);

  cgc.set_su2struct(nbond,locspin,bb);
  nten=cgc.get_nten();
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }

  iden[0]=1;
  iden[1]=0;
  iden[2]=0;
  iden[3]=1;
  bdim2[0]=2;
  bdim2[1]=2;
  tmp1.copy(2,bdim2,iden);
  tmp2.make_singlet(1);
  tmp.direct_product(tmp2,tmp1);
  tmp.mergeindex(1,2);
  tmp*=sqrt(2);
  
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    if(check){
      tmp1.make_cgc(1,1,angm[1]);
      tmp1.mergeindex(0,1);
      tcgc[i].contract(tmp,1,tmp1,0);
      tcgc[i].shift(1,0);
      tarr[i].copy(nbond,bdim,tele);
    }
  }
  this->make_standard_cgc();
  delete []angm;
  delete []bdim;
  delete []cdim;
  delete []tele;
  delete []bb;
}

//--------------------------------------------------------------------------------------
void tensor_su2::make_spinhalf_permutation_rightend(){
//--------------------------------------------------------------------------------------
  //direction -1 out going, direction 1 in going
  //order: down, right, up, left
  su2bond *bb;
  int i,*angm,*bdim,*cdim,bdim2[2];
  double *tele,iden[4];
  tensor tmp1,tmp2,tmp3,tmp;
  bool check;

  clean();
  nbond=3;
  locspin=0;
  bb=new su2bond[nbond];
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  tele=new double[1];
  tele[0]=1;
  bdim[0]=1;
  bdim[1]=1;
  bdim[2]=1;
  bdim[3]=1;
  
  angm[0]=0;
  angm[1]=2;
  bb[0].set_su2bond(2,1,angm,bdim);
  angm[0]=1;
  bb[1].set_su2bond(1,1,angm,bdim);
  angm[0]=1;
  bb[2].set_su2bond(1,-1,angm,bdim);

  cgc.set_su2struct(nbond,locspin,bb);
  nten=cgc.get_nten();
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }

  iden[0]=1;
  iden[1]=0;
  iden[2]=0;
  iden[3]=1;
  bdim2[0]=2;
  bdim2[1]=2;
  tmp1.copy(2,bdim2,iden);
  tmp2.make_singlet(1);
  tmp.direct_product(tmp2,tmp1);
  tmp.exchangeindex(1,2);
  tmp.mergeindex(0,1);
  tmp*=-sqrt(2);
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    if(check){
      tmp1.make_cgc(1,1,angm[0]);
      tmp1.mergeindex(0,1);
      tcgc[i].contract(tmp1,0,tmp,0);
      tarr[i].copy(nbond,bdim,tele);
    }
  }
  this->make_standard_cgc();
  this->conjugate(0);
  this->shift(2,0);
  this->make_standard_cgc();
  delete []angm;
  delete []bdim;
  delete []cdim;
  delete []tele;
  delete []bb;
}

//--------------------------------------------------------------------------------------
void tensor_su2::make_spinhalf_permutation_rightend_v2(){
//--------------------------------------------------------------------------------------
  //direction -1 out going, direction 1 in going
  //order: down, right, up, left
  su2bond *bb;
  int i,*angm,*bdim,*cdim,bdim2[2];
  double *tele,iden[4];
  tensor tmp1,tmp2,tmp3,tmp;
  bool check;

  clean();
  nbond=3;
  locspin=0;
  bb=new su2bond[nbond];
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  tele=new double[1];
  tele[0]=1;
  bdim[0]=1;
  bdim[1]=1;
  bdim[2]=1;
  bdim[3]=1;
  
  angm[0]=0;
  angm[1]=2;
  bb[0].set_su2bond(2,1,angm,bdim);
  angm[0]=1;
  bb[1].set_su2bond(1,1,angm,bdim);
  angm[0]=1;
  bb[2].set_su2bond(1,-1,angm,bdim);

  cgc.set_su2struct(nbond,locspin,bb);
  nten=cgc.get_nten();
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }

  iden[0]=1;
  iden[1]=0;
  iden[2]=0;
  iden[3]=1;
  bdim2[0]=2;
  bdim2[1]=2;
  tmp1.copy(2,bdim2,iden);
  tmp2.make_singlet(1);
  tmp.direct_product(tmp2,tmp1);
  tmp.exchangeindex(1,2);
  tmp.mergeindex(0,1);
  tmp*=-sqrt(2);
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    if(check){
      tmp1.make_cgc(1,1,angm[0]);
      tmp1.mergeindex(0,1);
      tcgc[i].contract(tmp1,0,tmp,0);
      tarr[i].copy(nbond,bdim,tele);
    }
  }
  this->make_standard_cgc();
  /*
  this->conjugate(0);
  this->shift(2,0);
  this->make_standard_cgc();
  */
  delete []angm;
  delete []bdim;
  delete []cdim;
  delete []tele;
  delete []bb;
}

//--------------------------------------------------------------------------------------
void tensor_su2::make_spinor_end(int physpn){
//--------------------------------------------------------------------------------------
  su2bond *bb;
  int *angm,*bdim;
  double *tele;
  //direction -1 out going, direction 1 in going
  //order: down, horizontal, up
  clean();
  nbond=3;
  locspin=0;
  bb=new su2bond[nbond];
  angm=new int[nbond];
  bdim=new int[nbond];
  tele=new double[1];
  angm[0]=2;
  bdim[0]=1;
  bb[0].set_su2bond(1,1,angm,bdim);
  angm[0]=physpn;
  bdim[0]=1;
  bb[1].set_su2bond(1,1,angm,bdim);
  angm[0]=physpn;
  bdim[0]=1;
  bb[2].set_su2bond(1,-1,angm,bdim);

  cgc.set_su2struct(nbond,locspin,bb);
  nten=cgc.get_nten();
  tarr=new tensor[1];
  tcgc=new tensor[1];
  parr=new tensor*[1];
  pcgc=new tensor*[1];
  parr[0]=&(tarr[0]);
  pcgc[0]=&(tcgc[0]);

  bdim[0]=1;
  bdim[1]=1;
  bdim[2]=1;
  tele[0]=1;
  tarr[0].copy(nbond,bdim,tele);
  tcgc[0].make_cgc(2,physpn,physpn);
  delete []angm;
  delete []bdim;
  delete []tele;
  delete []bb;
}

//--------------------------------------------------------------------------------------
void tensor_su2::make_spinhalf_on_vec_left(tensor_su2& vec, tensor_su2& sigma, int flag){
//--------------------------------------------------------------------------------------
  tensor_su2 tmp,stmp;
  tmp=vec;
  tmp.exchangeindex(0,1);
  if(flag==0){
    this->contract(tmp,0,sigma,0);
    this->exchangeindex(1,3);
    this->fuse(2,3);
  }
  else if(flag==1){
    this->contract(sigma,2,tmp,0);
    this->exchangeindex(1,2);
    this->fuse(0,1);
  }
  this->make_standard_cgc();
}

//--------------------------------------------------------------------------------------
void tensor_su2::make_spinhalf_on_vec_right(tensor_su2& vec, tensor_su2& sigma, int flag){
//--------------------------------------------------------------------------------------
  tensor_su2 tmp,stmp;
  if(flag==0){
    this->contract(vec,0,sigma,0);
    this->exchangeindex(1,3);
    this->exchangeindex(0,1);
    this->fuse(2,3);
  }
  else if(flag==1){
    this->contract(sigma,2,vec,0);
    this->exchangeindex(1,2);
    this->fuse(0,1);
    this->exchangeindex(0,1);
  }
  this->make_standard_cgc();
}

//--------------------------------------------------------------------------------------
void tensor::direct_product(tensor& t1,tensor& t2){
//--------------------------------------------------------------------------------------
  int i,j,k,nb1,nb2,*bdim1,*bdim2,*bb;
  clean();
  nb1=t1.get_nbond();
  nb2=t2.get_nbond();
  nbond=nb1+nb2;
  bonddim=new int[nbond];
  bb=new int[nbond];
  bdim1=new int[nb1];
  bdim2=new int[nb2];
  for(i=0;i<nb1;i++){
    bonddim[i]=t1.get_bonddim(i);
    bdim1[i]=bonddim[i];
  }
  for(i=0;i<nb2;i++){
    bonddim[i+nb1]=t2.get_bonddim(i);
    bdim2[i]=bonddim[i+nb1];
  }
  nelement=1;
  for(i=0;i<nbond;i++)
    nelement*=bonddim[i];
  telement=new double[nelement];
  for(i=0;i<nelement;i++){
    get_bond_index(i,nbond,bonddim,bb);
    get_tensor_index(j,nb1,bdim1,bb);
    get_tensor_index(k,nb2,bdim2,&(bb[nb1]));
    telement[i]=t1.get_telement(j)*t2.get_telement(k);
  }
  delete []bdim1;
  delete []bdim2;
  delete []bb;
}

//--------------------------------------------------------------------------------------
tensor& tensor::contract_dmrg_permutation(tensor& uu, tensor& vv, tensor& vec, tensor& op1, tensor& op2, int flag){
//--------------------------------------------------------------------------------------
  tensor tmp1,tmp2,tmp3,tmp4;
  int bd1,bd2,bd3,bd4,bd5,bd6,zero=0;
  double fac1,fac2;
  bd1=op1.get_bonddim(0);
  bd2=op1.get_bonddim(1);
  bd3=op1.get_bonddim(2);
  bd4=op2.get_bonddim(0);
  bd5=op2.get_bonddim(1);
  bd6=op2.get_bonddim(2);
  if(bd1==1&&bd2==1&&bd3==1&&bd4==1&&bd5==1&&bd6==1||flag==0&&bd1==1&&bd2==1&&bd3==2&&bd4==1&&bd5==1&&bd6==2||flag==1&&bd1==2&&bd2==1&&bd3==1&&bd4==2&&bd5==1&&bd6==1){
    tmp1=uu;
    tmp2=vec;
    tmp3=vv;
    tmp1.mergeindex(0,1);
    tmp2.mergeindex(0,1);
    tmp3.mergeindex(0,1);
    tmp4.contract(tmp1,0,tmp2,0);
    this->contract(tmp4,1,tmp3,0);
    this->separateindex(0,this->get_bonddim(0),1);
    fac1=op1.inner_prod(op2);
    (*this)*=fac1;
    return *this;
  }
  else if(flag==0){
    tmp1=vv;
    tmp1.exchangeindex(0,1);
    tmp2.contract(op1,0,tmp1,0);
    tmp2.exchangeindex(1,2);
    tmp2.mergeindex(0,1);
    tmp1=vec;
    tmp1.mergeindex(1,2);
    tmp3.contract(tmp1,1,tmp2,0);
    tmp3.mergeindex(0,1);
    tmp1.contract(uu,1,op2,0);
    tmp1.exchangeindex(1,2);
    tmp1.mergeindex(2,3);
    this->contract(tmp1,2,tmp3,0);
  }
  else if(flag==1){
    tmp1.contract(uu,0,op2,1);
    tmp1.exchangeindex(1,2);
    tmp1.mergeindex(0,1);
    tmp2=vec;
    tmp2.mergeindex(0,1);
    tmp3.contract(tmp1,0,tmp2,0);
    tmp3.mergeindex(1,2);
    tmp1=op1;
    tmp1.exchangeindex(0,1);
    tmp2.contract(tmp1,0,vv,0);
    tmp2.exchangeindex(1,2);
    tmp2.mergeindex(0,1);
    this->contract(tmp3,1,tmp2,0);
  }
  return *this;
  //this->print();
}


//--------------------------------------------------------------------------------------
void tensor_su2::permutation(tensor_su2& uu, tensor_su2& vv, tensor_su2& vec, tensor_su2& op1, tensor_su2& op2, 
                             const LookupTable_5& fac_permutation_left, const LookupTable_5& fac_permutation_rght,int flag){
//--------------------------------------------------------------------------------------
  int a0,a1,a2,b0,b1,b2,c0,c1,c2,d0,d1,d2,e0,e1,e2,nmoma0,nmoma1,nmoma2,nmomb0,nmomb1,nmomb2,nmomc0,nmomc1,nmomc2,nmomd0,nmomd1,nmomd2,nmome0,nmome1,nmome2,i0,i1,i2,j0,j1,j2,k0,k1,k2,l0,l1,l2,s0,s1,s2,m,n,p,q,r,t,i;
  su2bond *bb;
  tensor tmp1,tmp2,tmp3;
  int bdim[5],angm[3];
  double fac;
  clean();
  nbond=3;
  locspin=0;
  bb=new su2bond[nbond];
  uu.get_su2bond(2,bb[0]);
  if(flag==0)op2.get_su2bond(1,bb[1]);
  else if(flag==1)op1.get_su2bond(2,bb[1]);
  vv.get_su2bond(2,bb[2]);
  bb[2].invert_bonddir();
  cgc.set_su2struct(nbond,locspin,bb);
  delete []bb;
  nten=cgc.get_nten();
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }
  nmoma0=uu.get_nmoment(0);
  nmoma1=uu.get_nmoment(1);
  nmoma2=uu.get_nmoment(2);
  nmomb0=vv.get_nmoment(0);
  nmomb1=vv.get_nmoment(1);
  nmomb2=vv.get_nmoment(2);
  nmome0=vec.get_nmoment(0);
  nmome1=vec.get_nmoment(1);
  nmome2=vec.get_nmoment(2);
  nmomc0=op1.get_nmoment(0);
  nmomc1=op1.get_nmoment(1);
  nmomc2=op1.get_nmoment(2);
  nmomd0=op2.get_nmoment(0);
  nmomd1=op2.get_nmoment(1);
  nmomd2=op2.get_nmoment(2);

  for(i0=0;i0<nmoma0;i0++){
    a0=uu.get_angularmoment(0,i0);
    for(i1=0;i1<nmoma1;i1++){
      a1=uu.get_angularmoment(1,i1);
      for(i2=0;i2<nmoma2;i2++){
	a2=uu.get_angularmoment(2,i2);
	angm[0]=a0;
	angm[1]=a1;
	angm[2]=a2;
	if(uu.check_angularmoments(angm)==false)continue;

	for(j0=0;j0<nmomb0;j0++){
	  b0=vv.get_angularmoment(0,j0);
	  for(j1=0;j1<nmomb1;j1++){
	    b1=vv.get_angularmoment(1,j1);
	    for(j2=0;j2<nmomb2;j2++){
	      b2=vv.get_angularmoment(2,j2);
	      angm[0]=b0;
	      angm[1]=b1;
	      angm[2]=b2;
	      if(vv.check_angularmoments(angm)==false)continue;

	      for(k0=0;k0<nmomc0;k0++){
		c0=op1.get_angularmoment(0,k0);
		for(k1=0;k1<nmomc1;k1++){
		  c1=op1.get_angularmoment(1,k1);
		  for(k2=0;k2<nmomc2;k2++){
		    c2=op1.get_angularmoment(2,k2);
		    angm[0]=c0;
		    angm[1]=c1;
		    angm[2]=c2;
		    if(op1.check_angularmoments(angm)==false)continue;

		    for(l0=0;l0<nmomd0;l0++){
		      d0=op2.get_angularmoment(0,l0);
		      for(l1=0;l1<nmomd1;l1++){
			d1=op2.get_angularmoment(1,l1);
			for(l2=0;l2<nmomd2;l2++){
			  d2=op2.get_angularmoment(2,l2);
			  angm[0]=d0;
			  angm[1]=d1;
			  angm[2]=d2;
			  if(op2.check_angularmoments(angm)==false)continue;

			  for(s0=0;s0<nmome0;s0++){
			    e0=vec.get_angularmoment(0,s0);
			    for(s1=0;s1<nmome1;s1++){
			      e1=vec.get_angularmoment(1,s1);
			      for(s2=0;s2<nmome2;s2++){
				e2=vec.get_angularmoment(2,s2);
				angm[0]=e0;
				angm[1]=e1;
				angm[2]=e2;
				if(vec.check_angularmoments(angm)==false)continue;
			  
				if(flag==0&&(a1!=d0||b1!=c0||c2!=d2||c1!=e1||a0!=e0||b0!=e2))continue;
				if(flag==1&&(a0!=d1||b0!=c1||c0!=d0||d2!=e1||a1!=e0||b1!=e2))continue;

				m=i0+i1*nmoma0+i2*nmoma0*nmoma1;
				n=j0+j1*nmomb0+j2*nmomb0*nmomb1;
				p=k0+k1*nmomc0+k2*nmomc0*nmomc1;
				q=l0+l1*nmomd0+l2*nmomd0*nmomd1;
				r=s0+s1*nmome0+s2*nmome0*nmome1;
				if(uu.is_null(m)||vv.is_null(n)||vec.is_null(r)||op1.is_null(p)||op2.is_null(q))continue;
				tmp1.contract_dmrg_permutation(*(uu.get_parr(m)),*(vv.get_parr(n)),*(vec.get_parr(r)),*(op1.get_parr(p)),*(op2.get_parr(q)),flag);
				if(flag==0){
				  i=i2+l1*nmoma2+j2*nmoma2*nmomd1;
				  t=c1/2+((c2-1)/2)*3+(d1/2)*6;
				}
				else if(flag==1){
				  i=i2+k2*nmoma2+j2*nmoma2*nmomc2;
				  t=c2/2+((c0-1)/2)*3+(d2/2)*6;
				}
				if(flag==0)
				  fac=fac_permutation_left(a0,a2,b0,b2,t);
				else if(flag==1)
				  fac=fac_permutation_rght(a1,a2,b1,b2,t);
				tmp1*=fac;
				if(tarr[i].is_null()){
				  tarr[i]=tmp1;
				  if(flag==0)tcgc[i].make_cgc(a2,d1,b2);
				  else if(flag==1)tcgc[i].make_cgc(a2,c2,b2);
				}
				else
				  tarr[i]+=tmp1;
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

//--------------------------------------------------------------------------------------
void tensor::contract_dmrg_operator_initial(tensor& t1, tensor& t2, tensor& t3, tensor& endt, int flag){
//--------------------------------------------------------------------------------------
  if(flag==0&&(t1.bonddim[1]!=t3.bonddim[0]||t2.bonddim[1]!=t3.bonddim[2]||t1.bonddim[0]!=endt.bonddim[0]||t2.bonddim[0]!=endt.bonddim[1])){
    cout<<"tensor::contract_dmrg_operator_initial bonddim not consistent"<<endl;
    exit(0);
  }
  else if(flag==1&&(t1.bonddim[0]!=t3.bonddim[0]||t2.bonddim[0]!=t3.bonddim[2]||t1.bonddim[1]!=endt.bonddim[0]||t2.bonddim[1]!=endt.bonddim[1])){
    cout<<"tensor::contract_dmrg_operator_initial bonddim not consistent"<<endl;
    exit(0);
  }
  tensor tmp1,tmp2;
  if(flag==0)
    tmp1.contract(endt,0,t1,0);
  else if(flag==1){
    tmp1.contract(endt,0,t1,1);
    tmp1.shift(2,0);
  }
  this->contract_dmrg_operator_initial(tmp1,t2,t3,flag);
}

//--------------------------------------------------------------------------------------
void tensor_su2::operator_initial(tensor_su2& uu, tensor_su2& vv, tensor_su2& op, tensor_su2& endt, 
                                  const LookupTable_4& fac_operator_onsite_left, const LookupTable_4& fac_operator_onsite_rght,int flag){//end tensor is not identity
//--------------------------------------------------------------------------------------
  //build three index operator from current site
  int a0,a1,a2,b0,b1,b2,c0,c1,c2,d0,d1,nmoma0,nmoma1,nmoma2,nmomb0,nmomb1,nmomb2,nmomc0,nmomc1,nmomc2,nmomd0,nmomd1,i0,i1,i2,j0,j1,j2,k0,k1,k2,l0,l1,m,n,p,q,i;
  su2bond *bb;
  tensor tmp1;
  int angm[3];
  double fac;
  clean();
  if(op.get_nbond()==0)return;
  nbond=3;
  locspin=0;
  bb=new su2bond[nbond];
  uu.get_su2bond(2,bb[0]);
  op.get_su2bond(1,bb[1]);
  vv.get_su2bond(2,bb[2]);
  bb[2].invert_bonddir();
  cgc.set_su2struct(nbond,locspin,bb);
  delete []bb;
  nten=cgc.get_nten();
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }
  nmoma0=uu.get_nmoment(0);
  nmoma1=uu.get_nmoment(1);
  nmoma2=uu.get_nmoment(2);
  nmomb0=vv.get_nmoment(0);
  nmomb1=vv.get_nmoment(1);
  nmomb2=vv.get_nmoment(2);
  nmomc0=op.get_nmoment(0);
  nmomc1=op.get_nmoment(1);
  nmomc2=op.get_nmoment(2);
  nmomd0=endt.get_nmoment(0);
  nmomd1=endt.get_nmoment(1);

  for(i0=0;i0<nmoma0;i0++){
    a0=uu.get_angularmoment(0,i0);
    for(i1=0;i1<nmoma1;i1++){
      a1=uu.get_angularmoment(1,i1);
      for(i2=0;i2<nmoma2;i2++){
	a2=uu.get_angularmoment(2,i2);
	angm[0]=a0;
	angm[1]=a1;
	angm[2]=a2;
	if(uu.check_angularmoments(angm)==false)continue;

	for(j0=0;j0<nmomb0;j0++){
	  b0=vv.get_angularmoment(0,j0);
	  for(j1=0;j1<nmomb1;j1++){
	    b1=vv.get_angularmoment(1,j1);
	    for(j2=0;j2<nmomb2;j2++){
	      b2=vv.get_angularmoment(2,j2);
	      angm[0]=b0;
	      angm[1]=b1;
	      angm[2]=b2;
	      if(vv.check_angularmoments(angm)==false)continue;

	      for(k0=0;k0<nmomc0;k0++){
		c0=op.get_angularmoment(0,k0);
		for(k1=0;k1<nmomc1;k1++){
		  c1=op.get_angularmoment(1,k1);
		  for(k2=0;k2<nmomc2;k2++){
		    c2=op.get_angularmoment(2,k2);
		    angm[0]=c0;
		    angm[1]=c1;
		    angm[2]=c2;
		    if(op.check_angularmoments(angm)==false)continue;

		    for(l0=0;l0<nmomd0;l0++){
		      d0=endt.get_angularmoment(0,l0);
		      for(l1=0;l1<nmomd1;l1++){
			d1=endt.get_angularmoment(1,l1);
			  angm[0]=d0;
			  angm[1]=d1;
			  if(endt.check_angularmoments(angm)==false)continue;

			  if(flag==0&&(a0!=d0||b0!=d1||c0!=a1||c2!=b1)||flag==1&&(a1!=d0||b1!=d1||c0!=a0||c2!=b0))continue;
			  m=i0+i1*nmoma0+i2*nmoma0*nmoma1;
			  n=j0+j1*nmomb0+j2*nmomb0*nmomb1;
			  p=k0+k1*nmomc0+k2*nmomc0*nmomc1;
			  i=i2+k1*nmoma2+j2*nmoma2*nmomc1;
			  q=l0+l1*nmomd0;
			  if(uu.is_null(m)||vv.is_null(n)||op.is_null(p)||endt.is_null(q))continue;
			  tmp1.contract_dmrg_operator_initial(*(uu.get_parr(m)),*(vv.get_parr(n)),*(op.get_parr(p)),*(endt.get_parr(q)),flag);
			  if(flag==0)
			    fac=fac_operator_onsite_left(a2,c1,b2,a0);
			  else if(flag==1)
			    fac=fac_operator_onsite_rght(a2,c1,b2,a1);
			  tmp1*=fac;
			  if(tarr[i].is_null()){
			    tarr[i]=tmp1;
			    tcgc[i].make_cgc(a2,c1,b2);
			  }
			  else
			    tarr[i]+=tmp1;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

//------------------------------------------------------------------------------
void dmrg_su2::prepare_translationysquare_operators(){
//------------------------------------------------------------------------------
  int i;
  su2bond bb[4];
  tensor_su2 tmp1,tmp,ringbak[8];

  ring[2].get_su2bond(1,bb[0]);
  tmp.operator_tensor_product_identity(ring[2],bb[0]);
  tmp.exchangeindex(0,4);
  tmp.exchangeindex(2,3);
  tmp.fuse(0,1,2,3);
  tmp.shift(2,0);
  tmp.make_standard_cgc();
  ringbak[2]=tmp;
  for(i=0;i<2;i++){
    tmp.operator_tensor_product_identity(ring[i],bb[0]);
    tmp.exchangeindex(0,4);
    tmp.exchangeindex(0,1);
    tmp.fuse(0,1,2,3);
    tmp.shift(2,0);
    tmp.make_standard_cgc();
    ringbak[i]=tmp;
  }
  for(i=0;i<2;i++){
    tmp.contract(ring[2],0,ring[i],0);
    tmp.exchangeindex(1,3);
    tmp.exchangeindex(2,3);
    tmp.fuse(0,1,2,3);
    tmp.make_standard_cgc();
    tmp1.contract(ringbak[i],2,tmp,1);
    ringbak[i]=tmp1;
  }

  tmp.operator_tensor_product_identity(ring[5],bb[0]);
  tmp.shift(4,0);
  tmp.exchangeindex(2,4);
  tmp.fuse(0,1,2,3);
  tmp.exchangeindex(1,2);
  tmp.make_standard_cgc();
  ringbak[5]=tmp;
  for(i=0;i<2;i++){
    tmp.operator_tensor_product_identity(ring[3+i],bb[0]);
    tmp.exchangeindex(1,4);
    tmp.fuse(0,1,2,3);
    tmp.exchangeindex(1,2);
    tmp.make_standard_cgc();
    ringbak[3+i]=tmp;
  }
  for(i=0;i<2;i++){
    tmp.contract(ring[5],1,ring[3+i],1);
    tmp.shift(3,0);
    tmp.exchangeindex(0,1);
    tmp.fuse(0,1,2,3);
    tmp.make_standard_cgc();
    tmp1.contract(tmp,1,ringbak[3+i],0);
    ringbak[3+i]=tmp1;
  }
  tmp.contract(ring[6],2,ring[6],0);
  tmp.fuse(1,2);
  tmp.make_standard_cgc();
  ringbak[6]=tmp;
  tmp.contract(ring[7],2,ring[7],0);
  tmp.fuse(1,2);
  tmp.make_standard_cgc();
  ringbak[7]=tmp;
  for(i=0;i<8;i++){
    ring[i]=ringbak[i];
    ring[i].print();
  }
}
