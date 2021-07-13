#include "tensor_su2.hpp"
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

using namespace std;
extern "C"{
  void dsort2_(int*,double*,int*);
  //int memcmp(char*,char*,int);
  //void memcpy(char*,char*,int);
}
extern int comm_rank,psize,myrank;
extern tensor ***spin_op,*cgc_coef_singlet,*identity;
extern double **spin_op_trace,****fac_operator_onsite_left,****fac_operator_onsite_rght,*****fac_operator_transformation_left,*****fac_operator_transformation_rght,****fac_operator_pairup_left,****fac_operator_pairup_rght,***fac_hamilt_vec;
extern "C"{
  void dgemm_(char*,char*,int*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);
}
void obtain_symmetric_matrix_eigenvector(double*,double*,int);
bool sum_direct_product(tensor&,tensor&,tensor&,tensor&);
void sum_direct_product(tensor&,tensor&,tensor&,tensor&,tensor&,tensor&);

extern int max_dcut;
//--------------------------------------------------------------------------------------
tensor_su2::tensor_su2(){
//--------------------------------------------------------------------------------------
  nten=0;
  nbond=0;
  locspin=0;
  tarr=NULL;
  parr=NULL;
  tcgc=NULL;
  pcgc=NULL;
}

//--------------------------------------------------------------------------------------
tensor_su2::~tensor_su2(){
//--------------------------------------------------------------------------------------
  clean();
}

//--------------------------------------------------------------------------------------
void tensor_su2::clean(){
//--------------------------------------------------------------------------------------
  if(tarr!=NULL)delete []tarr;
  if(tcgc!=NULL)delete []tcgc;
  if(parr!=NULL)delete []parr;
  if(pcgc!=NULL)delete []pcgc;
  cgc.clean();
  tarr=NULL;
  tcgc=NULL;
  parr=NULL;
  pcgc=NULL;
  nbond=0;
  nten=0;
  locspin=0;
  nten=0;
}

//--------------------------------------------------------------------------------------
void tensor_su2::print(){
//--------------------------------------------------------------------------------------
  cgc.print();
  int i,j,*angm,*bdim,*cdim;
  bool check;
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    //if(!parr[i]->is_zero()&&!pcgc[i]->is_zero()){
    if(check&&!parr[i]->is_null()&&!pcgc[i]->is_null()){
      cout<<"****i="<<i<<"****"<<endl;
      for(j=0;j<nbond;j++)
	cout<<"angm["<<j<<"]="<<angm[j]<<endl;
      parr[i]->print();
      pcgc[i]->print();
    }
  }
  delete []angm;
  delete []bdim;
  delete []cdim;
}

//--------------------------------------------------------------------------------------
su2struct& tensor_su2::get_cgc(){
//--------------------------------------------------------------------------------------
  return cgc;
}

//--------------------------------------------------------------------------------------
tensor_su2& tensor_su2::operator = (double a){
//--------------------------------------------------------------------------------------
  int i;
  for(i=0;i<nten;i++)
    if(!parr[i]->is_null()){
      *(parr[i])=a;
    }
}

//--------------------------------------------------------------------------------------
tensor_su2& tensor_su2::operator = (tensor_su2& t1){
//--------------------------------------------------------------------------------------
  int i;
  clean();
  cgc=t1.cgc;
  nbond=t1.get_nbond();
  nten=t1.get_nten();
  locspin=t1.get_locspin();
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
    if(!t1.is_null(i)){
      tarr[i]=*(t1.get_parr(i));
      tcgc[i]=*(t1.get_pcgc(i));
    }
  }
  return *this;
}

//--------------------------------------------------------------------------------------
bool tensor_su2::operator != (tensor_su2& t1){
//--------------------------------------------------------------------------------------
  int i;
  if(cgc!=t1.cgc) return true;
  if(nbond!=t1.get_nbond()) return true;
  if(nten!=t1.get_nten()) return true;
  if(locspin!=t1.get_locspin()) return true;
  for(i=0;i<nten;i++){
    if(!parr[i]->is_null()||!t1.is_null(i)){
      if(*(parr[i])!=*(t1.get_parr(i))){
	parr[i]->print();
	t1.get_parr(i)->print();
	return true;
      }
      if(*(pcgc[i])!=*(t1.get_pcgc(i))){
	return true;
      }
    }
  }
  return false;
}

//--------------------------------------------------------------------------------------
bool tensor_su2::operator == (tensor_su2& t1){
//--------------------------------------------------------------------------------------
  if(*this!=t1)return false;
  else return true;
}

//--------------------------------------------------------------------------------------
tensor_su2& tensor_su2::operator += (tensor_su2& t1){
//--------------------------------------------------------------------------------------
  int i,*angm,*bdim,*cdim;
  bool check;
  double nor1,nor2;
  if(cgc!=t1.cgc){
    cout<<"tensor_su2::operator += can not operate"<<endl;
    t1.print();
    print();
    exit(0);    
  }
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  for(i=0;i<nten;i++){
    check=get_tensor_argument(i,angm,bdim,cdim);
    if(!parr[i]->is_null()||!t1.get_parr(i)->is_null()){
      if(t1.get_parr(i)->is_null()||t1.get_parr(i)->is_zero()||t1.get_pcgc(i)->is_zero())
	continue;
      else if(parr[i]->is_null()||parr[i]->is_zero()||pcgc[i]->is_zero()){
	*(parr[i])=*(t1.get_parr(i));
	*(pcgc[i])=*(t1.get_pcgc(i));
      }
      else if(!sum_direct_product(*(parr[i]),*(pcgc[i]),*(t1.get_parr(i)),*(t1.get_pcgc(i)))){
	cout<<"tensor_su2::operator += can not perform"<<endl;
	this->print();
	t1.print();
	exit(0);
      }
    }
  }
  delete []angm;
  delete []bdim;
  delete []cdim;
  return *this;
}

//--------------------------------------------------------------------------------------
tensor_su2& tensor_su2::operator -= (tensor_su2& t1){
//--------------------------------------------------------------------------------------
  int i,*angm,*bdim,*cdim;
  bool check;
  tensor tmp;
  double nor1,nor2;
  if(cgc!=t1.cgc){
    cout<<"tensor_su2::operator += can not operate"<<endl;
    t1.print();
    print();
    exit(0);    
  }
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  for(i=0;i<nten;i++){
    check=get_tensor_argument(i,angm,bdim,cdim);
    if(!parr[i]->is_null()||!t1.get_parr(i)->is_null()){
      if(t1.get_parr(i)->is_null()||t1.get_parr(i)->is_zero()||t1.get_pcgc(i)->is_zero())
	continue;
      else if(parr[i]->is_null()||parr[i]->is_zero()||pcgc[i]->is_zero()){
	*(parr[i])=*(t1.get_parr(i));
	*(pcgc[i])=*(t1.get_pcgc(i));
	*(parr[i])*=-1;
      }
      else{
	tmp=*(t1.get_parr(i));
	tmp*=-1;
	if(!sum_direct_product(*(parr[i]),*(pcgc[i]),tmp,*(t1.get_pcgc(i)))){
	  cout<<"tensor_su2::operator -= can not perform"<<endl;
	  this->print();
	  t1.print();
	  exit(0);
	}
      }
    }
  }
  delete []angm;
  delete []bdim;
  delete []cdim;
  return *this;
}

//--------------------------------------------------------------------------------------
tensor_su2& tensor_su2::operator *= (double fac){
//--------------------------------------------------------------------------------------
  int i;
  for(i=0;i<nten;i++)
    if(!parr[i]->is_null())
      *(parr[i])*=fac;
  return *this;
}

//--------------------------------------------------------------------------------------
tensor_su2& tensor_su2::operator /= (double fac){
//--------------------------------------------------------------------------------------
  int i;
  for(i=0;i<nten;i++)
    if(!parr[i]->is_null())
      *(parr[i])/=fac;
  return *this;
}

//--------------------------------------------------------------------------------------
void tensor_su2::cgc_make_scalar_operator(){
//--------------------------------------------------------------------------------------
  if(nbond!=2||locspin!=0){
    cout<<"tensor_su2::cgc_check_scalar_operator, not a scalar operator"<<endl;
    cout<<"nbond="<<nbond<<"\tlocspin="<<locspin<<endl;
    exit(0);
  }
  int i;
  double fac;
  for(i=0;i<nten;i++){
    if(!parr[i]->is_null()){
      if(pcgc[i]->is_zero())
	*(parr[i])=0;	
      else{ 
	fac=pcgc[i]->rescale();
	if(pcgc[i]->is_identity())
	  *(parr[i])*=fac;
	else if(pcgc[i]->is_minus_identity())
	  *(parr[i])*=-fac;
	else if(pcgc[i]->get_bonddim(0)!=pcgc[i]->get_bonddim(1)){
	  //else if(pcgc[i]->get_bonddim(0)!=pcgc[i]->get_bonddim(1)&&parr[i]->is_zero()){
	  pcgc[i]->clean();
	  parr[i]->clean();
	}
	else{
	  cout<<"tensor_su2::cgc_make_scalar_operator, cgc is not identity"<<endl;
	  pcgc[i]->print();
	  parr[i]->print();
	  exit(0);
	}
      }
    }
  }
}

//--------------------------------------------------------------------------------------
void tensor_su2::cgc_make_cgc(){
//--------------------------------------------------------------------------------------
  int i,*angm,*bdim,*cdim;
  bool check;
  if(nbond!=3){
    cout<<"tensor_su2::cgc_make_cgc, can not perform make_cgc"<<endl;
    exit(0);
  }
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  for(i=0;i<nten;i++){
    if(!parr[i]->is_null()){
      check=cgc.get_tensor_argument(i,angm,bdim,cdim);
      pcgc[i]->make_cgc(angm[0],angm[1],angm[2]);
    }
  }
  delete []angm;
  delete []bdim;
  delete []cdim;
  cgc.set_bonddir(0,1);
  cgc.set_bonddir(1,1);
  cgc.set_bonddir(2,-1);
}

//--------------------------------------------------------------------------------------
double tensor_su2::normalize_vector(){
//--------------------------------------------------------------------------------------
  int i,j,*angm,*bdim,*cdim,*flag;
  bool check;
  double nor,nor1,nor2;
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  flag=new int[nten];
  nor=0;
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    if(check){
      flag[i]=1;
      nor1=parr[i]->inner_prod(*(parr[i]));
      nor2=pcgc[i]->inner_prod(*(pcgc[i]));      
      nor+=nor1*nor2;
    }
    else flag[i]=0;
  }

  nor=sqrt(nor);
  for(i=0;i<nten;i++)
    if(flag[i])
      *(parr[i])/=nor;
  if(comm_rank==0)cout<<"nor="<<nor<<endl;
  delete []angm;
  delete []bdim;
  delete []cdim;
  delete []flag;
  return nor;
}

//--------------------------------------------------------------------------------------
void tensor_su2::multiply_singular_value(int leg,double *w){
//--------------------------------------------------------------------------------------
  int i,j,nmom,*angm,*bdim,*cdim;
  double **ww;
  bool check;
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  nmom=cgc.get_nmoment(leg);
  ww=new double*[nmom];
  j=0;
  for(i=0;i<nmom;i++){
    ww[i]=&(w[j]);
    j+=cgc.get_bonddim(leg,i);
  }
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    if(check){
      j=cgc.get_angularmoment_index(leg,angm[leg]);
      parr[i]->multiply_singular_value(leg,ww[j]);
    }
  }
  delete []ww;
  delete []angm;
  delete []bdim;
  delete []cdim;
}

//--------------------------------------------------------------------------------------
void tensor_su2::devide_singular_value(int leg,double *w){
//--------------------------------------------------------------------------------------
  int i,j,nmom,*angm,*bdim,*cdim;
  double **ww;
  bool check;
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  nmom=cgc.get_nmoment(leg);
  ww=new double*[nmom];
  j=0;
  for(i=0;i<nmom;i++){
    ww[i]=&(w[j]);
    j+=cgc.get_bonddim(leg,i);
  }
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    if(check){
      j=cgc.get_angularmoment_index(leg,angm[leg]);
      parr[i]->devide_singular_value(leg,ww[j]);
    }
  }
  delete []ww;
  delete []angm;
  delete []bdim;
  delete []cdim;
}

//--------------------------------------------------------------------------------------
void tensor_su2::diagonalize(){
//--------------------------------------------------------------------------------------
  int m,n,a,i,j,nmom1,nmom2;
  double *eig,*aa;  
  if(nbond==2){
    nmom1=get_nmoment(0);
    nmom2=get_nmoment(1);
    if(nmom1!=nmom2){
      cout<<"tensor_su2::diagonalize can not perform"<<endl;
      cout<<"nmom1!=nmom2"<<endl;
      exit(0);
    }
    for(i=0;i<nmom1;i++){
      m=cgc.get_bonddim(0,i);
      n=cgc.get_bonddim(1,i);
      a=cgc.get_angularmoment(0,i);
      cout<<"ang="<<a<<endl;
      //parr[i+i*nmom1]->print();
      //if(a!=0)continue;
      if(m!=n){
	cout<<"tensor_su2::diagonalize can not perform"<<endl;
	cout<<"m!=n"<<endl;
	exit(0);
      }
      eig=new double[m];
      aa=new double[m*m];
      parr[i+i*nmom1]->get_telement(aa);
      obtain_symmetric_matrix_eigenvector(aa,eig,m);
      for(j=0;j<m;j++)
	if(j<10)
	  cout<<"eig["<<j<<"]=\t"<<setw(10)<<eig[j]<<endl;
      delete []aa;
      delete []eig;
    }
  }
  else{
    cout<<"tensor_su2::diagonalize can not perform"<<endl;
    cout<<"nbond!=2"<<endl;
    exit(0);
  }
}

//--------------------------------------------------------------------------------------
int tensor_su2::get_nmoment(int i){
//--------------------------------------------------------------------------------------
  return cgc.get_nmoment(i);
}
//--------------------------------------------------------------------------------------
int tensor_su2::get_bonddir(int i){
//--------------------------------------------------------------------------------------
  return cgc.get_bonddir(i);
}
//--------------------------------------------------------------------------------------
int tensor_su2::get_angularmoment(int i,int j){
//--------------------------------------------------------------------------------------
  return cgc.get_angularmoment(i,j);
}
//--------------------------------------------------------------------------------------
int tensor_su2::get_bonddim(int i,int j){
//--------------------------------------------------------------------------------------
  return cgc.get_bonddim(i,j);
}
//--------------------------------------------------------------------------------------
int tensor_su2::get_cgcdim(int i,int j){
//--------------------------------------------------------------------------------------
  return cgc.get_cgcdim(i,j);
}

//--------------------------------------------------------------------------------------
void tensor_su2::get_su2bond(int i, su2bond& b){
//--------------------------------------------------------------------------------------
  cgc.get_su2bond(i,b);
}

//--------------------------------------------------------------------------------------
int tensor_su2::get_angularmoment_index(int i0, int angm){
//--------------------------------------------------------------------------------------
  return  cgc.get_angularmoment_index(i0,angm);
}

//--------------------------------------------------------------------------------------
int tensor_su2::get_tensor_index(int *angm){
//--------------------------------------------------------------------------------------
  return cgc.get_tensor_index(angm);
}

//--------------------------------------------------------------------------------------
bool tensor_su2::get_tensor_argument(int i0, int* angm, int* bdim, int* cdim){
//--------------------------------------------------------------------------------------
  return cgc.get_tensor_argument(i0,angm,bdim,cdim);
}

//--------------------------------------------------------------------------------------
tensor* tensor_su2::get_parr(int i){
//--------------------------------------------------------------------------------------
  if(i<nten)
    return parr[i];
  else{
    cout<<"tensor_su2::get_parr tensor index is too large"<<endl;
    exit(0);
  }
}

//--------------------------------------------------------------------------------------
tensor* tensor_su2::get_pcgc(int i){
//--------------------------------------------------------------------------------------
  if(i<nten)
    return pcgc[i];
  else{
    cout<<"tensor_su2::get_pcgc tensor index is too large"<<endl;
    exit(0);
  }
}

//--------------------------------------------------------------------------------------
void tensor_su2::take_conjugate(){
//--------------------------------------------------------------------------------------
  cgc.take_conjugate();
}

//--------------------------------------------------------------------------------------
void tensor_su2::get_nelement(int &nele1,int &nele2){
//--------------------------------------------------------------------------------------
  int i;

  nele1=0;
  nele2=0;
  for(i=0;i<nten;i++)
    if(!parr[i]->is_null()){
      nele1+=parr[i]->get_nelement();
      nele2+=pcgc[i]->get_nelement();
    }
}

//--------------------------------------------------------------------------------------
void tensor_su2::get_telement(double *tele1, double* tele2){
//--------------------------------------------------------------------------------------
  int i,nele1,nele2;

  nele1=0;
  nele2=0;
  for(i=0;i<nten;i++)
    if(!parr[i]->is_null()){
      parr[i]->get_telement(&(tele1[nele1]));
      pcgc[i]->get_telement(&(tele2[nele2]));
      nele1+=parr[i]->get_nelement();
      nele2+=pcgc[i]->get_nelement();
    }
}

//--------------------------------------------------------------------------------------
void tensor_su2::get_telement(double *tele1, double* tele2, int* iflag){
//--------------------------------------------------------------------------------------
  int i,nele1,nele2;

  nele1=0;
  nele2=0;
  for(i=0;i<nten;i++){
    iflag[i]=0;
    if(!parr[i]->is_null()){
      iflag[i]=1;
      parr[i]->get_telement(&(tele1[nele1]));
      pcgc[i]->get_telement(&(tele2[nele2]));
      nele1+=parr[i]->get_nelement();
      nele2+=pcgc[i]->get_nelement();
    }
  }
}

//--------------------------------------------------------------------------------------
void tensor_su2::take_conjugate(int ind){
//--------------------------------------------------------------------------------------
  if(ind>=nbond){
    cout<<"tensor_su2::take_conjugate on bond, bond index wrong\t"<<ind<<endl;
    exit(0);
  }
  int i,j,nmom,*angm,*bdim,*cdim;
  tensor *conj,tmp;
  bool check;

  nmom=cgc.get_nmoment(ind);
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  for(i=0;i<nten;i++){
    check=get_tensor_argument(i,angm,bdim,cdim);
    if(check&&!pcgc[i]->is_null()){
      tmp.contract(*(pcgc[i]),ind,cgc_coef_singlet[angm[ind]],0);
      tmp.shift(0,(ind+1)%nbond);
      //tmp*=sqrt(angm[ind]+1);
      *(pcgc[i])=tmp;
    }
  }
  cgc.take_conjugate(ind);
  delete []angm;
  delete []bdim;
  delete []cdim;
}

//--------------------------------------------------------------------------------------
void tensor_su2::conjugate(int ind){
//--------------------------------------------------------------------------------------
  if(ind>=nbond){
    cout<<"tensor_su2::take_conjugate on bond, bond index wrong\t"<<ind<<endl;
    exit(0);
  }
  int i,j,nmom,*angm,*bdim,*cdim;
  tensor *conj,tmp;
  bool check;

  nmom=cgc.get_nmoment(ind);
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  for(i=0;i<nten;i++){
    check=get_tensor_argument(i,angm,bdim,cdim);
    if(check&&!pcgc[i]->is_null()){
      tmp.contract(*(pcgc[i]),ind,cgc_coef_singlet[angm[ind]],0);
      tmp.shift(0,(ind+1)%nbond);
      tmp*=sqrt(angm[ind]+1);
      *(pcgc[i])=tmp;
    }
  }
  cgc.take_conjugate(ind);
  delete []angm;
  delete []bdim;
  delete []cdim;
}

//--------------------------------------------------------------------------------------
bool tensor_su2::is_null(int i){
//--------------------------------------------------------------------------------------
  return parr[i]->is_null();
}

//--------------------------------------------------------------------------------------
bool tensor_su2::is_null(){
//--------------------------------------------------------------------------------------
  if(nbond==0)return true;
  else return false;
}

//--------------------------------------------------------------------------------------
bool tensor_su2::check_angularmoments(int* angm){
//--------------------------------------------------------------------------------------
  return cgc.check_angularmoments(angm);
}

//--------------------------------------------------------------------------------------
double tensor_su2::inner_prod(tensor_su2& t1){
//--------------------------------------------------------------------------------------
  int i;
  double prod,prod1,prod2;
  if(cgc!=t1.cgc){
    cout<<"tensor_su2::inner_prod t1 has different su2 structure"<<endl;
    exit(0);
  }
  prod=0;
  for(i=0;i<nten;i++){
    if(!parr[i]->is_null()){
      prod1=parr[i]->inner_prod(*(t1.get_parr(i)));
      prod2=pcgc[i]->inner_prod(*(t1.get_pcgc(i)));
      if(fabs(prod2-1.)>1.e-12){
	//cout<<"inner_prod wrong\t"<<prod2<<endl;
	//exit(0);
      }
      prod+=prod1;
      //prod+=prod1*prod2;
    }
  }
  return prod;
}

//--------------------------------------------------------------------------------------
double tensor_su2::ss_inner_prod(tensor_su2& t1){
//--------------------------------------------------------------------------------------
  int i,*angm,*bdim,*cdim,a0,a1,a2;
  double prod,prod1,prod2;
  bool check;
  if(cgc!=t1.cgc){
    cout<<"tensor_su2::inner_prod t1 has different su2 structure"<<endl;
    exit(0);
  }
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  prod=0;
  for(i=0;i<nten;i++){
    check=get_tensor_argument(i,angm,bdim,cdim);
    if(check&&!parr[i]->is_null()){
      a0=angm[0];
      a1=angm[1];
      a2=angm[2];
      prod1=parr[i]->inner_prod(*(t1.get_parr(i)));
      prod2=fac_hamilt_vec[a0][a1][a2];
      prod+=prod1*prod2;
    }
  }
  delete []angm;
  delete []bdim;
  delete []cdim;
  return prod;
}

//--------------------------------------------------------------------------------------
double tensor_su2::take_trace(){
//--------------------------------------------------------------------------------------
  int i,*angm,*bdim,*cdim;
  double prod,prod1,prod2;
  bool check;
  if(nbond!=2){
    cout<<"tensor_su2::take_trace  wrong su2 tensor"<<endl;
    exit(0);
  }
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  prod=0;
  for(i=0;i<nten;i++){
    check=get_tensor_argument(i,angm,bdim,cdim);
    if(check&&angm[0]==angm[1]&&!parr[i]->is_null()){
      prod1=parr[i]->take_trace();
      //prod2=pcgc[i]->take_trace();
      prod2=1;
      prod+=prod1*prod2;
    }
  }
  delete []angm;
  delete []bdim;
  delete []cdim;
  return prod;
}

//--------------------------------------------------------------------------------------
void tensor_su2::left2right_vectran(){
//--------------------------------------------------------------------------------------
  //this function change uu projector from j1,j2 fuse to j3, to j2,j3 fuse to j1
  int i,j,ishift,*angm,*bdim,*cdim,*angm1;
  bool check;
  tensor **parr1,**pcgc1,tmp1,tmp2;
  su2struct cgc1;
  double nor;
  //three index, shift(1,0)
  ishift=2;
  cgc1=cgc;
  cgc1.shift(1,0);
  parr1=new tensor*[nten];
  pcgc1=new tensor*[nten];
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  angm1=new int[nbond];
  for(i=0;i<nten;i++){
    check=get_tensor_argument(i,angm,bdim,cdim);
    for(j=0;j<nbond;j++)
      angm1[(j+ishift)%nbond]=angm[j];
    j=cgc1.get_tensor_index(angm1);
    parr1[j]=parr[i];
    pcgc1[j]=pcgc[i];
    if(!parr1[j]->is_null()){
      parr1[j]->shift(1,0);
      pcgc1[j]->make_cgc(angm1[0],angm1[1],angm1[2]);
    }
  }
  cgc=cgc1;
  delete []parr;
  delete []pcgc;
  parr=parr1;
  pcgc=pcgc1;
  delete []angm;
  delete []angm1;
  delete []bdim;
  delete []cdim;
  cgc.set_bonddir(0,1);
  cgc.set_bonddir(1,1);
  cgc.set_bonddir(2,-1);
}

//--------------------------------------------------------------------------------------
void tensor_su2::right2left_vectran(){
//--------------------------------------------------------------------------------------
  //this function change uu projector from j1,j2 fuse to j3, to j2,j3 fuse to j1
  int i,j,ishift,*angm,*bdim,*cdim,*angm1;
  bool check;
  tensor **parr1,**pcgc1,tmp1,tmp2;
  su2struct cgc1;
  double nor;
  //three index, shift(2,0)
  ishift=1;
  cgc1=cgc;
  cgc1.shift(2,0);
  parr1=new tensor*[nten];
  pcgc1=new tensor*[nten];
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  angm1=new int[nbond];
  for(i=0;i<nten;i++){
    check=get_tensor_argument(i,angm,bdim,cdim);
    for(j=0;j<nbond;j++)
      angm1[(j+ishift)%nbond]=angm[j];
    j=cgc1.get_tensor_index(angm1);
    parr1[j]=parr[i];
    pcgc1[j]=pcgc[i];
    if(!parr1[j]->is_null()){
      parr1[j]->shift(2,0);
      pcgc1[j]->make_cgc(angm1[0],angm1[1],angm1[2]);
    }
  }
  cgc=cgc1;
  delete []parr;
  delete []pcgc;
  parr=parr1;
  pcgc=pcgc1;
  delete []angm;
  delete []angm1;
  delete []bdim;
  delete []cdim;
  cgc.set_bonddir(0,1);
  cgc.set_bonddir(1,1);
  cgc.set_bonddir(2,-1);
}

//--------------------------------------------------------------------------------------
void tensor_su2::reflection_vectran(){
//--------------------------------------------------------------------------------------
  int i,j,*angm,*bdim,*cdim;
  bool check,check2;
  tensor tmp1,tmp2;
  double nor;
  this->exchangeindex(0,1);
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  for(i=0;i<nten;i++){
    check=get_tensor_argument(i,angm,bdim,cdim);
    if(check&&!parr[i]->is_null()){
      tmp1.make_cgc(angm[0],angm[1],angm[2]);
      check2=pcgc[i]->is_proportional_to(tmp1,nor);
      if(check2){
	(*parr[i])*=nor;
	(*pcgc[i])=tmp1;
      }
      else{
	cout<<"reflection_vectran couldn't perform"<<endl;
      }
    }
  }
  delete []angm;
  delete []bdim;
  delete []cdim;
}

//--------------------------------------------------------------------------------------
void tensor_su2::fuse(su2bond& b1, su2bond& b2){
//--------------------------------------------------------------------------------------
  int i,j,k,l,n1,n2,n3,m1,m2,m3,ms,me,d1,d2,d3;
  int *bdim,*pos;
  su2bond *bb;
  clean();
  nbond=3;
  if(b1.get_bonddir()!=b2.get_bonddir()){
    cout<<"tensor_su2::fuse two bonds had different directions"<<endl;
    exit(0);
  }
  bb=new su2bond[nbond];
  bb[0]=b1;
  bb[1]=b2;
  bb[2].fuse(bb[0],bb[1]);

  locspin=0;
  cgc.set_su2struct(3,0,bb);
  nten=cgc.get_nten();

  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];

  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }

  n1=bb[0].get_nmoment();
  n2=bb[1].get_nmoment();
  n3=bb[2].get_nmoment();

  pos=new int[n3];
  bdim=new int[nbond];

  for(i=0;i<n3;i++)
    pos[i]=0;

  for(j=0;j<n2;j++){
    m2=bb[1].get_angularmoment(j);
    d2=bb[1].get_bonddim(j);
    bdim[1]=d2;
    for(i=0;i<n1;i++){
      m1=bb[0].get_angularmoment(i);
      d1=bb[0].get_bonddim(i);
      bdim[0]=d1;
      ms=abs(m1-m2);
      me=m1+m2;
      for(m3=ms;m3<=me;m3+=2){
	k=bb[2].get_angularmoment_index(m3);
	d3=bb[2].get_bonddim(k);
	bdim[2]=d3;
	l=i+j*n1+k*n1*n2;
	pcgc[l]->make_cgc(m1,m2,m3);
	parr[l]->shift_set_identity(3,pos[k],bdim);
	pos[k]+=d1*d2;
      }
    }
  }
  for(i=0;i<n3;i++){
    if(bb[2].get_bonddim(i)!=pos[i]){
      cout<<"tensor_su2::fuse, wrong in fuse two su2 bond"<<endl;
      exit(0);
    }
  }
  delete []pos;
  delete []bdim;
  delete []bb;
}

//--------------------------------------------------------------------------------------
void tensor_su2::fuse(int ind1, int ind2){
//--------------------------------------------------------------------------------------
  int i,j,k,l,m,*angm,*bdim,*cdim,nbond1,nten1,nmom,*angm1,*angm2;
  bool check;
  tensor **parr1,**pcgc1,*tarr1,*tcgc1,ta,tb,tc,td,te,tf;
  tensor_su2 tmp;
  su2struct cgc1;
  su2bond *bb,*bb2;
  if(ind2!=ind1+1||cgc.get_bonddir(ind1)!=cgc.get_bonddir(ind2)||ind2>=nbond){
    cout<<"tensor_su2::fuse: can not merge indices that are not next to each other or not have the same bond direction"<<endl;
    exit(0);
  }
  nbond1=nbond-1;
  bb=new su2bond[nbond1];
  bb2=new su2bond[2];
  get_su2bond(ind1,bb2[0]);
  get_su2bond(ind2,bb2[1]);
  bb2[0].invert_bonddir();
  bb2[1].invert_bonddir();
  tmp.fuse(bb2[0],bb2[1]);
  bb[ind1].fuse(bb2[0],bb2[1]);
  for(i=0;i<ind1;i++)
    get_su2bond(i,bb[i]);
  for(i=ind2;i<nbond-1;i++)
    get_su2bond(i+1,bb[i]);
  cgc1.set_su2struct(nbond1,locspin,bb);
  nten1=cgc1.get_nten();
  parr1=new tensor*[nten1];
  pcgc1=new tensor*[nten1];
  tarr1=new tensor[nten1];
  tcgc1=new tensor[nten1];
  for(i=0;i<nten1;i++){
    parr1[i]=&(tarr1[i]);
    pcgc1[i]=&(tcgc1[i]);
  }
  nmom=bb[ind1].get_nmoment();
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  angm1=new int[nbond1];
  angm2=new int[3];
  for(i=0;i<nten;i++){
    check=get_tensor_argument(i,angm,bdim,cdim);
    if(check==false)continue;
    for(j=0;j<ind1;j++)
      angm1[j]=angm[j];
    for(j=ind2;j<nbond-1;j++)
      angm1[j]=angm[j+1];
    angm2[0]=angm[ind1];
    angm2[1]=angm[ind2];
    for(j=0;j<nmom;j++){
      angm1[ind1]=bb[ind1].get_angularmoment(j);
      angm2[2]=angm1[ind1];
      k=cgc1.get_tensor_index(angm1);
      l=tmp.get_tensor_index(angm2);
      if(cgc1.check_angularmoments(angm1)==false)continue;
      if(this->is_null(i)||tmp.is_null(l))continue;
      te=*(parr[i]);
      tf=*(pcgc[i]);
      te.mergeindex(ind1,ind2);
      tf.mergeindex(ind1,ind2);
      ta=*(tmp.get_parr(l));
      tb=*(tmp.get_pcgc(l));
      ta.mergeindex(0,1);
      tb.mergeindex(0,1);
      tc.contract(ta,0,te,ind1);
      td.contract(tb,0,tf,ind1);
      tc.shift(0,ind1);
      td.shift(0,ind1);
      if(parr1[k]->is_null()){
	tarr1[k]=tc;
	tcgc1[k]=td;
      }
      else{
	if(!sum_direct_product(tarr1[k],tcgc1[k],tc,td)){
	  cout<<"tensor_su2::fuse, su2 tensor fuse: wrong"<<endl;
	  exit(0);
	}
      }
    }
  }
  nbond=nbond1;
  nten=nten1;
  cgc=cgc1;
  delete []tarr;
  delete []tcgc;
  delete []parr;
  delete []pcgc;
  parr=parr1;
  pcgc=pcgc1;
  tarr=tarr1;
  tcgc=tcgc1;
  delete []angm;
  delete []angm1;
  delete []angm2;
  delete []bdim;
  delete []cdim;
  delete []bb;
  delete []bb2;
  //make_standard_cgc();
}

//--------------------------------------------------------------------------------------
void tensor_su2::fuse_to_multiplet(int ind1, int ind2, int angt){
//--------------------------------------------------------------------------------------
  int i,j,k,l,m,*angm,*bdim,*cdim,nbond1,nten1,nmom,*angm1,*angm2;
  bool check;
  tensor **parr1,**pcgc1,*tarr1,*tcgc1,ta,tb,tc,td,te,tf;
  tensor_su2 tmp;
  su2struct cgc1;
  su2bond *bb,*bb2;
  if(ind2!=ind1+1||cgc.get_bonddir(ind1)!=cgc.get_bonddir(ind2)||ind2>=nbond){
    cout<<"tensor_su2::fuse: can not merge indices that are not next to each other or not have the same bond direction"<<endl;
    exit(0);
  }
  nbond1=nbond-1;
  bb=new su2bond[nbond1];
  bb2=new su2bond[2];
  get_su2bond(ind1,bb2[0]);
  get_su2bond(ind2,bb2[1]);
  bb2[0].invert_bonddir();
  bb2[1].invert_bonddir();
  tmp.fuse_to_multiplet(bb2[0],bb2[1],angt);
  bb[ind1].fuse_to_multiplet(bb2[0],bb2[1],angt);
  for(i=0;i<ind1;i++)
    get_su2bond(i,bb[i]);
  for(i=ind2;i<nbond-1;i++)
    get_su2bond(i+1,bb[i]);
  cgc1.set_su2struct(nbond1,locspin,bb);
  nten1=cgc1.get_nten();
  parr1=new tensor*[nten1];
  pcgc1=new tensor*[nten1];
  tarr1=new tensor[nten1];
  tcgc1=new tensor[nten1];
  for(i=0;i<nten1;i++){
    parr1[i]=&(tarr1[i]);
    pcgc1[i]=&(tcgc1[i]);
  }
  nmom=bb[ind1].get_nmoment();
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  angm1=new int[nbond1];
  angm2=new int[3];
  for(i=0;i<nten;i++){
    check=get_tensor_argument(i,angm,bdim,cdim);
    if(check==false)continue;
    for(j=0;j<ind1;j++)
      angm1[j]=angm[j];
    for(j=ind2;j<nbond-1;j++)
      angm1[j]=angm[j+1];
    angm2[0]=angm[ind1];
    angm2[1]=angm[ind2];
    for(j=0;j<nmom;j++){
      angm1[ind1]=bb[ind1].get_angularmoment(j);
      angm2[2]=angm1[ind1];
      k=cgc1.get_tensor_index(angm1);
      l=tmp.get_tensor_index(angm2);
      if(cgc1.check_angularmoments(angm1)==false)continue;
      if(this->is_null(i)||tmp.is_null(l))continue;
      te=*(parr[i]);
      tf=*(pcgc[i]);
      te.mergeindex(ind1,ind2);
      tf.mergeindex(ind1,ind2);
      ta=*(tmp.get_parr(l));
      tb=*(tmp.get_pcgc(l));
      ta.mergeindex(0,1);
      tb.mergeindex(0,1);
      tc.contract(ta,0,te,ind1);
      td.contract(tb,0,tf,ind1);
      tc.shift(0,ind1);
      td.shift(0,ind1);
      if(parr1[k]->is_null()){
	tarr1[k]=tc;
	tcgc1[k]=td;
      }
      else{
	if(!sum_direct_product(tarr1[k],tcgc1[k],tc,td)){
	  cout<<"tensor_su2::fuse, su2 tensor fuse: wrong"<<endl;
	  exit(0);
	}
      }
    }
  }
  nbond=nbond1;
  nten=nten1;
  cgc=cgc1;
  delete []tarr;
  delete []tcgc;
  delete []parr;
  delete []pcgc;
  parr=parr1;
  pcgc=pcgc1;
  tarr=tarr1;
  tcgc=tcgc1;
  delete []angm;
  delete []angm1;
  delete []angm2;
  delete []bdim;
  delete []cdim;
  delete []bb;
  delete []bb2;
  make_standard_cgc();
}

//--------------------------------------------------------------------------------------
void tensor_su2::make_standard_cgc(){
//--------------------------------------------------------------------------------------
  int i,j,k,l,*angm,*bdim,*cdim;
  bool check;
  tensor tmp;
  double nor;
  if(nbond!=3&&nbond!=2){
    cout<<"tensor_su2::make_standard_cgc not necessary nbond="<<nbond<<endl;    
    return;
    //exit(0);
  }
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  for(i=0;i<nten;i++){
    check=get_tensor_argument(i,angm,bdim,cdim);
    if(check==false)continue;
    if(nbond==3)
      tmp.make_cgc(angm[0],angm[1],angm[2]);
    else if(nbond==2&&angm[0]==angm[1])
      tmp.make_identity(angm[0]);
    else{
      cout<<"tensor_su2::make_standard_cgc cgc is not correct"<<endl;
      pcgc[i]->print();
      //parr[i]->print();
      exit(0);
    }
    if(pcgc[i]->is_proportional_to(tmp,nor)){
      *(pcgc[i])=tmp;
      *(parr[i])*=nor;
    }    
    else if(pcgc[i]->is_zero()){
      *(pcgc[i])=tmp;
      *(parr[i])*=0;
    }    
    else if(parr[i]->is_zero()){
      *(pcgc[i])=tmp;
      *(parr[i])*=0;
    }    
    else{
      cout<<"tensor_su2::make_standard_cgc cgc is not correct"<<endl;
      pcgc[i]->print();
      tmp.print();
      exit(0);
    }
  }
  delete []angm;
  delete []bdim;
  delete []cdim;
}

//--------------------------------------------------------------------------------------
void tensor_su2::fuse(su2bond& b1, su2bond& b2, int dir){
//--------------------------------------------------------------------------------------
  int i,j,k,l,n0,n1,n2,n3,m0,m1,m2,m3,ms,me,d0,d1,d2,d3;
  int *bdim,*pos,*bdim2;
  su2bond *bb;
  clean();
  nbond=3;

  bb=new su2bond[nbond];
  bb[1]=b1;
  bb[2]=b2;
  bb[0].fuse(bb[1],bb[2],dir);

  locspin=0;
  cgc.set_su2struct(3,0,bb);
  nten=cgc.get_nten();

  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];

  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }

  n0=bb[0].get_nmoment();
  n1=bb[1].get_nmoment();
  n2=bb[2].get_nmoment();

  pos=new int[n0];
  bdim=new int[nbond];
  bdim2=new int[nbond];

  for(i=0;i<n0;i++)
    pos[i]=0;

  for(k=0;k<n2;k++){
    m2=bb[2].get_angularmoment(k);
    d2=bb[2].get_bonddim(k);
    bdim[2]=d2;
    for(j=0;j<n1;j++){
      m1=bb[1].get_angularmoment(j);
      d1=bb[1].get_bonddim(j);
      bdim[1]=d1;
      ms=abs(m1-m2);
      me=m1+m2;
      for(m0=ms;m0<=me;m0+=2){
	i=bb[0].get_angularmoment_index(m0);
	d0=bb[0].get_bonddim(i);
	bdim[0]=d0;
	l=i+j*n0+k*n0*n1;
	pcgc[l]->make_cgc(m0,m1,m2);
	bdim2[0]=bdim[1];
	bdim2[1]=bdim[2];
	bdim2[2]=bdim[0];
	parr[l]->shift_set_identity(3,pos[i],bdim2);
	parr[l]->shift(2,0);
	pos[i]+=d1*d2;
      }
    }
  }
  for(i=0;i<n0;i++){
    if(bb[0].get_bonddim(i)!=pos[i]){
      cout<<"tensor_su2::fuse, wrong in fuse two su2 bond"<<endl;
      exit(0);
    }
  }
  delete []pos;
  delete []bdim;
  delete []bdim2;
  delete []bb;
}

//--------------------------------------------------------------------------------------
void tensor_su2::fuse(int ind1, int ind2, int ind3, int ind4){
//--------------------------------------------------------------------------------------
  int i,j,k,l,m,n,p,i0,*angm,*bdim,*cdim,nbond1,nten1,nmom1,nmom2,*angm1,*angm2,*angm3,*angmf1,*angmf2,*angmt,*id;
  int nb1,nb2;
  bool check;
  tensor **parr1,**pcgc1,*tarr1,*tcgc1,ta,tb,tc,td,te,tf,tg,th,tmp1,tmp3,tmp4;
  tensor_su2 *tt1,*tt2;
  su2struct cgc1;
  su2bond *bb1,*bb2,*bb,*bbf1,*bbf2;
  if(!(ind2>ind1&&ind3>ind2&&ind4>ind3)){
    cout<<"can not fuse from ind1 to ind2, from ind3 to ind4"<<endl;
    exit(0);
  }
  nb1=ind2-ind1+1;
  nb2=ind4-ind3+1;
  bb1=new su2bond[nb1+1];
  bb2=new su2bond[nb2+1];
  bbf1=new su2bond[nb1];
  bbf2=new su2bond[nb2];
  tt1=new tensor_su2[nb1-1];
  tt2=new tensor_su2[nb2-1];
  angmf1=new int[nb1-1];
  angmf2=new int[nb2-1];
  id=new int[nb1+nb2];

  for(i=ind1;i<=ind2;i++)
    get_su2bond(i,bb1[i-ind1]);
  for(i=ind3;i<=ind4;i++)
    get_su2bond(i,bb2[i-ind3]);

  bbf1[0]=bb1[0];
  for(i=0;i<nb1-1;i++){
    tt1[i].fuse(bbf1[i],bb1[i+1]);
    tt1[i].get_su2bond(2,bbf1[i+1]);
    bbf1[i+1].invert_bonddir();
  }
  bbf2[0]=bb2[0];
  for(i=0;i<nb2-1;i++){
    tt2[i].fuse(bbf2[i],bb2[i+1]);
    tt2[i].get_su2bond(2,bbf2[i+1]);
    bbf2[i+1].invert_bonddir();
  }
  nbond1=nbond-(ind2-ind1)-(ind4-ind3);
  bb=new su2bond[nbond1];
  for(i=0;i<ind1;i++)
    get_su2bond(i,bb[i]);
  for(i=ind1+1;i<ind3-(ind2-ind1);i++)
    get_su2bond(i+(ind2-ind1),bb[i]);
  for(i=ind3-(ind2-ind1)+1;i<nbond1;i++)
    get_su2bond(i+(ind2-ind1)+(ind4-ind3),bb[i]);
  tt1[nb1-2].get_su2bond(2,bb[ind1]);
  tt2[nb2-2].get_su2bond(2,bb[ind3-(ind2-ind1)]);
  bb[ind1].invert_bonddir();
  bb[ind3-(ind2-ind1)].invert_bonddir();
  cgc1.set_su2struct(nbond1,locspin,bb);
  delete []bbf1;
  delete []bbf2;
  delete []bb1;
  delete []bb2;
  delete []bb;
  nten1=cgc1.get_nten();
  parr1=new tensor*[nten1];
  pcgc1=new tensor*[nten1];
  tarr1=new tensor[nten1];
  tcgc1=new tensor[nten1];
  for(i=0;i<nten1;i++){
    parr1[i]=&(tarr1[i]);
    pcgc1[i]=&(tcgc1[i]);
  }
  nmom1=1;
  for(i=0;i<nb1-1;i++)
    nmom1*=tt1[i].get_nmoment(2);
  nmom2=1;
  for(i=0;i<nb2-1;i++)
    nmom2*=tt2[i].get_nmoment(2);

  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  angm1=new int[nb1+1];
  angm2=new int[nb2+1];
  angm3=new int[nbond1];
  angmf1=new int[nb1-1];
  angmf2=new int[nb2-1];
  angmt=new int[3];
  for(i=0;i<nten;i++){
    check=get_tensor_argument(i,angm,bdim,cdim);
    if(check==false||this->is_null(i))continue;
    for(j=0;j<ind1;j++)
      angm3[j]=angm[j];
    for(j=ind1+1;j<ind3-(ind2-ind1);j++)
      angm3[j]=angm[j+(ind2-ind1)];
    for(j=ind3-(ind2-ind1)+1;j<nbond1;j++)
      angm3[j]=angm[j+(ind2-ind1)+(ind4-ind3)];
    for(j=0;j<nb1;j++)
      angm1[j]=angm[ind1+j];
    for(j=0;j<nb2;j++)
      angm2[j]=angm[ind3+j];
    for(j=0;j<nmom1;j++){
      i0=j;
      for(l=0;l<nb1-1;l++){
	angmf1[l]=tt1[l].get_angularmoment(2,i0%tt1[l].get_nmoment(2));
	i0/=tt1[l].get_nmoment(2);
      }
      angm3[ind1]=angmf1[nb1-2];
      angmt[0]=angm1[0];
      angmt[1]=angm1[1];
      angmt[2]=angmf1[0];
      id[0]=tt1[0].get_tensor_index(angmt);
      check=tt1[0].check_angularmoments(angmt);
      for(l=1;l<nb1-1;l++){
	angmt[0]=angmt[2];
	angmt[1]=angm1[l+1];
	angmt[2]=angmf1[l];
	id[l]=tt1[l].get_tensor_index(angmt);
	check=check&(tt1[l].check_angularmoments(angmt));
      }
      if(check==false)continue;
      ta=*(tt1[0].get_parr(id[0]));
      tb=*(tt1[0].get_pcgc(id[0]));
      ta.mergeindex(0,1);
      tb.mergeindex(0,1);
      for(l=1;l<nb1-1;l++){
	tmp1.contract(ta,1,*(tt1[l].get_parr(id[l])),0);
	tmp1.mergeindex(0,1);
	ta=tmp1;
	tmp1.contract(tb,1,*(tt1[l].get_pcgc(id[l])),0);
	tmp1.mergeindex(0,1);
	tb=tmp1;
      }
      for(k=0;k<nmom2;k++){
	i0=k;
	for(l=0;l<nb2-1;l++){
	  angmf2[l]=tt2[l].get_angularmoment(2,i0%tt2[l].get_nmoment(2));
	  i0/=tt2[l].get_nmoment(2);
	}
	angm3[ind3-(ind2-ind1)]=angmf2[nb2-2];
	angmt[0]=angm2[0];
	angmt[1]=angm2[1];
	angmt[2]=angmf2[0];
	id[0]=tt2[0].get_tensor_index(angmt);
	check=tt2[0].check_angularmoments(angmt);
	for(l=1;l<nb2-1;l++){
	  angmt[0]=angmt[2];
	  angmt[1]=angm2[l+1];
	  angmt[2]=angmf2[l];
	  id[l]=tt2[l].get_tensor_index(angmt);
	  check=check&(tt2[l].check_angularmoments(angmt));
	}
	if(check==false)continue;
	tc=*(tt2[0].get_parr(id[0]));
	td=*(tt2[0].get_pcgc(id[0]));
	tc.mergeindex(0,1);
	td.mergeindex(0,1);
	for(l=1;l<nb2-1;l++){
	  tmp1.contract(tc,1,*(tt2[l].get_parr(id[l])),0);
	  tmp1.mergeindex(0,1);
	  tc=tmp1;
	  tmp1.contract(td,1,*(tt2[l].get_pcgc(id[l])),0);
	  tmp1.mergeindex(0,1);
	  td=tmp1;
	}
	m=cgc1.get_tensor_index(angm3);
	if(cgc1.check_angularmoments(angm3)==false)continue;
	te=*(parr[i]);
	tf=*(pcgc[i]);
	for(p=1;p<nb2;p++){
	  te.mergeindex(ind3,ind3+1);
	  tf.mergeindex(ind3,ind3+1);
	}
	for(p=1;p<nb1;p++){
	  te.mergeindex(ind1,ind1+1);
	  tf.mergeindex(ind1,ind1+1);
	}
	tmp3.contract(ta,0,te,ind1);
	tmp4.contract(tb,0,tf,ind1);
	tg.contract(tc,0,tmp3,ind3-ind2);
	th.contract(td,0,tmp4,ind3-ind2);
	tg.shift(0,ind3-(ind2-ind1));
	th.shift(0,ind3-(ind2-ind1));
	if(parr1[m]->is_null()){
	  tarr1[m]=tg;
	  tcgc1[m]=th;
	}
	else{
	  if(!sum_direct_product(tarr1[m],tcgc1[m],tg,th)){
	    cout<<"tensor_su2::fuse, su2 tensor fuse: wrong"<<endl;
	    exit(0);
	  }
	}
      }
    }
  }
  nbond=nbond1;
  nten=nten1;
  cgc=cgc1;
  delete []tarr;
  delete []tcgc;
  delete []parr;
  delete []pcgc;
  parr=parr1;
  pcgc=pcgc1;
  tarr=tarr1;
  tcgc=tcgc1;
  delete []angm;
  delete []angm1;
  delete []angm2;
  delete []angm3;
  delete []angmf1;
  delete []angmf2;
  delete []angmt;
  delete []bdim;
  delete []cdim;
  delete []id;
  delete []tt1;
  delete []tt2;
}

//--------------------------------------------------------------------------------------
void tensor_su2::operator_tensor_product_identity(tensor_su2& op, su2bond& bd_idn){
//--------------------------------------------------------------------------------------
  int i,j,*angm,*bdim,*cdim,nbond1,nten1;
  su2bond *bb;
  bool check;
  tensor iden;
  clean();
  nbond1=op.get_nbond();
  nbond=nbond1+2;
  bb=new su2bond[nbond];
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  locspin=op.get_locspin();
  for(i=0;i<nbond1;i++)
    op.get_su2bond(i,bb[i]);
  bb[nbond1]=bd_idn;
  bb[nbond1+1]=bd_idn;
  bb[nbond1].set_bonddir(bb[nbond1-1].get_bonddir());
  bb[nbond1+1].set_bonddir(-bb[nbond1-1].get_bonddir());
  cgc.set_su2struct(nbond,locspin,bb);
  nten=cgc.get_nten();
  delete []bb;
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    if(check&&angm[nbond1]==angm[nbond1+1]){
      j=op.get_tensor_index(angm);
      if(op.get_parr(j)->is_null())continue;
      iden.make_identity(bdim[nbond1]-1);
      tarr[i].tensor_product(*(op.get_parr(j)),iden);
      tcgc[i].tensor_product(*(op.get_pcgc(j)),identity[angm[nbond1]]);
    }
  }
  delete []angm;
  delete []bdim;
  delete []cdim;
}

//--------------------------------------------------------------------------------------
void tensor_su2::operator_tensor_product_identity(tensor_su2& op, tensor_su2& ten_idn){
//--------------------------------------------------------------------------------------
  int i,j,k,*angm,*bdim,*cdim,nbond1,nten1;
  su2bond *bb;
  bool check;
  clean();
  nbond1=op.get_nbond();
  nbond=nbond1+2;
  bb=new su2bond[nbond];
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  locspin=op.get_locspin();
  for(i=0;i<nbond1;i++)
    op.get_su2bond(i,bb[i]);
  ten_idn.get_su2bond(0,bb[nbond1]);;
  ten_idn.get_su2bond(1,bb[nbond1+1]);;
  bb[nbond1].set_bonddir(bb[nbond1-1].get_bonddir());
  bb[nbond1+1].set_bonddir(-bb[nbond1-1].get_bonddir());
  cgc.set_su2struct(nbond,locspin,bb);
  nten=cgc.get_nten();
  delete []bb;
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    if(check&&angm[nbond1]==angm[nbond1+1]){
      j=op.get_tensor_index(angm);
      k=ten_idn.get_tensor_index(&(angm[nbond1]));
      if(op.get_parr(j)->is_null()||ten_idn.get_parr(k)->is_null())continue;
      tarr[i].tensor_product(*(op.get_parr(j)),*(ten_idn.get_parr(k)));
      tcgc[i].tensor_product(*(op.get_pcgc(j)),*(ten_idn.get_pcgc(k)));
    }
  }
  delete []angm;
  delete []bdim;
  delete []cdim;
}

//--------------------------------------------------------------------------------------
void tensor_su2::fuse_to_multiplet(su2bond& b1, su2bond& b2, int angm){
//--------------------------------------------------------------------------------------
  int i,j,k,l,n1,n2,n3,m1,m2,m3,ms,me,d1,d2,d3;
  int *bdim,*pos;
  su2bond *bb;
  clean();
  nbond=3;

  if(b1.get_bonddir()!=b2.get_bonddir()){
    cout<<"tensor_su2::fuse two bonds had different directions"<<endl;
    exit(0);
  }
  bb=new su2bond[nbond];
  bb[0]=b1;
  bb[1]=b2;
  bb[2].fuse_to_multiplet(bb[0],bb[1],angm);

  locspin=0;
  cgc.set_su2struct(3,0,bb);
  nten=cgc.get_nten();

  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];

  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }

  n1=bb[0].get_nmoment();
  n2=bb[1].get_nmoment();
  n3=bb[2].get_nmoment();

  pos=new int[n3];
  bdim=new int[nbond];

  for(i=0;i<n3;i++)
    pos[i]=0;

  for(j=0;j<n2;j++){
    m2=bb[1].get_angularmoment(j);
    d2=bb[1].get_bonddim(j);
    bdim[1]=d2;
    for(i=0;i<n1;i++){
      m1=bb[0].get_angularmoment(i);
      d1=bb[0].get_bonddim(i);
      bdim[0]=d1;
      ms=abs(m1-m2);
      me=m1+m2;
      for(m3=ms;m3<=me;m3+=2){
	if(m3!=angm)continue;
	k=bb[2].get_angularmoment_index(m3);
	d3=bb[2].get_bonddim(k);
	bdim[2]=d3;
	l=i+j*n1+k*n1*n2;
	pcgc[l]->make_cgc(m1,m2,m3);
	parr[l]->shift_set_identity(3,pos[k],bdim);
	pos[k]+=d1*d2;
      }
    }
  }
  for(i=0;i<n3;i++){
    if(bb[2].get_bonddim(i)!=pos[i]){
      cout<<"tensor_su2::fuse, wrong in fuse two su2 bond"<<endl;
      exit(0);
    }
  }
  delete []pos;
  delete []bdim;
  delete []bb;
}

//--------------------------------------------------------------------------------------
void tensor_su2::fuse_to_singlet(su2bond& b1, su2bond& b2){
//--------------------------------------------------------------------------------------
  int i,j,k,l,n1,n2,n3,m1,m2,m3,ms,me,d1,d2,d3;
  int *bdim,*pos;
  su2bond *bb;
  clean();
  nbond=2;

  if(b1.get_bonddir()!=b2.get_bonddir()){
    cout<<"tensor_su2::fuse two bonds had different directions"<<endl;
    exit(0);
  }
  bb=new su2bond[nbond];
  bdim=new int[2];
  bb[0]=b1;
  bb[1]=b2;

  locspin=0;
  cgc.set_su2struct(2,0,bb);
  nten=cgc.get_nten();

  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];

  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }

  n1=bb[0].get_nmoment();
  n2=bb[1].get_nmoment();

  for(j=0;j<n2;j++){
    m2=bb[1].get_angularmoment(j);
    bdim[1]=bb[1].get_bonddim(j);
    for(i=0;i<n1;i++){
      m1=bb[0].get_angularmoment(i);
      bdim[0]=bb[0].get_bonddim(i);
      if(m1==m2){
	l=i+j*n1;
	pcgc[l]->make_singlet(m1);
	parr[l]->alloc_space(2,bdim);
	(*parr[l])=1;
      }
    }
  }
  delete []bdim;
  delete []bb;
}

//--------------------------------------------------------------------------------------
void tensor_su2::contract(tensor_su2& t1, int i1, tensor_su2& t2, int i2){
//--------------------------------------------------------------------------------------
  su2bond *bb,*bb2;
  int i,j,k,l,i0,nbond1,nbond2,nmom;
  int *angm,*bdim,*cdim,*angm1,*angm2;
  tensor tmp1,tmp2;
  bool check,check1,check2;
  double nor1,nor2;

  bb2=new su2bond[2];
  clean();
  t1.get_su2bond(i1,bb2[0]);
  t2.get_su2bond(i2,bb2[1]);
  bb2[1].invert_bonddir();
  if(bb2[0]!=bb2[1]){
    if(comm_rank==0){
      cout<<"tensor_su2::contract, two su2 bonds can not be contracted, check su2bond parameters"<<endl;
      bb2[0].print();
      bb2[1].print();
      t1.print();
      t2.print();
    }
    exit(0);
  }
  locspin=t1.get_locspin()+t2.get_locspin();
  nbond1=t1.get_nbond();
  nbond2=t2.get_nbond();
  nbond=nbond1+nbond2-2;
  bb=new su2bond[nbond];
  for(i=0;i<nbond1-1;i++)
    t1.get_su2bond((i1+1+i)%nbond1,bb[i]);
  for(i=0;i<nbond2-1;i++)
    t2.get_su2bond((i2+1+i)%nbond2,bb[i+nbond1-1]);
  cgc.set_su2struct(nbond,locspin,bb);
  nten=cgc.get_nten();
  angm=new int[nbond];
  angm1=new int[nbond1];
  angm2=new int[nbond2];
  bdim=new int[nbond];
  cdim=new int[nbond];
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }
  nmom=bb2[0].get_nmoment();
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    for(k=0;k<nbond1-1;k++)
      angm1[(i1+1+k)%nbond1]=angm[k];
    for(k=0;k<nbond2-1;k++)
      angm2[(i2+1+k)%nbond2]=angm[k+nbond1-1];
    for(j=0;j<nmom;j++){
      angm1[i1]=bb2[0].get_angularmoment(j);
      angm2[i2]=angm1[i1];
      check1=t1.check_angularmoments(angm1);
      check2=t2.check_angularmoments(angm2);
      if(check==false||check1==false||check2==false)continue;
      k=t1.get_tensor_index(angm1);
      l=t2.get_tensor_index(angm2);
      if(t1.get_parr(k)->is_null()||t2.get_parr(l)->is_null())continue;
      if(tarr[i].is_null()){
	tarr[i].contract(*(t1.get_parr(k)),i1,*(t2.get_parr(l)),i2);
	tcgc[i].contract(*(t1.get_pcgc(k)),i1,*(t2.get_pcgc(l)),i2);
      }
      else{
	tmp1.contract(*(t1.get_parr(k)),i1,*(t2.get_parr(l)),i2);
	tmp2.contract(*(t1.get_pcgc(k)),i1,*(t2.get_pcgc(l)),i2);
	if(!sum_direct_product(tarr[i],tcgc[i],tmp1,tmp2)){
	  cout<<"tensor_su2::contract, su2 tensor contraction: not getting unitary tensors"<<endl;
	  tarr[i].print();
	  tcgc[i].print();
	  tmp1.print();
	  tmp2.print();
	  exit(0);
	}
      }
    }
  }
  delete []angm;
  delete []angm1;
  delete []angm2;
  delete []bdim;
  delete []cdim;
  delete []bb;
  delete []bb2;
}

//--------------------------------------------------------------------------------------
void tensor_su2::shift(int i0, int i1){
//--------------------------------------------------------------------------------------
  int i,j,ishift,*angm,*bdim,*cdim,*angm1;
  bool check;
  tensor **parr1,**pcgc1;
  su2struct cgc1;
  if(i0==i1)return;
  if(i1>i0)ishift=i1-i0;
  else ishift=nbond-(i0-i1);
  cgc1=cgc;
  cgc1.shift(i0,i1);
  parr1=new tensor*[nten];
  pcgc1=new tensor*[nten];
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  angm1=new int[nbond];
  for(i=0;i<nten;i++){
    check=get_tensor_argument(i,angm,bdim,cdim);
    for(j=0;j<nbond;j++)
      angm1[(j+ishift)%nbond]=angm[j];
    j=cgc1.get_tensor_index(angm1);
    parr1[j]=parr[i];
    pcgc1[j]=pcgc[i];
    if(!parr1[j]->is_null()){
      parr1[j]->shift(i0,i1);
      pcgc1[j]->shift(i0,i1);
    }
  }
  cgc=cgc1;
  delete []parr;
  delete []pcgc;
  parr=parr1;
  pcgc=pcgc1;
  delete []angm;
  delete []angm1;
  delete []bdim;
  delete []cdim;
}

//--------------------------------------------------------------------------------------
void tensor_su2::exchangeindex(int ind1, int ind2){
//--------------------------------------------------------------------------------------
  int i,j,*angm,*bdim,*cdim;
  bool check;
  tensor **parr1,**pcgc1;
  su2struct cgc1;

  if(ind1==ind2)return;

  cgc1=cgc;
  cgc1.exchangeindex(ind1,ind2);

  parr1=new tensor*[nten];
  pcgc1=new tensor*[nten];
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  for(i=0;i<nten;i++){
    check=get_tensor_argument(i,angm,bdim,cdim);
    j=angm[ind1];
    angm[ind1]=angm[ind2];
    angm[ind2]=j;
    j=cgc1.get_tensor_index(angm);
    parr1[j]=parr[i];
    pcgc1[j]=pcgc[i];
    if(!parr1[j]->is_null()){
      parr1[j]->exchangeindex(ind1,ind2);
      pcgc1[j]->exchangeindex(ind1,ind2);
    }
  }
  cgc=cgc1;
  delete []parr;
  delete []pcgc;
  parr=parr1;
  pcgc=pcgc1;
  delete []angm;
  delete []bdim;
  delete []cdim;
}

//--------------------------------------------------------------------------------------
void tensor_su2::set_tensor_su2(int nb, int locsp, su2bond* bb, tensor **pa){
//--------------------------------------------------------------------------------------
  bool check;
  int i,j,*angm,*bdim,*cdim;
  clean();
  nbond=nb;
  locspin=locsp;
  cgc.set_su2struct(nbond,locspin,bb);
  nten=cgc.get_nten();
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }  
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  j=0;
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    if(check){
      tarr[i]=*(pa[j]);
      tcgc[i].make_identity(angm[0]);
      j++;
    }
  }
  delete []angm;
  delete []bdim;
  delete []cdim;
}

//--------------------------------------------------------------------------------------
void tensor_su2::set_tensor_su2(su2struct& cgc1, double *tele1, double* tele2,int* iflag){
//--------------------------------------------------------------------------------------
  bool check;
  int i,j1,j2,*angm,*bdim,*cdim;
  clean();
  cgc=cgc1;
  nbond=cgc.get_nbond();
  nten=cgc.get_nten();
  locspin=cgc.get_locspin();
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }  
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  j1=0;
  j2=0;
  for(i=0;i<nten;i++){
    if(iflag[i]){
      check=cgc.get_tensor_argument(i,angm,bdim,cdim);
      parr[i]->copy(nbond,bdim,&(tele1[j1]));
      pcgc[i]->copy(nbond,cdim,&(tele2[j2]));
      j1+=parr[i]->get_nelement();
      j2+=pcgc[i]->get_nelement();
    }
  }
  delete []angm;
  delete []bdim;
  delete []cdim;
}

//--------------------------------------------------------------------------------------
void tensor_su2::set_tensor_su2(su2struct& cgc1, double *tele1, double* tele2){
//--------------------------------------------------------------------------------------
  bool check;
  int i,j1,j2,*angm,*bdim,*cdim;
  clean();
  cgc=cgc1;
  nbond=cgc.get_nbond();
  nten=cgc.get_nten();
  locspin=cgc.get_locspin();
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }  
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  j1=0;
  j2=0;
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    if(check){
      parr[i]->copy(nbond,bdim,&(tele1[j1]));
      pcgc[i]->copy(nbond,cdim,&(tele2[j2]));
      j1+=parr[i]->get_nelement();
      j2+=pcgc[i]->get_nelement();
    }
  }
  delete []angm;
  delete []bdim;
  delete []cdim;
}

//--------------------------------------------------------------------------------------
void tensor_su2::set_tensor_su2(su2struct& cgc1, double *tele1){
//--------------------------------------------------------------------------------------
  bool check;
  int i,j1,j2,*angm,*bdim,*cdim;
  clean();
  cgc=cgc1;
  nbond=cgc.get_nbond();
  nten=cgc.get_nten();
  locspin=cgc.get_locspin();
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }  
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  j1=0;
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    if(check){
      parr[i]->copy(nbond,bdim,&(tele1[j1]));
      j1+=parr[i]->get_nelement();
    }
  }
  if(nbond==3)
    cgc_make_cgc();
  else if(nbond==2&&get_bonddir(0)==1&&get_bonddir(1)==-1)
    cgc_make_identity_cgc();
  else{
    cout<<"can not do set_tensor_su2"<<endl;
    this->print();
    exit(0);
  }
  delete []angm;
  delete []bdim;
  delete []cdim;
}

//--------------------------------------------------------------------------------------
void tensor_su2::cgc_make_identity_cgc(){
//--------------------------------------------------------------------------------------
  if(nbond!=2){
    cout<<"tensor_su2::cgc_make_identity_cgc, can not perform"<<endl;
    exit(0);
  }
  int i,*angm,*bdim,*cdim;
  bool check;
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  for(i=0;i<nten;i++){
    check=cgc.get_tensor_argument(i,angm,bdim,cdim);
    if(check)
      pcgc[i]->make_identity(angm[0]);
  }
  delete []angm;
  delete []bdim;
  delete []cdim;
}

//--------------------------------------------------------------------------------------
int tensor_su2::svd(tensor_su2& tu, double p1, tensor_su2& tv, double p2, double *wout){
//--------------------------------------------------------------------------------------
  int i,j,k,m,n,n0,a1,a2,m1,m2,nmom1,nmom2,nmom3,*dc,*dim,*bdim,*isort,*amom;
  int nm1,nm2;
  tensor tmp,*uu,*vv,**puarr,**pvarr;
  tensor_su2 proju,projv,tmp_su2,left,rght,tmp1,tmp2;
  double **ww,*wsort;
  su2bond bb[2],bb2[4],bb3[2];
  //print();
  nm1=tu.get_nmoment(2);
  tu.get_su2bond(0,bb2[0]);
  tu.get_su2bond(1,bb2[1]);
  tv.get_su2bond(0,bb2[2]);
  tv.get_su2bond(1,bb2[3]);

  bb[0].fuse(bb2[0],bb2[1]);
  bb[1].fuse(bb2[2],bb2[3]);

  nmom1=bb[0].get_nmoment();
  nmom2=bb[1].get_nmoment();
  uu=new tensor[nmom1];
  vv=new tensor[nmom1];
  dc=new int[nmom1];
  dim=new int[nmom1];
  bdim=new int[nmom1+2];
  amom=new int[nmom1];
  ww=new double*[nmom1];
  for(i=0;i<nmom1;i++)
    ww[i]=new double[max_dcut];
  wsort=new double[max_dcut*nmom1];
  isort=new int[max_dcut*nmom1];
  for(i=0;i<nmom1;i++){
    dc[i]=0;
    dim[i]=0;
    bdim[i]=0;
  }

  n=0;
  n0=0;
  for(i=0;i<nmom1;i++){
    a1=bb[0].get_angularmoment(i);
    m1=bb[0].get_bonddim(i);
    for(j=0;j<nmom2;j++){
      a2=bb[1].get_angularmoment(j);
      m2=bb[1].get_bonddim(j);
      if(a1!=a2)continue;
      bdim[0]=m1;
      bdim[1]=m2;
      //tmp.copy(2,bdim,parr[0]->getptr(n0));
      tmp.copy(2,bdim,parr[i+j*nmom1]->getptr());
      tmp.svd(uu[i],p1,vv[i],p2,ww[i],dc[i],1);
      for(k=0;k<dc[i];k++){
	wsort[n+k]=ww[i][k];
	isort[n+k]=i;
	//cout<<i<<"\t"<<k<<"\t"<<wsort[n+k]<<endl;
      }
      n0+=m1*m2;
      n+=dc[i];
    }
  }
  dsort2_(&n,wsort,isort);

  if(max_dcut>n)k=0;
  else k=n-max_dcut;
  for(i=n-1;i>=k;i--){
    for(j=0;j<nmom1;j++)
      if(isort[i]==j){
	dim[j]++;
	break;
      }
  }
  nmom3=0;
  for(i=0;i<nmom1;i++)
    if(dim[i]!=0)
      nmom3++;
  //cout<<"nmom3="<<nmom3<<"\tnmom1="<<nmom1<<"\tnmom2="<<nmom2<<endl;
  puarr=new tensor*[nmom3];
  pvarr=new tensor*[nmom3];
  j=0;
  m=0;
  m1=max_dcut;
  for(i=0;i<nmom1;i++){
    if(dim[i]>0&&dim[i]<dc[i]){
      uu[i].direct_subtract(1,dc[i]-dim[i],tmp);
      vv[i].direct_subtract(1,dc[i]-dim[i],tmp);
    }
    if(dim[i]>0){
      puarr[j]=&(uu[i]);
      pvarr[j]=&(vv[i]);
      bdim[j]=dim[i];
      amom[j]=bb[0].get_angularmoment(i);
      for(k=0;k<dim[i];k++)
	wout[m+k]=ww[i][dc[i]-dim[i]+k];
      for(k=0;k<dc[i]-dim[i];k++)
	wout[m1+k]=ww[i][k];
      m+=dim[i];
      m1+=dc[i]-dim[i];
      j++;
    }
    else{
      for(k=0;k<dc[i];k++)
	wout[m1+k]=ww[i][k];
      m1+=dc[i];
    }
  }
  for(k=m1;k<max_dcut*2;k++)
    wout[k]=0;

  bb3[1].set_su2bond(nmom3,-1,amom,bdim);
  bb3[0]=bb[0];
  bb3[0].invert_bonddir();
  proju.set_tensor_su2(2,0,bb3,puarr);
  bb3[0]=bb[1];
  bb3[0].invert_bonddir();
  projv.set_tensor_su2(2,0,bb3,pvarr);
  left.fuse(bb2[0],bb2[1]);
  rght.fuse(bb2[2],bb2[3]);
  tu.contract(left,2,proju,0);
  tv.contract(rght,2,projv,0);
  nm2=tu.get_nmoment(2);
  /*
  tmp1.overlap_initial(tu,left,0);
  tmp2.overlap_initial(tv,rght,1);
  tmp2.take_conjugate(0);
  tmp_su2.contract(tmp1,0,tmp2,0);
  tmp_su2.normalize_vector();
  //this->print();
  //cout<<"+++++++++++++++"<<endl;
  //tmp_su2.print();
  if(*this!=tmp_su2)
    cout<<"svd compare tmp_su2 and vec not consistent"<<endl;
  else
    cout<<"svd compare tmp_su2 and vec consistent"<<endl;
  exit(0);
  projv.take_conjugate(1);
  tmp_su2.contract(proju,1,projv,1);
  if(*this!=tmp_su2)
    cout<<"svd compare tmp_su2 and vec not consistent"<<endl;
  else
    cout<<"svd compare tmp_su2 and vec consistent"<<endl;
  exit(0);
  */
  delete []uu;
  delete []vv;
  delete []dc;
  delete []dim;
  delete []bdim;
  delete []amom;
  for(i=0;i<nmom1;i++)
    delete []ww[i];
  delete []ww;
  delete []wsort;
  delete []isort;
  delete []puarr;
  delete []pvarr;
  if(nm1!=nm2)return 1;
  else return 0;
}

//--------------------------------------------------------------------------------------
bool sum_direct_product(tensor& ta1,tensor& tb1,tensor& ta2,tensor& tb2){
//--------------------------------------------------------------------------------------
  double nor1,nor2;
  if(tb1.is_zero()){
    ta1=ta2;
    tb1=tb2;
  }
  else if(tb2.is_zero())
    return true;
  else if(tb1==tb2)
    ta1+=ta2;
  else if(ta1==ta2)
    tb1+=tb2;
  else if(tb1.is_proportional_to(tb2,nor1)){
    ta1*=nor1;
    ta1+=ta2;
  }
  else if(ta1.is_proportional_to(ta2,nor1)){
    tb1*=nor1;
    tb1+=tb2;
  }
  else{
    cout<<"sum_direct_product(1,2,3,4) can not perform"<<endl;
    ta1.print();
    tb1.print();
    ta2.print();
    tb2.print();
    exit(0);
    return false;
  }
  return true;
}

//--------------------------------------------------------------------------------------
void sum_direct_product(tensor& ta1, tensor& tb1, tensor& ta2, tensor& tb2, tensor& ta3, tensor& tb3){
//--------------------------------------------------------------------------------------
  double nor1,nor2;
  bool pass1,pass2;
  pass1=sum_direct_product(ta1,tb1,ta2,tb2);  
  pass2=sum_direct_product(ta1,tb1,ta3,tb3);  
  if(!pass1||!pass2){
    cout<<"sum_direct_product can not sum two direct product tensors"<<endl;
    ta1.print();
    ta2.print();
    ta3.print();
    tb1.print();
    tb2.print();
    tb3.print();
    exit(0);
  }
}

//--------------------------------------------------------------------------------------
void tensor_su2::make_spinor_start(int physpn){
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
  angm[0]=physpn;
  bdim[0]=1;
  bb[0].set_su2bond(1,-1,angm,bdim);
  angm[0]=2;
  bdim[0]=1;
  bb[1].set_su2bond(1,-1,angm,bdim);
  angm[0]=physpn;
  bdim[0]=1;
  bb[2].set_su2bond(1,1,angm,bdim);

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
  tcgc[0].make_cgc(physpn,2,physpn);
  delete []angm;
  delete []bdim;
  delete []tele;
  delete []bb;
}

//--------------------------------------------------------------------------------------
void tensor_su2::makeup_input_vector(){
//--------------------------------------------------------------------------------------
  int i,j,k,nmom0,nmom1,a0,a1,bdim[2];
  nmom0=cgc.get_nmoment(0);
  nmom1=cgc.get_nmoment(1);
  for(i=0;i<nmom0;i++){
    a0=cgc.get_angularmoment(0,i);
    bdim[0]=cgc.get_bonddim(0,i);
    for(j=0;j<nmom1;j++){
      a1=cgc.get_angularmoment(1,j);
      bdim[1]=cgc.get_bonddim(1,j);
      if(a0==a1&&parr[i+j*nmom0]->is_null()){
	parr[i+j*nmom0]->alloc_space(2,bdim);
	*(pcgc[i+j*nmom0])=cgc_coef_singlet[a0];
      }
    }
  }
}

//--------------------------------------------------------------------------------------
void tensor_su2::initialize_input_vector(){
//--------------------------------------------------------------------------------------
  int i,j,k,nmom0,nmom1,a0,a1,bdim[2];
  su2bond bb[2];
  nmom0=cgc.get_nmoment(0);
  nmom1=cgc.get_nmoment(1);
  for(i=0;i<nmom0;i++){
    a0=cgc.get_angularmoment(0,i);
    bdim[0]=cgc.get_bonddim(0,i);
    for(j=0;j<nmom1;j++){
      a1=cgc.get_angularmoment(1,j);
      bdim[1]=cgc.get_bonddim(1,j);
      if(a0==a1){
	(*parr[i+j*nmom0])=1.;
	pcgc[i+j*nmom0]->make_singlet(a0);
      }
    }
  }
}

//--------------------------------------------------------------------------------------
void tensor_su2::direct_sum(int ind, tensor_su2& t1, tensor_su2& t2){
//--------------------------------------------------------------------------------------
  su2struct cgc1,cgc2;
  int i,j,k,l,*angm,*bdim,*cdim,nten1,nten2;  
  bool check;
  clean();
  cgc1=t1.get_cgc();
  cgc2=t2.get_cgc();
  cgc.direct_sum(ind,cgc1,cgc2);
  nbond=cgc.get_nbond();
  nten=cgc.get_nten();
  locspin=cgc.get_locspin();
  parr=new tensor*[nten];
  pcgc=new tensor*[nten];
  tarr=new tensor[nten];
  tcgc=new tensor[nten];
  for(i=0;i<nten;i++){
    parr[i]=&(tarr[i]);
    pcgc[i]=&(tcgc[i]);
  }
  angm=new int[nbond];
  bdim=new int[nbond];
  cdim=new int[nbond];
  nten1=t1.get_nten();
  nten2=t2.get_nten();
  for(i=0;i<nten;i++){
    check=get_tensor_argument(i,angm,bdim,cdim);
    if(check==false)continue;
    j=t1.get_tensor_index(angm);
    k=t2.get_tensor_index(angm);
    if(j==nten1&&k!=nten2){
      tarr[i]=*(t2.get_parr(k));
      tcgc[i]=*(t2.get_pcgc(k));
    }
    else if(j!=nten1&&k==nten2){
      tarr[i]=*(t1.get_parr(j));
      tcgc[i]=*(t1.get_pcgc(j));
    }
    else if(j!=nten1&&k!=nten2){
      tarr[i].direct_sum(ind,*(t1.get_parr(j)),*(t2.get_parr(k)));
      tcgc[i]=*(t1.get_pcgc(j));
    }
    else if(j==nten1&&k==nten2){
      cout<<"tensor_su2::direct_sum wrong with cgc"<<endl;
      exit(0);
    }
  }
  delete []angm;
  delete []bdim;
  delete []cdim;
}

//--------------------------------------------------------------------------------------
void tensor_su2::overlap_initial(tensor_su2& uu, tensor_su2& vv,int flag){
//--------------------------------------------------------------------------------------
  int a0,a1,a2,b0,b1,b2,c0,c1,c2,nmoma0,nmoma1,nmoma2,nmomb0,nmomb1,nmomb2,nmomc0,nmomc1,nmomc2,i0,i1,i2,j0,j1,j2,k0,k1,k2,m,n,p,i;
  su2bond *bb;
  tensor tmp1,tmp2,tmp3;
  int bdim[5],angm[3];
  double fac;
  clean();
  nbond=2;
  locspin=0;
  bb=new su2bond[nbond];
  uu.get_su2bond(2,bb[0]);
  vv.get_su2bond(2,bb[1]);
  bb[1].invert_bonddir();
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

	      if(a0==b0&&a2==b2&&a1==b1){
		m=i0+i1*nmoma0+i2*nmoma0*nmoma1;
		n=j0+j1*nmomb0+j2*nmomb0*nmomb1;
		i=i2+j2*nmoma2;
		if(uu.is_null(m)||vv.is_null(n))continue;
		if(tarr[i].is_null()){
		  tarr[i].contract_dmrg_overlap_initial(*(uu.get_parr(m)),*(vv.get_parr(n)),flag);
		  tcgc[i].make_identity(a2);
		}
		else{
		  tmp1.contract_dmrg_overlap_initial(*(uu.get_parr(m)),*(vv.get_parr(n)),flag);
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

//--------------------------------------------------------------------------------------
void tensor_su2::overlap_transformation(tensor_su2& uu, tensor_su2& vv, tensor_su2& op,int flag){
//--------------------------------------------------------------------------------------
  tensor_su2 tmp;
  clean();
  if(op.get_nbond()==0)return;
  if(flag==0)
    tmp.contract(op,0,uu,0);
  else if(flag==1){
    tmp.contract(op,0,uu,1);
    tmp.shift(2,0);
  }
  this->overlap_initial(tmp,vv,flag);
}

//--------------------------------------------------------------------------------------
void tensor_su2::operator_initial(tensor_su2& uu, tensor_su2& vv, tensor_su2& op, int flag){
//--------------------------------------------------------------------------------------
  //build three index operator from current site
  int a0,a1,a2,b0,b1,b2,c0,c1,c2,nmoma0,nmoma1,nmoma2,nmomb0,nmomb1,nmomb2,nmomc0,nmomc1,nmomc2,i0,i1,i2,j0,j1,j2,k0,k1,k2,m,n,p,i;
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
		    if(flag==0&&(a0!=b0||c0!=a1||c2!=b1)||flag==1&&(a1!=b1||c0!=a0||c2!=b0))continue;
		    m=i0+i1*nmoma0+i2*nmoma0*nmoma1;
		    n=j0+j1*nmomb0+j2*nmomb0*nmomb1;
		    p=k0+k1*nmomc0+k2*nmomc0*nmomc1;
		    i=i2+k1*nmoma2+j2*nmoma2*nmomc1;
		    if(uu.is_null(m)||vv.is_null(n)||op.is_null(p))continue;
		    tmp1.contract_dmrg_operator_initial(*(uu.get_parr(m)),*(vv.get_parr(n)),*(op.get_parr(p)),flag);
		    if(flag==0)
		      fac=fac_operator_onsite_left[a2][c1][b2][a0];
		    else if(flag==1)
		      fac=fac_operator_onsite_rght[a2][c1][b2][a1];
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

//--------------------------------------------------------------------------------------
void tensor_su2::operator_transformation(tensor_su2& uu, tensor_su2& vv, tensor_su2& op,int flag){
//--------------------------------------------------------------------------------------
  int a0,a1,a2,b0,b1,b2,c0,c1,c2,nmoma0,nmoma1,nmoma2,nmomb0,nmomb1,nmomb2,nmomc0,nmomc1,nmomc2,i0,i1,i2,j0,j1,j2,k0,k1,k2,m,n,p,i;
  su2bond *bb;
  tensor tmp1,tmp2,tmp3;
  int bdim[5],angm[3];
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
		    if(flag==0&&(c0!=a0||c2!=b0||a1!=b1))continue;
		    if(flag==1&&(c0!=a1||c2!=b1||a0!=b0))continue;
		    m=i0+i1*nmoma0+i2*nmoma0*nmoma1;
		    n=j0+j1*nmomb0+j2*nmomb0*nmomb1;
		    p=k0+k1*nmomc0+k2*nmomc0*nmomc1;
		    i=i2+k1*nmoma2+j2*nmoma2*nmomc1;
		    if(uu.is_null(m)||vv.is_null(n)||op.is_null(p))continue;
		    tmp1.contract_dmrg_operator_transformation(*(uu.get_parr(m)),*(vv.get_parr(n)),*(op.get_parr(p)),flag);
		    if(flag==0)
		      fac=fac_operator_transformation_left[a2][c1][b2][a0][b0];
		    else if(flag==1)
		      fac=fac_operator_transformation_rght[a2][c1][b2][a1][b1];
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

//--------------------------------------------------------------------------------------
void tensor_su2::operator_pairup(tensor_su2& uu, tensor_su2& vv, tensor_su2& op1, tensor_su2& op2, int flag){
//--------------------------------------------------------------------------------------
  int a0,a1,a2,b0,b1,b2,c0,c1,c2,d0,d1,d2,nmoma0,nmoma1,nmoma2,nmomb0,nmomb1,nmomb2,nmomc0,nmomc1,nmomc2,nmomd0,nmomd1,nmomd2,i0,i1,i2,j0,j1,j2,k0,k1,k2,l0,l1,l2,m,n,p,q,i;
  su2bond *bb;
  tensor tmp1,tmp2,tmp3;
  int bdim[5],angm[3];
  double fac;
  clean();
  if(op1.get_nbond()==0||op2.get_nbond()==0)return;
  nbond=2;
  locspin=0;
  bb=new su2bond[nbond];
  uu.get_su2bond(2,bb[0]);
  vv.get_su2bond(2,bb[1]);
  bb[1].invert_bonddir();
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
	      if(a2!=b2)continue;

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

			  if(flag==0&&(c0!=a0||c2!=b0||d0!=a1||d2!=b1||c1!=d1))continue;
			  if(flag==1&&(c0!=a1||c2!=b1||d0!=a0||d2!=b0||c1!=d1))continue;

			  m=i0+i1*nmoma0+i2*nmoma0*nmoma1;
			  n=j0+j1*nmomb0+j2*nmomb0*nmomb1;
			  p=k0+k1*nmomc0+k2*nmomc0*nmomc1;
			  q=l0+l1*nmomd0+l2*nmomd0*nmomd1;
			  i=i2+j2*nmoma2;
			  if(uu.is_null(m)||vv.is_null(n)||op1.is_null(p)||op2.is_null(q))continue;
			  tmp1.contract_dmrg_operator_pairup(*(uu.get_parr(m)),*(vv.get_parr(n)),*(op1.get_parr(p)),*(op2.get_parr(q)),flag);
			  if(flag==0)
			    fac=fac_operator_pairup_left[a2][a0][b0][c1];
			  else if(flag==1)
			    fac=fac_operator_pairup_rght[a2][a1][b1][c1];
			  tmp1*=fac;
			  if(tarr[i].is_null()){
			    tarr[i]=tmp1;
			    tcgc[i].make_identity(a2);
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

//--------------------------------------------------------------------------------------
void tensor_su2::hamiltonian_vector_multiplication(tensor_su2& vec, tensor_su2& op1, tensor_su2& op2){
//--------------------------------------------------------------------------------------
  int a0,a1,a2,b0,b1,b2,c0,c1,c2,nmoma0,nmoma1,nmoma2,nmomb0,nmomb1,nmomb2,nmomc0,nmomc1,nmomc2,i0,i1,i2,j0,j1,j2,k0,k1,k2,m,n,p,i,j,angm[3];
  su2bond *bb;
  tensor tmp,tmp1;
  tensor_su2 tmpsu2;
  double fac;
  clean();
  if(op1.get_nbond()==0||op2.get_nbond()==0)return;
  nbond=2;
  locspin=0;
  if(op1.get_nbond()==2&&op2.get_nbond()==2){
    tmpsu2.contract(op1,0,vec,0);
    this->contract(tmpsu2,1,op2,0);
    return;
  }
  bb=new su2bond[nbond];
  op1.get_su2bond(2,bb[0]);
  op2.get_su2bond(2,bb[1]);
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
  nmoma0=op1.get_nmoment(0);
  nmoma1=op1.get_nmoment(1);
  nmoma2=op1.get_nmoment(2);
  nmomb0=op2.get_nmoment(0);
  nmomb1=op2.get_nmoment(1);
  nmomb2=op2.get_nmoment(2);

  for(i0=0;i0<nmoma0;i0++){
    a0=op1.get_angularmoment(0,i0);
    for(i1=0;i1<nmoma1;i1++){
      a1=op1.get_angularmoment(1,i1);
      for(i2=0;i2<nmoma2;i2++){
	a2=op1.get_angularmoment(2,i2);
	angm[0]=a0;
	angm[1]=a1;
	angm[2]=a2;
	if(op1.check_angularmoments(angm)==false)continue;

	for(j0=0;j0<nmomb0;j0++){
	  b0=op2.get_angularmoment(0,j0);
	  for(j1=0;j1<nmomb1;j1++){
	    b1=op2.get_angularmoment(1,j1);
	    for(j2=0;j2<nmomb2;j2++){
	      b2=op2.get_angularmoment(2,j2);
	      angm[0]=b0;
	      angm[1]=b1;
	      angm[2]=b2;
	      if(op2.check_angularmoments(angm)==false)continue;

	      if(a0==b0&&a2==b2&&a1==b1){
		m=i0+i1*nmoma0+i2*nmoma0*nmoma1;
		n=j0+j1*nmomb0+j2*nmomb0*nmomb1;
		i=i2+j2*nmoma2;
		p=i0+j0*nmoma0;
		if(vec.is_null(p)||op1.is_null(m)||op2.is_null(n))continue;
		tmp1.contract(*(vec.get_parr(p)),0,*(op1.get_parr(m)),0);
		if(tarr[i].is_null()){
		  tarr[i].contract_dmrg_overlap_initial(tmp1,*(op2.get_parr(n)),0);
		  tarr[i]*=fac_hamilt_vec[a0][a1][a2];
		  tcgc[i]=cgc_coef_singlet[a2];
		}
		else{
		  tmp.contract_dmrg_overlap_initial(tmp1,*(op2.get_parr(n)),0);
		  tmp*=fac_hamilt_vec[a0][a1][a2];
		  tarr[i]+=tmp;
		}
	      }
	    }
	  }
	}
      }
    }
  }
}
