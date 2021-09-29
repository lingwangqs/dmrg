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
extern "C"{
  double ran_();
  void isort_(int*,int*);
}
extern int max_dcut,bdry,myrank,psize,jobid,memory_flag,disk,comm_rank;
extern double delta,qdelta;
extern double t1,t2,preparetime,lanczostime,svdtime;
extern tensor ***spin_op,***cgc_coef_left,***cgc_coef_rght,*cgc_coef_singlet,*identity;
extern double **spin_op_trace ;

//------------------------------------------------------------------------------
dmrg_su2::dmrg_su2(int xx, int yy, int sec, int r, int r2, int ex){
//------------------------------------------------------------------------------
  int i,j,k,l,p,m,x1,y1,x2,y2,x3,y3,x4,y4,max_dcut_ww,flag,max_angm;
  lx=xx;
  ly=yy;
  ns=lx*ly;
  physpn=1;
  phdim=physpn+1;
  read=r;
  read2=r2;
  exci=ex;
  totspin=sec;
  if(read2>read)max_dcut_ww=read2;
  else max_dcut_ww=read;
  if(max_dcut_ww<max_dcut)max_dcut_ww=max_dcut;
  max_angm=50;
  nfree=1;
  j=1;
  for(i=0;i<ns;i++){
    j*=phdim;
    if(j<=r&&r!=0||r==0&&j<=max_dcut)nfree++;
    else if(j>r&&r!=0||r==0&&j>max_dcut)break;
  }
  
  // Initialize lookup table storage
  // 5D sparse tables
  fac_permutation_left = LookupTable_5(0.0, max_angm,max_angm,max_angm,max_angm,18) ;
	fac_permutation_rght = LookupTable_5(0.0, max_angm,max_angm,max_angm,max_angm,18) ;
  fac_operator_transformation_left = LookupTable_5(0.0, max_angm,max_angm,max_angm,max_angm,max_angm) ;
	fac_operator_transformation_rght = LookupTable_5(0.0, max_angm,max_angm,max_angm,max_angm,max_angm) ;
  // 4D sparse tables
  fac_operator_onsite_left = LookupTable_4(0.0, max_angm,max_angm,max_angm,max_angm) ;
  fac_operator_onsite_rght = LookupTable_4(0.0, max_angm,max_angm,max_angm,max_angm) ; 
  fac_operator_pairup_left = LookupTable_4(0.0, max_angm,max_angm,max_angm,max_angm) ; 
  fac_operator_pairup_rght = LookupTable_4(0.0, max_angm,max_angm,max_angm,max_angm) ; 
  // 3D Dense table
  // Allocate and initialize fac_hamilt_vec to 0.0 
  fac_hamilt_vec.resize(max_angm) ;
  for(i=0;i<max_angm;i++) {
      fac_hamilt_vec[i].resize(max_angm)  ;
      for(j=0;j<max_angm;j++) {
          fac_hamilt_vec[i][j].resize(max_angm, 0.0) ; 
      }
  }  

  
  nfree=2;
  //if(ly==12&&max_dcut==6000&&totspin==2)nfree=4*ly;
  ww=new double*[exci+1];
  for(i=0;i<exci+1;i++)
    ww[i]=new double[max_dcut_ww*phdim];
  wtmp=new double[max_dcut_ww*phdim];
  hh=new tensor_su2[ns];
  uu=new tensor_su2[ns];
  opr=new tensor_su2*[ns];
  plqopr=new tensor_su2*[ns];
  qfll=new int[ns];
  qfrr=new int[ns];
  fll=new int*[ns];
  frr=new int*[ns];
  hmap=new int*[ns];
  plqpos=new int*[ns];
  plqflg=new int*[ns];
  for(i=0;i<ns;i++){
    opr[i]=new tensor_su2[ns];
    plqopr[i]=new tensor_su2[ns*2];
    fll[i]=new int[ns];
    frr[i]=new int[ns];
    hmap[i]=new int[ns];
    plqpos[i]=new int[4];
    plqflg[i]=new int[ns];
  }
  max_exci=exci;
  orth=new tensor_su2*[max_exci];
  ovlp=new tensor_su2*[max_exci];
  overlapvec=new tensor_su2[max_exci];
  for(i=0;i<max_exci;i++){
    orth[i]=new tensor_su2[ns];//orthogonal wavefunction
    ovlp[i]=new tensor_su2[ns];//overlap with orth wavefunc
  }
  //build up bond map
  for(i=0;i<ns;i++)
    for(j=0;j<ns;j++)
      hmap[i][j]=0;

  for(m=0;m<ns;m++){
    x1=m/ly;
    y1=m%ly;
    x2=(x1+0)%lx;;
    y2=(y1+1)%ly;
    x3=(x1+1)%lx;
    y3=(y1+1)%ly;
    x4=(x1+1)%lx;
    y4=(y1+0)%ly;
    if(x1%2==0)
      i=x1*ly+y1;
    else if(x1%2==1)
      i=x1*ly+(ly-1-y1);
    if(x2%2==0)
      j=x2*ly+y2;
    else if(x2%2==1)
      j=x2*ly+(ly-1-y2);
    if(x3%2==0)
      k=x3*ly+y3;
    else if(x3%2==1)
      k=x3*ly+(ly-1-y3);
    if(x4%2==0)
      l=x4*ly+y4;
    else if(x4%2==1)
      l=x4*ly+(ly-1-y4);
    plqpos[m][0]=i;
    plqpos[m][1]=j;
    plqpos[m][2]=k;
    plqpos[m][3]=l;    
    p=4;
    isort_(&p,plqpos[m]);
    //if(comm_rank==0)cout<<m<<"\t"<<plqpos[m][0]<<"\t"<<plqpos[m][1]<<"\t"<<plqpos[m][2]<<"\t"<<plqpos[m][3]<<endl;
  }
  for(k=0;k<4;k++){
    for(i=0;i<ns;i++){
      x1=i/ly;
      if(x1%2==0)
	y1=i%ly;
      else if(x1%2==1)
	y1=ly-1-(i%ly);
      flag=1;
      if(bdry==3){//shastry-sutherland model
	if(k==0){
	  x2=x1+1;
	  y2=y1;
	}
	else if(k==1){
	  x2=x1;
	  y2=(y1+1)%ly;
	}
	else if(k==2&&x1%2==0&&y1%2==0){
	  x2=x1+1;
	  y2=(y1+1)%ly;
	}
	else if(k==3&&x1%2==1&&y1%2==0){
	  x2=x1+1;
	  y2=(y1+ly-1)%ly;
	}
	else flag=0;
      }
      if(flag&&x2>=0&&x2<lx&&y2>=0&&y2<ly){
	if(x2%2==0)
	  j=x2*ly+y2;
	else if(x2%2==1)
	  j=x2*ly+(ly-1-y2);
	if(k==1&&(x1==0&&(x1%2==y1%2)||x1==lx-1&&(x1%2!=y1%2))){
	  hmap[i][j]=3;
	  hmap[j][i]=3;
	  if(comm_rank==0)cout<<"hmap3\t"<<i<<"\t"<<j<<endl;
	}
	else if(k<2){
	  hmap[i][j]=2;
	  hmap[j][i]=2;
	  if(comm_rank==0)cout<<"hmap2\t"<<i<<"\t"<<j<<endl;
	}
	else if(k>=2){
	  if(fabs(qdelta-1)<1.e-5){
	    hmap[i][j]=0;
	    hmap[j][i]=0;
	  }
	  else{
	    hmap[i][j]=1;
	    hmap[j][i]=1;
	  }
	  if(comm_rank==0)cout<<"hmap1\t"<<i<<"\t"<<j<<endl;
	}
      }
    }
  }

  //build up plqflg map
  for(i=0;i<ns;i++){
    qfll[i]=-1;
    qfrr[i]=-1;
    for(j=0;j<ns;j++){
      if((j/ly)%2!=(j%ly)%2&&(j/ly!=lx-1)&&i>=plqpos[j][0]&&i<=plqpos[j][3])
	plqflg[i][j]=1;
      else plqflg[i][j]=0;
      if(plqflg[i][j]&&i==plqpos[j][1]){
	qfll[i]=j;
	//if(comm_rank==0&&qfll[i]!=-1)
	//cout<<i<<"\tqfll\t"<<j<<"\t"<<plqpos[j][0]<<"\t"<<plqpos[j][1]<<"\t"<<plqpos[j][2]<<"\t"<<plqpos[j][3]<<endl;
      }
      if(plqflg[i][j]&&i==plqpos[j][2]){
	qfrr[i]=j;
	//if(comm_rank==0&&qfrr[i]!=-1)
	//cout<<i<<"\tqfrr\t"<<j<<"\t"<<plqpos[j][0]<<"\t"<<plqpos[j][1]<<"\t"<<plqpos[j][2]<<"\t"<<plqpos[j][3]<<endl;
      }
    }
  }
  //build up flagmap: at each i position, what is needed for later calculation
  for(i=0;i<ns;i++){
    for(j=0;j<ns;j++){
      fll[i][j]=0;
      frr[i][j]=0;
    }
    for(j=0;j<i;j++)//if there is bond acrossing i, all left j operator will be stored
      for(k=i+1;k<ns;k++)
	if(hmap[j][k]!=0)
	  fll[i][j]=1;
    for(j=0;j<ns;j++)//if there is plquette j such that for site i
      if(plqflg[i][j]&&i>plqpos[j][0]&&i<plqpos[j][1])
	fll[i][plqpos[j][0]]=1;
    for(j=i+1;j<ns;j++)//if there is bond acrossing i, all rght j operator will be stored
      for(k=0;k<i;k++)
	if(hmap[j][k]!=0)
	  frr[i][j]=1;    
    for(j=0;j<ns;j++)//if there is plquette j such that for site i
      if(plqflg[i][j]&&i>plqpos[j][2]&&i<plqpos[j][3])
	frr[i][plqpos[j][3]]=1;
  }
  makeup_clebsch_gordan_coefficient_tensors(max_angm);
  initialize_local_operators();
  //test_clebsch_gordan_coefficient();
}

//------------------------------------------------------------------------------
dmrg_su2::~dmrg_su2(){
//------------------------------------------------------------------------------
  int i,j,k,l;
  for(i=0;i<ns;i++){
    delete []opr[i];
    delete []plqopr[i];
    delete []fll[i];
    delete []frr[i];
    delete []hmap[i];
    delete []plqpos[i];
    delete []plqflg[i];
  }
  delete []uu;
  delete []hh;
  delete []opr;
  delete []plqopr;
  for(i=0;i<max_exci;i++){
    delete []ovlp[i];
    delete []orth[i];
  }
  delete []ovlp;
  delete []orth;
  delete []overlapvec;
  delete []qfll;
  delete []qfrr;
  delete []fll;
  delete []frr;
  delete []hmap;
  delete []plqpos;
  delete []plqflg;
  for(i=0;i<exci+1;i++)
    delete []ww[i];
  delete []ww;
  delete []wtmp;
  //delete global variables
  int max_angm=50;
  for(i=0;i<max_angm;i++){
    for(j=0;j<max_angm;j++){
      delete []cgc_coef_left[i][j];
      delete []cgc_coef_rght[i][j];
    }
    delete []cgc_coef_left[i];
    delete []cgc_coef_rght[i];
  }
  delete []cgc_coef_left;
  delete []cgc_coef_rght;
  delete []cgc_coef_singlet;
  delete []identity;
}

//------------------------------------------------------------------------------
void dmrg_su2::initialize_local_operators(){
//------------------------------------------------------------------------------
  tensor_su2 tmp1,tmp2,tmp3,tmp;
  int i;
  ring[0].make_spinhalf_permutation_rightmove(ring[2]);
  ring[1].make_spinhalf_permutation_leftmove();
  ring[3]=ring[0];
  ring[4]=ring[1];
  ring[5]=ring[2];
  ring[3].conjugate(1);
  ring[3].conjugate(2);
  ring[3].shift(2,0);
  ring[3].make_standard_cgc();

  ring[4].conjugate(1);
  ring[4].conjugate(2);
  ring[4].shift(2,0);
  ring[4].make_standard_cgc();
  ring[5].conjugate(1);
  ring[5].conjugate(2);
  ring[5].shift(2,0);
  ring[5].make_standard_cgc();

  ring[6].make_spinhalf_permutation_leftend();
  ring[7].make_spinhalf_permutation_rightend();
  tmp1.make_spinhalf_Qterm(qterm[0],qterm[1],qterm[2],qterm[3]);
  /*
  if(comm_rank==0){
    qterm[0].print();
    qterm[1].print();
    qterm[2].print();
    qterm[3].print();
  }
  */
  for(i=0;i<4;i++){
    qterm[4+i]=qterm[i];
    qterm[4+i].conjugate(1);
    qterm[4+i].conjugate(2);
    qterm[4+i].shift(2,0);
    qterm[4+i].make_standard_cgc();
  }
  sigma[0].make_spinor_start(physpn);
  sigma[1].make_spinor_start(physpn);
  if(physpn==2){
    sigma[0]*=sqrt(2);
    sigma[1]*=-sqrt(2);
    tmp1=sigma[0];
    tmp2.contract(tmp1,2,tmp1,0);
    tmp2.fuse(1,2);
    tmp2/=sqrt(3.);
    sigma[0].direct_sum(1,tmp2,tmp1);
    tmp1*=-1.;
    sigma[1].direct_sum(1,tmp2,tmp1);
    sigma[0]/=sqrt(2.);
    sigma[1]/=sqrt(2.);
  }
  else if(physpn==1){
    sigma[0]*=sqrt(3)/2;
    sigma[1]*=-sqrt(3)/2;
  }
  else if(physpn==3){
    sigma[0]/=2.*0.25819888974716109775;
    sigma[1]/=-2.*0.25819888974716109775;
    tmp1=sigma[0];
    tmp2.contract(tmp1,2,tmp1,0);
    tmp2.fuse(1,2);
    tmp3.contract(tmp2,2,tmp1,0);
    tmp3.fuse(1,2);
    tmp3*=sqrt(16./243.);
    tmp2*=sqrt(116./243.);
    tmp.direct_sum(1,tmp3,tmp1);
    sigma[0].direct_sum(1,tmp,tmp2);
    tmp*=-1;
    sigma[1].direct_sum(1,tmp,tmp2);
    sigma[0]*=sqrt(27./160.);
    sigma[1]*=sqrt(27./160.);
  }
}

//------------------------------------------------------------------------------
void dmrg_su2::lanczos_solve_eigenvector_idmrg(int il, int ir,tensor_su2& vec){
//------------------------------------------------------------------------------
  lanczos_su2 lan;
  int mlanc=50,nele1,nele2;
  vec.get_nelement(nele1,nele2);
  if(mlanc>nele1)mlanc=nele1+1;
  cout<<"mlanc="<<mlanc<<endl;
  lan.initialize_lanczos(vec,mlanc);
  lan.lanczos1(il,ir,vec,minlan);
  gs_enr[exci]=lan.get_eigval();
  if(comm_rank==0)cout<<il<<" gs_enr["<<exci<<"]="<<gs_enr[exci]<<endl;
}

//------------------------------------------------------------------------------
void dmrg_su2::hamiltonian_vector_multiplication_idmrg(int il, int ir, tensor_su2& vec1, tensor_su2& vec2){
//------------------------------------------------------------------------------

  tensor_su2 tmp,htmp,*vec2arr,hmlt;
  int i,j,k,l,m,m0,m1,m2,m3,m4,m5,m6,m7,*count;
  double fac,overlap;
  double t3,t4,*timer;
  static int cnt=0;
  count=new int[psize];
  timer=new double[psize];
  vec2arr=new tensor_su2[psize];
  for(i=0;i<psize;i++){
    timer[i]=0;
    count[i]=0;
  }
  if(hh[il].get_nbond()==0&&hh[ir].get_nbond()==0){
    vec2=vec1;
    vec2=0;
  }
  else if(hh[il].get_nbond()!=0&&hh[ir].get_nbond()!=0){
    vec2.contract(hh[il],0,vec1,0);
    tmp.contract(vec1,1,hh[ir],0);
    vec2+=tmp;
  }
  else if(hh[il].get_nbond()==0&&hh[ir].get_nbond()!=0)
    vec2.contract(vec1,1,hh[ir],0);
  else if(hh[il].get_nbond()!=0&&hh[ir].get_nbond()==0)
    vec2.contract(hh[il],0,vec1,0);
  m1=(il+1)*(ns-ir);
  m2=m1+exci;
  m3=m2+ns*2;
#pragma omp parallel for default(shared) private(i,j,k,l,m,m0,overlap,fac,tmp,myrank,t3,t4) schedule(dynamic,1)
  for(m0=0;m0<m3;m0++){
    myrank=omp_get_thread_num();
    if(m0<m1){
      m=m0;
      j=m/(ns-ir);
      k=m%(ns-ir)+ir;
      if(hmap[j][k]!=0||j==il&&k==ir&&il+1!=ir){
	tmp.hamiltonian_vector_multiplication(vec1,opr[il][j],opr[ir][k], fac_hamilt_vec);
	if(hmap[j][k]==1)fac=1-qdelta;
	else if(hmap[j][k]==2)fac=delta+qdelta/4.;
	else if(hmap[j][k]==3)fac=delta;
	else if(j==il&&k==ir&&il+1!=ir)fac=1;
	tmp*=fac;
	if(tmp.get_nbond()==0){
	  cout<<"something wrong with idmrg"<<endl;
	  exit(0);
	}
	if(vec2arr[myrank].get_nbond()==0)
	  vec2arr[myrank]=tmp;
	else vec2arr[myrank]+=tmp;
	count[myrank]++;
      }
    }
    else if(m0>=m1&&m0<m2){
      if(il==ir-1){
	j=m0-m1;
	overlap=vec1.inner_prod(overlapvec[j]);
	overlap*=gs_enr[j];
	tmp=overlapvec[j];
	tmp*=-overlap;
	if(vec2arr[myrank].get_nbond()==0)
	  vec2arr[myrank]=tmp;
	else vec2arr[myrank]+=tmp;
	count[myrank]++;
      }
    }
    else if(m0>=m2&&m0<m3){
      m=m0-m2;
      j=(m0-m2)%ns;
      k=(m0-m2)/ns;
      if(il!=ir-1||plqflg[il][j]!=1||plqflg[ir][j]!=1)continue;
      if(il>=plqpos[j][0]&&ir<=plqpos[j][1])
	tmp.hamiltonian_vector_multiplication(vec1,opr[il][plqpos[j][0]],plqopr[ir][m], fac_hamilt_vec);
      else if(il>=plqpos[j][1]&&ir<=plqpos[j][2])
	tmp.hamiltonian_vector_multiplication(vec1,plqopr[il][m],plqopr[ir][m], fac_hamilt_vec);
      else if(il>=plqpos[j][2]&&ir<=plqpos[j][3])
	tmp.hamiltonian_vector_multiplication(vec1,plqopr[il][m],opr[ir][plqpos[j][3]], fac_hamilt_vec);
      if(tmp.get_nbond()==0){
	cout<<"plq something wrong with idmrg"<<endl;
	cout<<plqpos[j][0]<<"\t"<<plqpos[j][1]<<"\t"<<plqpos[j][2]<<"\t"<<plqpos[j][3]<<"\t"<<il<<"\t"<<ir<<endl;
	exit(0);
      }
      tmp*=-qdelta;
      if(vec2arr[myrank].get_nbond()==0)
	vec2arr[myrank]=tmp;
      else vec2arr[myrank]+=tmp;
      count[myrank]++;
    }
  }
  for(i=0;i<psize;i++){
    if(vec2arr[i].get_nbond()!=0)
      vec2+=vec2arr[i];
  }
  delete []vec2arr;
  delete []timer;
  delete []count;
}

//------------------------------------------------------------------------------
void dmrg_su2::test_clebsch_gordan_coefficient(){
//------------------------------------------------------------------------------
  int i,j,k,max_angm=10;
  tensor tmp,tmp1,tmp2,sum;
  double nor;
  for(i=1;i<=1;i++){
    for(j=0;j<=2;j+=2){
      for(k=1;k<=3;k+=2){
	if(k<=i+j&&k>=abs(i-j)){
	  tmp1.contract(cgc_coef_left[i][j][k],2,cgc_coef_singlet[k],1);    
	  tmp2.contract(tmp1,1,cgc_coef_singlet[j],0);
	  if(!tmp2.is_proportional_to(cgc_coef_left[k][i][j],nor)){
	    tmp2.print();
	    cout<<"can not convert direction i="<<i<<"\tj="<<j<<"\tk="<<k<<endl;
	  }
	  else{
	    //tmp2.print();
	    //cout<<"can convert direction i="<<i<<"\tj="<<j<<"\tk="<<k<<endl;
	  }
	}
      }
    }
  }
  for(k=0;k<max_angm;k+=1){
    tmp.contract(cgc_coef_singlet[k],0,cgc_coef_singlet[k],0);
    tmp*=k+1;
    if(comm_rank==0){
      cout<<"k="<<k<<endl;
      tmp.print();
    }
  }
  return;
  i=2;
  j=2;
  for(k=0;k<=4;k+=2){
    tmp1.contract(cgc_coef_left[i][j][k],0,cgc_coef_singlet[i],1);
    tmp2.contract(tmp1,0,cgc_coef_singlet[j],1);
    tmp1.contract(tmp2,0,cgc_coef_singlet[k],0);
    if(tmp1.is_proportional_to(cgc_coef_left[i][j][k],nor)){
      if(comm_rank==0){
	cout<<"i="<<i<<"\tj="<<j<<"\tk="<<k<<"\tnor="<<nor*sqrt((i+1)*(j+1)*(k+1))<<endl;
      }
    }
  }
}

//------------------------------------------------------------------------------
void dmrg_su2::makeup_clebsch_gordan_coefficient_tensors(int max_angm){
//------------------------------------------------------------------------------

  int i,j,k,l,m,a,a0,a1,a2,a3;
  tensor tmp,tmp1,tmp2;
  double nor,sgn;
  bool check1,check2;
  cout<<"max_angm="<<max_angm<<endl;
  cgc_coef_left=new tensor**[max_angm];
  cgc_coef_rght=new tensor**[max_angm];
  for(i=0;i<max_angm;i++){
    cgc_coef_left[i]=new tensor*[max_angm];
    cgc_coef_rght[i]=new tensor*[max_angm];
    for(j=0;j<max_angm;j++){
      cgc_coef_left[i][j]=new tensor[max_angm];
      cgc_coef_rght[i][j]=new tensor[max_angm];
    }
  }
  cgc_coef_singlet=new tensor[max_angm];
  identity=new tensor[max_angm];
  for(i=0;i<max_angm;i++){
    cgc_coef_singlet[i].make_singlet(i);
    identity[i].make_identity(i);
    //if(comm_rank==0){tmp.contract(cgc_coef_singlet[i],1,cgc_coef_singlet[i],0);tmp1=tmp*(i+1);tmp1.print();}
  }
        
  for(i=0;i<=6;i++)
    for(j=0;j<max_angm;j++)
      for(k=0;k<max_angm;k++)
	if(k>=abs(i-j)&&k<=abs(i+j))
	  cgc_coef_rght[i][j][k].make_cgc(i,j,k);
  for(i=0;i<max_angm;i++)
    for(j=0;j<=6;j++)
      for(k=0;k<max_angm;k++)
	if(k>=abs(i-j)&&k<=abs(i+j)){
	  cgc_coef_left[i][j][k].make_cgc(i,j,k);
	  tmp=cgc_coef_left[i][j][k];
	  tmp.shift(2,0);
	  //if(comm_rank==0&&j==2&&i==k&&(i==3||i==2||i==1))tmp.print();
	}
  for(i=0;i<max_angm;i++)
    for(j=0;j<=6;j+=2)
      for(k=0;k<max_angm;k++)
	if(k>=abs(i-j)&&k<=abs(i+j)&&(i%2==k%2)){
	  tmp1.contract(cgc_coef_left[i][j][k],0,cgc_coef_singlet[i],1);
	  tmp2.contract(tmp1,0,cgc_coef_singlet[j],1);
	  tmp1.contract(tmp2,0,cgc_coef_singlet[k],1);
	  fac_hamilt_vec[i][j][k]=cgc_coef_left[i][j][k].inner_prod(tmp1)*sqrt(j+1);
	  //fac_hamilt_vec[i][j][k]=cgc_coef_left[i][j][k].inner_prod(tmp1);
	  //if(comm_rank==0)cout<<"i="<<i<<"\tj="<<j<<"\tk="<<k<<"\t"<<fac_hamilt_vec[i][j][k]<<"\t"<<sqrt(k+1)/sqrt(i+1)/sqrt(j+1)<<endl;
	}

 
  for(m=0;m<=6;m+=2)
    for(i=0;i<max_angm;i++)
      for(j=0;j<max_angm;j++)
	if(j>=abs(i-physpn)&&j<=abs(i+physpn))
	  for(k=0;k<max_angm;k++)
	    for(l=0;l<max_angm;l++)
	      if(l>=abs(k-physpn)&&l<=abs(k+physpn)){
		if(i==k&&l>=abs(j-m)&&l<=abs(j+m)&&(j%2==l%2)){
		  if(physpn>=abs(physpn-m)&&physpn<=abs(physpn+m)){
		    //left operator initialize
		    tmp1.contract_dmrg_operator_initial(cgc_coef_left[i][physpn][j],cgc_coef_left[k][physpn][l],cgc_coef_left[physpn][m][physpn],0);
		    check1=tmp1.is_proportional_to(cgc_coef_left[j][m][l],fac_operator_onsite_left(j,m,l,i));
		    if(check1==false&&comm_rank==0){
		      //cout<<"wrong operator_initialize_left\ti="<<i<<"\tj="<<j<<"\tk="<<k<<"\tl="<<l<<"\tm="<<m<<"\t"<<fac_operator_onsite_left[j][m][l][i]<<endl;
		    }
		    //rght operator initialize
		    tmp2.contract_dmrg_operator_initial(cgc_coef_rght[physpn][i][j],cgc_coef_rght[physpn][k][l],cgc_coef_left[physpn][m][physpn],1);
		    check2=tmp2.is_proportional_to(cgc_coef_left[j][m][l],fac_operator_onsite_rght(j,m,l,i));
		    if(check2==false&&comm_rank==0){
		      //cout<<"wrong operator_initialize_rght\ti="<<i<<"\tj="<<j<<"\tk="<<k<<"\tl="<<l<<"\tm="<<m<<"\t"<<fac_operator_onsite_rght[j][m][l][i]<<endl;
		    }
		    //else if(check2==true&&comm_rank==0)cout<<"crrct operator_initialize_rght\ti="<<i<<"\tj="<<j<<"\tk="<<k<<"\tl="<<l<<"\t"<<fac_operator_onsite_rght[j][m][l][i]<<endl;
		  }	    
		}  
	      }

  ////////////////////////////////////
  for(m=0;m<=6;m+=2)
    for(i=0;i<max_angm;i++)
      for(j=0;j<max_angm;j++)
	if(j>=abs(i-physpn)&&j<=abs(i+physpn))
	  for(k=0;k<max_angm;k++)
	    for(l=0;l<max_angm;l++)
	      if(l>=abs(k-physpn)&&l<=abs(k+physpn)){
		if(k>=abs(i-m)&&k<=abs(i+m)&&(k%2==i%2)&&l>=abs(j-m)&&l<=(j+m)&&(j%2==l%2)){
		  //left canonical form
		  tmp1.contract_dmrg_operator_transformation(cgc_coef_left[i][physpn][j],cgc_coef_left[k][physpn][l],cgc_coef_left[i][m][k],0);
      fac_operator_transformation_left(j,m,l,i,k) = 0.0 ;
		  check1=tmp1.is_proportional_to(cgc_coef_left[j][m][l],fac_operator_transformation_left(j,m,l,i,k));
		  if(check1==false&&comm_rank==0){
		    //cout<<"wrong operator_transformation_left\ti="<<i<<"\tj="<<j<<"\tk="<<k<<"\tl="<<l<<"\tm="<<m<<"\t"<<fac_operator_transformation_left[j][m][l][i][k]<<endl;
		  }
		  //rght canonical form
		  tmp2.contract_dmrg_operator_transformation(cgc_coef_rght[physpn][i][j],cgc_coef_rght[physpn][k][l],cgc_coef_left[i][m][k],1);
      fac_operator_transformation_rght(j,m,l,i,k) = 0.0 ;
		  check2=tmp2.is_proportional_to(cgc_coef_left[j][m][l],fac_operator_transformation_rght(j,m,l,i,k));
		  if(check2==false&&comm_rank==0){
		    //cout<<"wrong operator_transformation_rght\ti="<<i<<"\tj="<<j<<"\tk="<<k<<"\tl="<<l<<"\tm="<<m<<"\t"<<fac_operator_transformation_rght[j][m][l][i][k]<<endl;
		  }
		  //else if(check2==true&&comm_rank==0)cout<<"crrct operator_transformation_rght\ti="<<i<<"\tj="<<j<<"\tk="<<k<<"\tl="<<l<<"\tm="<<m<<"\t"<<fac_operator_transformation_rght[j][m][l][i][k]<<endl;
		}	
	      }
 
  for(m=0;m<=6;m+=2)
    for(i=0;i<max_angm;i++)
      for(j=0;j<max_angm;j++)
	if(j>=abs(i-physpn)&&j<=abs(i+physpn))
	  for(k=0;k<max_angm;k++)
	    for(l=0;l<max_angm;l++)
	      if(l>=abs(k-physpn)&&l<=abs(k+physpn)){
		if(k>=abs(i-m)&&k<=(i+m)&&j==l&&physpn>=abs(physpn-m)&&physpn<=abs(physpn+m)){
		  tmp.contract(cgc_coef_left[physpn][m][physpn],1,cgc_coef_singlet[m],1);
		  tmp.shift(0,2);
		  tmp*=sqrt(m+1);
		  tmp1.contract_dmrg_operator_pairup(cgc_coef_left[i][physpn][j],cgc_coef_left[k][physpn][l],cgc_coef_left[i][m][k],tmp,0);
		  check1=tmp1.is_proportional_to(identity[j],fac_operator_pairup_left(j,i,k,m));
		  if(check1==false&&comm_rank==0){
		    //cout<<"wrong operator_pairup_left\ti="<<i<<"\tj="<<j<<"\tk="<<k<<"\tl="<<l<<"\t"<<fac_operator_pairup_left[j][i][k][m]<<endl;
		  }
		  //else if(check1==true&&comm_rank==0)cout<<"crrct operator_pairup_rght\ti="<<i<<"\tj="<<j<<"\tk="<<k<<"\tl="<<l<<"\t"<<fac_operator_pairup_left[j][i][k][m]<<endl;
		  tmp2.contract_dmrg_operator_pairup(cgc_coef_rght[physpn][i][j],cgc_coef_rght[physpn][k][l],cgc_coef_left[i][m][k],tmp,1);
		  check2=tmp2.is_proportional_to(identity[j],fac_operator_pairup_rght(j,i,k,m));
		  if(check2==false&&comm_rank==0){
		    //cout<<"wrong operator_pairup_rght\ti="<<i<<"\tj="<<j<<"\tk="<<k<<"\tl="<<l<<"\t"<<fac_operator_pairup_rght[j][i][k][m]<<endl;
		  }
		  //else if(check2==true&&comm_rank==0)cout<<"crrct operator_pairup_rght\ti="<<i<<"\tj="<<j<<"\tk="<<k<<"\tl="<<l<<"\t"<<fac_operator_pairup_rght[j][i][k][m]<<endl;
		}
	      }
  ////////////////////////////////////
  cout<<"start permutation"<<endl;
    for(i=0;i<max_angm;i++)
      for(j=0;j<max_angm;j++)
	if(j>=abs(i-physpn)&&j<=abs(i+physpn)){
	  for(k=0;k<max_angm;k++)
	    for(l=0;l<max_angm;l++)
	      if(l>=abs(k-physpn)&&l<=abs(k+physpn)){
		for(m=0;m<18;m++){
		  a0=m;
		  a1=(a0%3)*2;a0/=3;
		  a2=(a0%2)*2+1;a0/=2;
		  a3=(a0%3)*2;
		  if(!(a2<=abs(a1+physpn)&&a2>=abs(a1-physpn)&&a3<=abs(a2+physpn)&&a3>=abs(a2-physpn)))continue;
		  if((i%2==k%2)&&(j%2==l%2)&&k<=abs(a1+i)&&k>=abs(a1-i)&&l<=abs(a3+j)&&l>=abs(a3-j)){
		    tmp1.contract_dmrg_permutation(cgc_coef_left[i][physpn][j],cgc_coef_left[k][physpn][l],cgc_coef_left[i][a1][k],cgc_coef_left[physpn][a1][a2],cgc_coef_left[physpn][a3][a2],0);
		    if(!tmp1.is_proportional_to(cgc_coef_left[j][a3][l],nor)){
		      if(comm_rank==0){
			//cout<<"contract permutation_left wrong\t"<<endl;
			tmp1.print();
			cgc_coef_left[j][a3][l].print();
		      }
		      exit(0);
		    }
		    else{
		      fac_permutation_left(i,j,k,l,m)=nor;
		      //if(comm_rank==0) cout<<"contract permutation_left right\t"<<endl;
		    }
		  }
		  if((i%2==k%2)&&(j%2==l%2)&&k<=abs(a3+i)&&k>=abs(a3-i)&&l<=abs(a1+j)&&l>=abs(a1-j)){
		    tmp2.contract_dmrg_permutation(cgc_coef_rght[physpn][i][j],cgc_coef_rght[physpn][k][l],cgc_coef_left[i][a3][k],cgc_coef_left[a2][physpn][a1],cgc_coef_left[a2][physpn][a3],1);
		    if(!tmp2.is_proportional_to(cgc_coef_left[j][a1][l],nor)){
		      if(comm_rank==0){
			//cout<<"contract permutation_rght wrong\t"<<endl;
			tmp1.print();
			cgc_coef_left[j][a1][l].print();
		      }
		      exit(0);
		    }
		    else{
		      fac_permutation_rght(i,j,k,l,m)=nor;
		      //if(comm_rank==0) cout<<"contract permutation_rght right\t"<<endl;
		    }
		  }
		}
	      }
	}
    return;
}
