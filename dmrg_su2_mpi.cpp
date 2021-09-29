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
}
extern int max_dcut,bdry,psize,jobid,memory_flag,comm_rank,comm_size,myrank;
extern double delta,qdelta;
extern double t1,t2,preparetime,lanczostime,svdtime;

//------------------------------------------------------------------------------
void dmrg_su2::do_idmrg(){
//------------------------------------------------------------------------------
  int i,j,k,m,n;
  tensor_su2 tmp,htmp,tconj,vec;
  su2bond *bb;
  double fac;
  cout<<comm_rank<<" begin idmrg"<<endl;

  bb=new su2bond[4];
  if(comm_rank==1)
    for(j=0;j<exci;j++){
      orth[j][ns/2].multiply_singular_value(2,ww[j]);
    }
  for(i=0;i<ns/2;i++){
    for(j=0;j<exci;j++){
      if(comm_rank==0){
	send_tensor_su2(orth[j][i],1);
	recv_tensor_su2(orth[j][ns-1-i],1);
      }
      else if(comm_rank==1){
	recv_tensor_su2(orth[j][i],0);
	send_tensor_su2(orth[j][ns-1-i],0);
      }
    }
  }

  i=0;
  j=1;
  bb[0].set_su2bond(1,1,&i,&j);
  i=physpn;
  j=1;
  bb[1].set_su2bond(1,1,&i,&j);
  uu[0].fuse(bb[0],bb[1]);
  i=physpn;
  j=1;
  bb[2].set_su2bond(1,1,&i,&j);
  i=totspin;
  j=1;
  bb[3].set_su2bond(1,1,&i,&j);
  uu[ns-1].fuse(bb[2],bb[3]);

  for(i=0;i<=ns/2-2;i++){
    if(comm_rank==0)	cout<<"rank0 i="<<i<<endl;
    prepare_site_operator_from_left(i,0);
    prepare_site_operator_from_right(ns-1-i,0);
    if(i==ns/2-2){
      for(j=0;j<exci;j++){
	ovlp[j][ns/2-1].contract(ovlp[j][ns/2-2],0,orth[j][ns/2-1],0);
	ovlp[j][ns/2].contract(orth[j][ns/2],1,ovlp[j][ns/2+1],0);
	ovlp[j][ns/2].shift(1,0);
      }
    }
    uu[i].get_su2bond(2,bb[0]);
    uu[ns-1-i].get_su2bond(2,bb[3]);
    bb[0].invert_bonddir();
    bb[3].invert_bonddir();
    uu[i+1].fuse(bb[0],bb[1]);
    uu[ns-2-i].fuse(bb[2],bb[3]);
    if((i+2)*2*physpn>=totspin&&i>=1||i==ns/2-2){
      if(comm_rank==0){
	prepare_input_vector_initial(i+1,ns-2-i,vec);
	prepare_site_operator_from_left(i+1);
	prepare_site_operator_from_right(ns-2-i);
	vec.initialize_input_vector();
	lanczos_solve_eigenvector_idmrg(i+1,ns-2-i,vec);
	vec.svd(uu[i+1],0,uu[ns-2-i],0,ww[exci]);
      }
      if(comm_rank==0){
	send_tensor_su2(uu[i+1],1);
	send_tensor_su2(uu[ns-2-i],1);
	MPI::COMM_WORLD.Send(ww[exci],max_dcut,MPI_DOUBLE,1,0);
      }
      else if(comm_rank==1){
	uu[i+1].clean();
	uu[ns-2-i].clean();
	recv_tensor_su2(uu[i+1],0);
	recv_tensor_su2(uu[ns-2-i],0);
	MPI::COMM_WORLD.Recv(ww[exci],max_dcut,MPI_DOUBLE,0,0);
      }
    }
  }

  delete []bb;
  MPI::COMM_WORLD.Barrier();
  meet_again_at_center(0);
  MPI::COMM_WORLD.Barrier();
  cout<<comm_rank<<" done idmrg"<<endl;
}

//------------------------------------------------------------------------------
void dmrg_su2::prepare_sweep(){
//------------------------------------------------------------------------------
  int i,j;
  tensor_su2 tconj;
  su2bond bb[4];
  {
    //reset boundary two spins to make SU(2)_v1 and current SU(2) consistent
    i=0;
    j=1;
    bb[0].set_su2bond(1,1,&i,&j);
    i=physpn;
    j=1;
    bb[1].set_su2bond(1,1,&i,&j);
    uu[0].fuse(bb[0],bb[1]);
    i=physpn;
    j=1;
    bb[2].set_su2bond(1,1,&i,&j);
    i=totspin;
    j=1;
    bb[3].set_su2bond(1,1,&i,&j);
    uu[ns-1].fuse(bb[2],bb[3]);
  }

  for(i=0;i<exci;i++){
    if(comm_rank==1){
      orth[i][ns/2].multiply_singular_value(2,ww[i]);
      send_tensor_su2(orth[i][ns/2],0);
      recv_tensor_su2(orth[i][ns/2-1],0);
    }
    else if(comm_rank==0){
      recv_tensor_su2(orth[i][ns/2],1);
      send_tensor_su2(orth[i][ns/2-1],1);
    }
  }
  for(i=0;i<ns/2-1;i++){
    if(comm_rank==0){
      cout<<"prepare sweep i="<<i<<endl;
      prepare_site_operator_from_left(i,0);
    }
    if(comm_rank==1){
      cout<<"prepare sweep i="<<ns-1-i<<endl;
      prepare_site_operator_from_right(ns-1-i,0);
    }
  }
  MPI::COMM_WORLD.Barrier();
  meet_again_at_center(0);
}

//------------------------------------------------------------------------------
void dmrg_su2::meet_again_at_center(int flag){
//------------------------------------------------------------------------------
  int i,j,k,src,des,tag;
  tensor_su2 vec,htmp,tconj;
  su2bond bb[4];

  synchronize_site(ns/2+1,1,1,0);
  MPI::COMM_WORLD.Barrier();
  if(comm_rank==1){
    des=0;
    send_tensor_su2(uu[ns/2],des);
    for(j=0;j<exci;j++){
      ovlp[j][ns/2].contract(orth[j][ns/2],1,ovlp[j][ns/2+1],0);
      ovlp[j][ns/2].shift(1,0);
      send_tensor_su2(ovlp[j][ns/2],des);
    }
  }
  else if(comm_rank==0){
    src=1;
    recv_tensor_su2(uu[ns/2],src);
    for(j=0;j<exci;j++){
      ovlp[j][ns/2-1].contract(ovlp[j][ns/2-2],0,orth[j][ns/2-1],0);
      recv_tensor_su2(ovlp[j][ns/2],src);
    }
  }
  MPI::COMM_WORLD.Barrier();
  if(comm_rank==0){
    if(flag==1)
      uu[ns/2].devide_singular_value(2,ww[exci]);
    else if(flag==0)//smoothly connect from non_mpi version
      uu[ns/2].multiply_singular_value(2,ww[exci]);
    prepare_input_vector(ns/2-1,ns/2,vec);
    prepare_site_operator_from_left(ns/2-1);
    prepare_site_operator_from_right(ns/2);
    lanczos_solve_eigenvector_idmrg(ns/2-1,ns/2,vec);
    vec.svd(uu[ns/2-1],0,uu[ns/2],0,ww[exci]);
  }
  MPI::COMM_WORLD.Barrier();
  if(comm_rank==0){
    des=1;
    send_tensor_su2(uu[ns/2],des);
    send_tensor_su2(uu[ns/2-1],des);
    MPI::COMM_WORLD.Send(ww[exci],max_dcut,MPI_DOUBLE,des,0);
    MPI::COMM_WORLD.Send(&minlan,1,MPI_DOUBLE,des,0);
  }
  else if(comm_rank==1){
    src=0;
    recv_tensor_su2(uu[ns/2],src);
    recv_tensor_su2(uu[ns/2-1],src);
    MPI::COMM_WORLD.Recv(ww[exci],max_dcut,MPI_DOUBLE,src,0);
    MPI::COMM_WORLD.Recv(&minlan,1,MPI_DOUBLE,src,0);
  }
  synchronize_site(ns/2-2,0,0,1);
  MPI::COMM_WORLD.Barrier();
  save_mps1();
  MPI::COMM_WORLD.Barrier();
  save_mps2();
  MPI::COMM_WORLD.Barrier();
  save_ww();
  save_enr();
  if(comm_rank==0)
    uu[ns/2-1].multiply_singular_value(2,ww[exci]);
  else if(comm_rank==1)
    uu[ns/2].multiply_singular_value(2,ww[exci]);
  MPI::COMM_WORLD.Barrier();
  if(comm_rank==0)cout<<"after meet"<<endl;
}

//------------------------------------------------------------------------------
void dmrg_su2::prepare_input_vector_initial(int il,int ir, tensor_su2& vec){
//------------------------------------------------------------------------------
  tensor_su2 tmp1,tmp2,left,rght,tmp,tmp3;
  double e,o;
  int k,j;
  su2bond bb[4];
  uu[il].get_su2bond(0,bb[0]);
  uu[il].get_su2bond(1,bb[1]);
  uu[ir].get_su2bond(0,bb[2]);
  uu[ir].get_su2bond(1,bb[3]);
  tmp1.fuse(bb[0],bb[1]);
  tmp2.fuse(bb[2],bb[3]);
  left.overlap_initial(uu[il],tmp1,0);
  rght.overlap_initial(uu[ir],tmp2,1);

  left.get_su2bond(0,bb[0]);
  rght.get_su2bond(0,bb[1]);
  bb[0].invert_bonddir();
  bb[1].invert_bonddir();
  tmp.fuse_to_singlet(bb[0],bb[1]);
  tmp3.contract(left,0,tmp,0);
  vec.contract(tmp3,1,rght,0);
  vec.makeup_input_vector();
  if(il+1!=ir)return;
  for(j=0;j<exci;j++){
    left.overlap_initial(ovlp[j][il],tmp1,0);
    left.take_conjugate(0);
    rght.overlap_initial(ovlp[j][ir],tmp2,1);
    overlapvec[j].contract(left,0,rght,0);
    overlapvec[j].makeup_input_vector();
  }
}

//------------------------------------------------------------------------------
void dmrg_su2::prepare_input_vector(int il,int ir, tensor_su2& vec){
//------------------------------------------------------------------------------
  tensor_su2 tmp1,tmp2,left,rght,tmp,tmp3;
  double e,o;
  int k,j;
  su2bond bb[4];
  uu[il].get_su2bond(0,bb[0]);
  uu[il].get_su2bond(1,bb[1]);
  uu[ir].get_su2bond(0,bb[2]);
  uu[ir].get_su2bond(1,bb[3]);
  tmp1.fuse(bb[0],bb[1]);
  tmp2.fuse(bb[2],bb[3]);
  left.overlap_initial(uu[il],tmp1,0);
  left.take_conjugate(0);
  rght.overlap_initial(uu[ir],tmp2,1);
  vec.contract(left,0,rght,0);
  vec.makeup_input_vector();
  if(il+1!=ir)return;
  for(j=0;j<exci;j++){
    left.overlap_initial(ovlp[j][il],tmp1,0);
    left.take_conjugate(0);
    rght.overlap_initial(ovlp[j][ir],tmp2,1);
    overlapvec[j].contract(left,0,rght,0);
    overlapvec[j].makeup_input_vector();
  }
}

//------------------------------------------------------------------------------
void dmrg_su2::sweep(){
//------------------------------------------------------------------------------
  if(comm_rank==0){
    sweep_to_left();
  }
  else if(comm_rank==1){
    sweep_to_right();
  }
  MPI::COMM_WORLD.Barrier();
  meet_again_at_center(1);
  MPI::COMM_WORLD.Barrier();
}

//------------------------------------------------------------------------------
void dmrg_su2::sweep_to_left(){
//------------------------------------------------------------------------------
  int i;
  tensor_su2 vec;
  for(i=ns/2;i>=nfree+2;i--){
    wavefunc_transformation(i-1,1);
    prepare_site_operator_from_right(i,1);
    prepare_input_vector(i-2,i-1,vec);
    prepare_site_operator_from_left(i-2);
    prepare_site_operator_from_right(i-1);
    lanczos_solve_eigenvector_idmrg(i-2,i-1,vec);
    if(i>nfree+2) vec.svd(uu[i-2],1,uu[i-1],0,wtmp);
    else if(i==nfree+2) vec.svd(uu[i-2],0,uu[i-1],1,wtmp);
  }
  for(i=nfree;i<=ns/2-2;i++){
    wavefunc_transformation(i+1,0);
    if(i!=ns/2-2)
      prepare_site_operator_from_left(i,1);
    else
      prepare_site_operator_from_left(i,0);
    if(i!=ns/2-2){
      prepare_input_vector(i+1,i+2,vec);
      prepare_site_operator_from_left(i+1);
      prepare_site_operator_from_right(i+2);
      lanczos_solve_eigenvector_idmrg(i+1,i+2,vec);
      vec.svd(uu[i+1],0,uu[i+2],1,wtmp);
    }
  }
}

//------------------------------------------------------------------------------
void dmrg_su2::sweep_to_right(){
//------------------------------------------------------------------------------
  int i;
  tensor_su2 vec;
  for(i=ns/2-1;i<ns-nfree-2;i++){
    wavefunc_transformation(i+1,0);
    prepare_site_operator_from_left(i,1);
    prepare_input_vector(i+1,i+2,vec);
    prepare_site_operator_from_left(i+1);
    prepare_site_operator_from_right(i+2);
    lanczos_solve_eigenvector_idmrg(i+1,i+2,vec);
    if(i<ns-nfree-3) vec.svd(uu[i+1],0,uu[i+2],1,wtmp);
    else if(i==ns-nfree-3) vec.svd(uu[i+1],1,uu[i+2],0,wtmp);
  }
  for(i=ns-nfree-1;i>=ns/2+1;i--){
    wavefunc_transformation(i-1,1);
    if(i!=ns/2+1)
      prepare_site_operator_from_right(i,1);
    else
      prepare_site_operator_from_right(i,0);
    if(i!=ns/2+1){
      prepare_input_vector(i-2,i-1,vec);
      prepare_site_operator_from_left(i-2);
      prepare_site_operator_from_right(i-1);
      lanczos_solve_eigenvector_idmrg(i-2,i-1,vec);
      vec.svd(uu[i-2],1,uu[i-1],0,wtmp);
    }
  }
}

//------------------------------------------------------------------------------
void dmrg_su2::prepare_site_operator_from_left(int i,int flag){
//------------------------------------------------------------------------------
  int j,k,l,m,rb,m0,m1,m2,m3,m4,m5;
  tensor_su2 htmp,vec,uu1,tconj,*hh2,uuspinor[3];
  double fac;
  double t3,t4;

  if(i==0){
    opr[i][i].operator_initial(uu[i],uu[i],sigma[0],fac_operator_onsite_left, fac_operator_onsite_rght, 0);
    for(j=0;j<exci;j++){
      ovlp[j][i].overlap_initial(orth[j][i],uu[i],0);
    }
    hh[i].clean();
    return;
  }

  hh2=new tensor_su2[psize];
  release_memory2(i+1,1);
  release_memory2(i+2,1);
  if(flag){
    read_memory(i+3,1);
    read_orth(i+1);
    read_orth(i+2);
  }
  if(hh[i-1].get_nbond()!=0)
    hh[i].overlap_transformation(uu[i],uu[i],hh[i-1],0);
  else hh[i].clean();
  m1=i;
  m2=m1+i;
  m3=m2+1;
  m4=m3+exci;
  m5=m4+ns*2;
  //m5=m4+ns;
#pragma omp parallel for default(shared) private(j,k,m,m0,myrank,fac,htmp,tconj) schedule(dynamic,1)
  for(m0=0;m0<m5;m0++){
    myrank=omp_get_thread_num();
    if(m0<m1){
      m=m0;
      k=m;
      if(fll[i][k]==1)
	opr[i][k].operator_transformation(uu[i],uu[i],opr[i-1][k],fac_operator_transformation_left,fac_operator_transformation_rght,0);
    }
    else if(m0>=m1&&m0<m2){
      m=m0-m1;
      k=m;
      if(hmap[i][k]!=0){
	htmp.operator_pairup(uu[i],uu[i],opr[i-1][k],sigma[1],fac_operator_pairup_left, fac_operator_pairup_rght, 0);
	for(j=0;j<ns;j++)
	  if(plqflg[i][j]&&plqpos[j][0]==k&&plqpos[j][1]==i){
	    plqopr[i][j]=htmp;
	    break;
	  }
	  //if(qfll[i]!=-1&&plqpos[qfll[i]][0]==k)
	  //plqopr[i][qfll[i]]=htmp;
	if(hmap[i][k]==1)fac=1-qdelta;
	else if(hmap[i][k]==2)fac=delta+qdelta/4.;
	else if(hmap[i][k]==3)fac=delta;
	htmp*=fac;
	if(hh2[myrank].get_nbond()==0)hh2[myrank]=htmp;
	else hh2[myrank]+=htmp;
      }
    }
    else if(m0>=m2&&m0<m3)
      opr[i][i].operator_initial(uu[i],uu[i],sigma[0],fac_operator_onsite_left,fac_operator_onsite_rght,0);
    else if(m0>=m3&&m0<m4){
      j=m0-m3;
      ovlp[j][i].overlap_transformation(orth[j][i],uu[i],ovlp[j][i-1],0);
      if(flag==1){//compute ovlp[j][i+1] and ovlp[j][i+1]
	ovlp[j][i+1].contract(ovlp[j][i],0,orth[j][i+1],0);
	ovlp[j][i+2].contract(orth[j][i+2],1,ovlp[j][i+3],0);
	ovlp[j][i+2].shift(1,0);
      }
    }
    else if(m0>=m4&&m0<m5){
      m=m0-m4;
      j=(m0-m4)%ns;
      k=(m0-m4)/ns;
      if(plqflg[i][j]==0)continue;
      if(i<plqpos[j][1])continue;
      if(k==0){
	if(i==plqpos[j][1])continue;
	//plqopr[i][m].operator_pairup(uu[i],uu[i],opr[i-1][plqpos[j][0]],sigma[1],0);
	else if(i==plqpos[j][2])
	  plqopr[i][m].operator_initial(uu[i],uu[i],sigma[0],plqopr[i-1][m],fac_operator_onsite_left,fac_operator_onsite_rght,0);
	else if(i==plqpos[j][3]){
	  plqopr[i][m].operator_pairup(uu[i],uu[i],plqopr[i-1][m],sigma[1],fac_operator_pairup_left, fac_operator_pairup_rght, 0);
	  plqopr[i][m]*=-qdelta;
	  if(hh2[myrank].get_nbond()==0)hh2[myrank]=plqopr[i][m];
	  else hh2[myrank]+=plqopr[i][m];
	}
	else if(i>plqpos[j][1]&&i<plqpos[j][2])
	  plqopr[i][m].overlap_transformation(uu[i],uu[i],plqopr[i-1][m],0);
	else if(i>plqpos[j][2]&&i<plqpos[j][3])
	  plqopr[i][m].operator_transformation(uu[i],uu[i],plqopr[i-1][m],fac_operator_transformation_left,fac_operator_transformation_rght,0);
      }
      else if(k==1){
	if(i==plqpos[j][1])
	  plqopr[i][m].permutation(uu[i],uu[i],opr[i-1][plqpos[j][0]],qterm[0],qterm[1],fac_permutation_left, fac_permutation_rght,0);
	else if(i==plqpos[j][2])
	  plqopr[i][m].permutation(uu[i],uu[i],plqopr[i-1][m],qterm[2],qterm[3],fac_permutation_left, fac_permutation_rght,0);
	else if(i==plqpos[j][3]){
	  plqopr[i][m].operator_pairup(uu[i],uu[i],plqopr[i-1][m],sigma[1],fac_operator_pairup_left, fac_operator_pairup_rght, 0);  
	  plqopr[i][m]*=-qdelta;
	  if(hh2[myrank].get_nbond()==0)hh2[myrank]=plqopr[i][m];
	  else hh2[myrank]+=plqopr[i][m];
	}
	else if(i>plqpos[j][1]&&i<plqpos[j][2]||i>plqpos[j][2]&&i<plqpos[j][3])
	  plqopr[i][m].operator_transformation(uu[i],uu[i],plqopr[i-1][m],fac_operator_transformation_left,fac_operator_transformation_rght,0);
      }
    }
  }
  save_memory(i-1,0);
  release_memory_all(i-1,0);
  for(j=0;j<psize;j++)
    if(hh2[j].get_nbond()!=0){
      if(hh[i].get_nbond()==0)
	hh[i]=hh2[j];
      else
	hh[i]+=hh2[j];
    }
  delete []hh2;
}

//------------------------------------------------------------------------------
void dmrg_su2::prepare_site_operator_from_right(int i,int flag){
//------------------------------------------------------------------------------
  int j,k,l,m,rb,m0,m1,m2,m3,m4,m5;
  tensor_su2 htmp,vec,tconj,*hh2,uuspinor[3];
  double fac;
  double t3,t4;
  if(i==ns-1){
    opr[i][i].operator_initial(uu[i],uu[i],sigma[1],fac_operator_onsite_left,fac_operator_onsite_rght, 1);
    for(j=0;j<exci;j++)
      ovlp[j][i].overlap_initial(orth[j][i],uu[i], 1);
    hh[i].clean();
    return;
  }
  hh2=new tensor_su2[psize];
  release_memory2(i-1,0);
  release_memory2(i-2,0);
  if(flag){
    read_memory(i-3,0);
    read_orth(i-1);
    read_orth(i-2);
  }
  if(hh[i+1].get_nbond()!=0)
    hh[i].overlap_transformation(uu[i],uu[i],hh[i+1],1);
  else hh[i].clean();
  m1=(ns-1-i);
  m2=m1+(ns-1-i);
  m3=m2+1;
  m4=m3+exci;
  m5=m4+ns*2;
  //m5=m4+ns;
#pragma omp parallel for default(shared) private(m,m0,j,k,myrank,fac,htmp,tconj) schedule(dynamic,1)
  for(m0=0;m0<m5;m0++){
    myrank=omp_get_thread_num();
    if(m0<m1){
      m=m0;
      k=i+1+m;
      if(frr[i][k]==1)
	opr[i][k].operator_transformation(uu[i],uu[i],opr[i+1][k],fac_operator_transformation_left,fac_operator_transformation_rght,1);
    }
    else if(m0>=m1&&m0<m2){
      m=m0-m1;
      k=i+1+m;
      if(hmap[i][k]!=0){
	htmp.operator_pairup(uu[i],uu[i],opr[i+1][k],sigma[0],fac_operator_pairup_left, fac_operator_pairup_rght, 1);
	for(j=0;j<ns;j++)
	  if(plqflg[i][j]&&plqpos[j][3]==k&&plqpos[j][2]==i){
	    plqopr[i][j]=htmp;
	    break;
	  }
	//if(qfrr[i]!=-1&&plqpos[qfrr[i]][3]==k)
	//plqopr[i][qfrr[i]]=htmp;
	if(hmap[i][k]==1)fac=1-qdelta;
	else if(hmap[i][k]==2)fac=delta+qdelta/4.;
	else if(hmap[i][k]==3)fac=delta;
	htmp*=fac;
	if(hh2[myrank].get_nbond()==0)hh2[myrank]=htmp;
	else hh2[myrank]+=htmp;
      }
    }
    else if(m0>=m2&&m0<m3)
      opr[i][i].operator_initial(uu[i],uu[i],sigma[1],fac_operator_onsite_left,fac_operator_onsite_rght, 1);
    else if(m0>=m3&&m0<m4){
      j=m0-m3;
      ovlp[j][i].overlap_transformation(orth[j][i],uu[i],ovlp[j][i+1],1);
      if(flag==1){
	ovlp[j][i-2].contract(ovlp[j][i-3],0,orth[j][i-2],0);
	ovlp[j][i-1].contract(orth[j][i-1],1,ovlp[j][i],0);
	ovlp[j][i-1].shift(1,0);
      }
    }
    else if(m0>=m4&&m0<m5){
      m=m0-m4;
      j=(m0-m4)%ns;
      k=(m0-m4)/ns;
      if(plqflg[i][j]==0)continue;
      if(i>plqpos[j][2])continue;
      if(k==0){
	if(i==plqpos[j][2])continue;
	  //plqopr[i][m].operator_pairup(uu[i],uu[i],opr[i+1][plqpos[j][3]],sigma[0],fac_operator_pairup_left, fac_operator_pairup_rght, 1);
	else if(i==plqpos[j][1])
	  plqopr[i][m].operator_initial(uu[i],uu[i],sigma[1],plqopr[i+1][m],fac_operator_onsite_left,fac_operator_onsite_rght, 1);
	else if(i==plqpos[j][0]){
	  plqopr[i][m].operator_pairup(uu[i],uu[i],plqopr[i+1][m],sigma[0],fac_operator_pairup_left, fac_operator_pairup_rght, 1);
	  plqopr[i][m]*=-qdelta;
	  if(hh2[myrank].get_nbond()==0)hh2[myrank]=plqopr[i][m];
	  else hh2[myrank]+=plqopr[i][m];
	}
	else if(i>plqpos[j][1]&&i<plqpos[j][2])
	  plqopr[i][m].overlap_transformation(uu[i],uu[i],plqopr[i+1][m],1);
	else if(i>plqpos[j][0]&&i<plqpos[j][1])
	  plqopr[i][m].operator_transformation(uu[i],uu[i],plqopr[i+1][m],fac_operator_transformation_left,fac_operator_transformation_rght,1);
      }
      else if(k==1){
	if(i==plqpos[j][2])
	  plqopr[i][m].permutation(uu[i],uu[i],opr[i+1][plqpos[j][3]],qterm[2+4],qterm[3+4],fac_permutation_left, fac_permutation_rght,1);
	else if(i==plqpos[j][1])
	  plqopr[i][m].permutation(uu[i],uu[i],plqopr[i+1][m],qterm[0+4],qterm[1+4],fac_permutation_left, fac_permutation_rght,1);
	else if(i==plqpos[j][0]){
	  plqopr[i][m].operator_pairup(uu[i],uu[i],plqopr[i+1][m],sigma[0],fac_operator_pairup_left, fac_operator_pairup_rght, 1);
	  plqopr[i][m]*=-qdelta;
	  if(hh2[myrank].get_nbond()==0)hh2[myrank]=plqopr[i][m];
	  else hh2[myrank]+=plqopr[i][m];
	}
	else if(i>plqpos[j][1]&&i<plqpos[j][2]||i>plqpos[j][0]&&i<plqpos[j][1])
	  plqopr[i][m].operator_transformation(uu[i],uu[i],plqopr[i+1][m],fac_operator_transformation_left,fac_operator_transformation_rght,1);
      }
    }
  }
  save_memory(i+1,1);
  release_memory_all(i+1,1);
  for(j=0;j<psize;j++)
    if(hh2[j].get_nbond()!=0){
      if(hh[i].get_nbond()==0)
	hh[i]=hh2[j];
      else
	hh[i]+=hh2[j];
    }
  delete []hh2;
}

//------------------------------------------------------------------------------
void dmrg_su2::prepare_site_operator_from_left(int i){
//------------------------------------------------------------------------------
  int j,k,l,m,rb,m0,m1,m2,m3,m4;
  tensor_su2 htmp,uu1,*hh2;
  double fac;
  su2bond bb[2];
  uu[i].get_su2bond(0,bb[0]);
  uu[i].get_su2bond(1,bb[1]);
  uu1.fuse(bb[0],bb[1]);  

  hh2=new tensor_su2[psize];
  if(hh[i-1].get_nbond()!=0)
    hh[i].overlap_transformation(uu1,uu1,hh[i-1],0);
  else hh[i].clean();
  m1=i;
  m2=m1+i;
  m3=m2+1;
  m4=m3+ns*2;
  //m4=m3+ns;
#pragma omp parallel for default(shared) private(j,k,m,m0,myrank,fac,htmp) schedule(dynamic,1)
  for(m0=0;m0<m4;m0++){
    myrank=omp_get_thread_num();
    if(m0<m1){
      m=m0;
      k=m;
      if(fll[i][k]==1)
	opr[i][k].operator_transformation(uu1,uu1,opr[i-1][k],fac_operator_transformation_left,fac_operator_transformation_rght,0);
    }
    else if(m0>=m1&&m0<m2){
      m=m0-m1;
      k=m;
      if(hmap[i][k]!=0){
	htmp.operator_pairup(uu1,uu1,opr[i-1][k],sigma[1],fac_operator_pairup_left, fac_operator_pairup_rght, 0);
	for(j=0;j<ns;j++)
	  if(plqflg[i][j]&&plqpos[j][0]==k&&plqpos[j][1]==i){
	    plqopr[i][j]=htmp;
	    break;
	  }
	//if(qfll[i]!=-1&&plqpos[qfll[i]][0]==k)
	//plqopr[i][qfll[i]]=htmp;
	if(hmap[i][k]==1)fac=1-qdelta;
	else if(hmap[i][k]==2)fac=delta+qdelta/4.;
	else if(hmap[i][k]==3)fac=delta;
	htmp*=fac;
	if(hh2[myrank].get_nbond()==0)hh2[myrank]=htmp;
	else hh2[myrank]+=htmp;
      }
    }
    else if(m0>=m2&&m0<m3)
      opr[i][i].operator_initial(uu1,uu1,sigma[0],fac_operator_onsite_left,fac_operator_onsite_rght, 0);
    else if(m0>=m3&&m0<m4){
      m=m0-m3;
      j=(m0-m3)%ns;
      k=(m0-m3)/ns;
      if(plqflg[i][j]==0)continue;
      if(i<plqpos[j][1])continue;
      if(k==0){
	if(i==plqpos[j][1])continue;
	//plqopr[i][m].operator_pairup(uu1,uu1,opr[i-1][plqpos[j][0]],sigma[1],0);
	else if(i==plqpos[j][2])
	  plqopr[i][m].operator_initial(uu1,uu1,sigma[0],plqopr[i-1][m],fac_operator_onsite_left,fac_operator_onsite_rght, 0);
	else if(i==plqpos[j][3]){
	  plqopr[i][m].operator_pairup(uu1,uu1,plqopr[i-1][m],sigma[1],fac_operator_pairup_left, fac_operator_pairup_rght, 0);
	  plqopr[i][m]*=-qdelta;
	  if(hh2[myrank].get_nbond()==0)hh2[myrank]=plqopr[i][m];
	  else hh2[myrank]+=plqopr[i][m];
	}
	else if(i>plqpos[j][1]&&i<plqpos[j][2])
	  plqopr[i][m].overlap_transformation(uu1,uu1,plqopr[i-1][m],0);
	else if(i>plqpos[j][2]&&i<plqpos[j][3])
	  plqopr[i][m].operator_transformation(uu1,uu1,plqopr[i-1][m],fac_operator_transformation_left,fac_operator_transformation_rght,0);
      }
      else if(k==1){
	if(i==plqpos[j][1])
	  plqopr[i][m].permutation(uu1,uu1,opr[i-1][plqpos[j][0]],qterm[0],qterm[1],fac_permutation_left, fac_permutation_rght,0);
	else if(i==plqpos[j][2])
	  plqopr[i][m].permutation(uu1,uu1,plqopr[i-1][m],qterm[2],qterm[3],fac_permutation_left, fac_permutation_rght,0);
	else if(i==plqpos[j][3]){
	  plqopr[i][m].operator_pairup(uu1,uu1,plqopr[i-1][m],sigma[1],fac_operator_pairup_left, fac_operator_pairup_rght, 0);
	  plqopr[i][m]*=-qdelta;
	  if(hh2[myrank].get_nbond()==0)hh2[myrank]=plqopr[i][m];
	  else hh2[myrank]+=plqopr[i][m];
	}
	else if(i>plqpos[j][1]&&i<plqpos[j][2]||i>plqpos[j][2]&&i<plqpos[j][3])
	  plqopr[i][m].operator_transformation(uu1,uu1,plqopr[i-1][m],fac_operator_transformation_left,fac_operator_transformation_rght,0);
      }
    }
  }
  for(j=0;j<psize;j++)
    if(hh2[j].get_nbond()!=0){
      if(hh[i].get_nbond()==0)
	hh[i]=hh2[j];
      else
	hh[i]+=hh2[j];
    }
  delete []hh2;
}

//------------------------------------------------------------------------------
void dmrg_su2::prepare_site_operator_from_right(int i){
//------------------------------------------------------------------------------
  int j,k,l,m,m0,m1,m2,m3,m4;
  tensor_su2 htmp,*hh2,uu1;
  double fac;
  su2bond bb[2];
  uu[i].get_su2bond(0,bb[0]);
  uu[i].get_su2bond(1,bb[1]);
  uu1.fuse(bb[0],bb[1]);

  hh2=new tensor_su2[psize];
  if(hh[i+1].get_nbond()!=0)
    hh[i].overlap_transformation(uu1,uu1,hh[i+1],1);
  else hh[i].clean();
  m1=(ns-1-i);
  m2=m1+(ns-1-i);
  m3=m2+1;
  m4=m3+ns*2;
  //m4=m3+ns;
#pragma omp parallel for default(shared) private(m,m0,j,k,myrank,fac,htmp) schedule(dynamic,1)
  for(m0=0;m0<m4;m0++){
    myrank=omp_get_thread_num();
    if(m0<m1){
      m=m0;
      k=i+1+m;
      if(frr[i][k]==1)
	opr[i][k].operator_transformation(uu1,uu1,opr[i+1][k],fac_operator_transformation_left,fac_operator_transformation_rght,1);
    }
    else if(m0>=m1&&m0<m2){
      m=m0-m1;
      k=i+1+m;
      if(hmap[i][k]!=0){
	htmp.operator_pairup(uu1,uu1,opr[i+1][k],sigma[0],fac_operator_pairup_left, fac_operator_pairup_rght, 1);
	for(j=0;j<ns;j++)
	  if(plqflg[i][j]&&plqpos[j][3]==k&&plqpos[j][2]==i){
	    plqopr[i][j]=htmp;
	    break;
	  }
	//if(qfrr[i]!=-1&&plqpos[qfrr[i]][3]==k)
	//plqopr[i][qfrr[i]]=htmp;
	if(hmap[i][k]==1)fac=1-qdelta;
	else if(hmap[i][k]==2)fac=delta+qdelta/4.;
	else if(hmap[i][k]==3)fac=delta;
	htmp*=fac;
	if(hh2[myrank].get_nbond()==0)hh2[myrank]=htmp;
	else hh2[myrank]+=htmp;
      }
    }
    else if(m0>=m2&&m0<m3)
      opr[i][i].operator_initial(uu1,uu1,sigma[1],fac_operator_onsite_left,fac_operator_onsite_rght, 1);
    else if(m0>=m3&&m0<m4){
      m=m0-m3;
      j=(m0-m3)%ns;
      k=(m0-m3)/ns;
      if(plqflg[i][j]==0)continue;
      if(i>plqpos[j][2])continue;
      if(k==0){
	if(i==plqpos[j][2])continue;
	//plqopr[i][m].operator_pairup(uu1,uu1,opr[i+1][plqpos[j][3]],sigma[0],1);
	else if(i==plqpos[j][1])
	  plqopr[i][m].operator_initial(uu1,uu1,sigma[1],plqopr[i+1][m],fac_operator_onsite_left,fac_operator_onsite_rght, 1);
	else if(i==plqpos[j][0]){
	  plqopr[i][m].operator_pairup(uu1,uu1,plqopr[i+1][m],sigma[0],fac_operator_pairup_left, fac_operator_pairup_rght, 1);
	  plqopr[i][m]*=-qdelta;
	  if(hh2[myrank].get_nbond()==0)hh2[myrank]=plqopr[i][m];
	  else hh2[myrank]+=plqopr[i][m];
	}
	else if(i>plqpos[j][1]&&i<plqpos[j][2])
	  plqopr[i][m].overlap_transformation(uu1,uu1,plqopr[i+1][m],1);
	else if(i>plqpos[j][0]&&i<plqpos[j][1])
	  plqopr[i][m].operator_transformation(uu1,uu1,plqopr[i+1][m],fac_operator_transformation_left,fac_operator_transformation_rght,1);
      }
      else if(k==1){
	if(i==plqpos[j][2])
	  plqopr[i][m].permutation(uu1,uu1,opr[i+1][plqpos[j][3]],qterm[2+4],qterm[3+4],fac_permutation_left, fac_permutation_rght,1);
	else if(i==plqpos[j][1])
	  plqopr[i][m].permutation(uu1,uu1,plqopr[i+1][m],qterm[0+4],qterm[1+4],fac_permutation_left, fac_permutation_rght,1);
	else if(i==plqpos[j][0]){
	  plqopr[i][m].operator_pairup(uu1,uu1,plqopr[i+1][m],sigma[0],fac_operator_pairup_left, fac_operator_pairup_rght, 1);
	  plqopr[i][m]*=-qdelta;
	  if(hh2[myrank].get_nbond()==0)hh2[myrank]=plqopr[i][m];
	  else hh2[myrank]+=plqopr[i][m];
	}
	else if(i>plqpos[j][1]&&i<plqpos[j][2]||i>plqpos[j][0]&&i<plqpos[j][1])
	  plqopr[i][m].operator_transformation(uu1,uu1,plqopr[i+1][m],fac_operator_transformation_left,fac_operator_transformation_rght,1);
      }
    }
  }
  for(j=0;j<psize;j++)
    if(hh2[j].get_nbond()!=0){
      if(hh[i].get_nbond()==0)
	hh[i]=hh2[j];
      else
	hh[i]+=hh2[j];
    }
  delete []hh2;
}

//------------------------------------------------------------------------------
void dmrg_su2::send_tensor_su2(tensor_su2& tu1,int des){
//------------------------------------------------------------------------------
  int i,j,nele1,nele2,nb,lc,icnt,*isend,dir,nc,ishift,*iflag,nten;
  double *ptr1,*ptr2;
  tensor_su2 *tptr;

  tptr=&(tu1);
  nb=tptr->get_nbond();
  lc=tptr->get_locspin();
  icnt=4;
  for(i=0;i<nb;i++){
    nc=tptr->get_nmoment(i);
    icnt+=nc*2+2;  
  }
  isend=new int[icnt];
  isend[0]=nb;
  isend[1]=lc;
  for(i=0;i<nb;i++){
    nc=tptr->get_nmoment(i);
    dir=tptr->get_bonddir(i);
    isend[2+i*2]=dir;
    isend[3+i*2]=nc;
  }
  ishift=2+nb*2;
  for(i=0;i<nb;i++){
    nc=isend[3+i*2];
    for(j=0;j<nc;j++)
      isend[ishift+j]=tptr->get_angularmoment(i,j);
    for(j=0;j<nc;j++)
      isend[ishift+nc+j]=tptr->get_bonddim(i,j);
    ishift+=nc*2;    
  }
  tptr->get_nelement(nele1,nele2);
  isend[icnt-2]=nele1;
  isend[icnt-1]=nele2;
  ptr1=new double[nele1];
  ptr2=new double[nele2];
  nten=tptr->get_nten();
  iflag=new int[nten];
  tptr->get_telement(ptr1,ptr2,iflag);
  //start communication
  MPI::COMM_WORLD.Send(&icnt,1,MPI_INT,des,0);
  MPI::COMM_WORLD.Send(isend,icnt,MPI_INT,des,1);
  MPI::COMM_WORLD.Send(ptr1,nele1,MPI_DOUBLE,des,2);
  MPI::COMM_WORLD.Send(ptr2,nele2,MPI_DOUBLE,des,3);
  MPI::COMM_WORLD.Send(iflag,nten,MPI_INT,des,4);
  delete []iflag;
  delete []isend;
  delete []ptr1;
  delete []ptr2;
}

//------------------------------------------------------------------------------
void dmrg_su2::recv_tensor_su2(tensor_su2& tu1,int src){
//------------------------------------------------------------------------------
  int i,j,k,nele1,nele2,nb,lc,icnt,*isend,dir,nc,*cval,*bdim,tag,ishift,*iflag,nten;
  double *ptr1,*ptr2;
  su2struct u1;
  su2bond *b1;
  tensor_su2 *tptr;
  tptr=&(tu1);
  MPI::COMM_WORLD.Recv(&icnt,1,MPI_INT,src,0);
  isend=new int[icnt];
  MPI::COMM_WORLD.Recv(isend,icnt,MPI_INT,src,1);
  nele1=isend[icnt-2];
  nele2=isend[icnt-1];
  ptr1=new double[nele1];
  ptr2=new double[nele2];
  MPI::COMM_WORLD.Recv(ptr1,nele1,MPI_DOUBLE,src,2);
  MPI::COMM_WORLD.Recv(ptr2,nele2,MPI_DOUBLE,src,3);
  nb=isend[0];
  lc=isend[1];
  b1=new su2bond[nb];
  ishift=2+nb*2;
  for(i=0;i<nb;i++){
    dir=isend[2+i*2];
    nc=isend[3+i*2];
    cval=new int[nc];
    bdim=new int[nc];
    for(j=0;j<nc;j++)
      cval[j]=isend[ishift+j];
    for(j=0;j<nc;j++)
      bdim[j]=isend[ishift+nc+j];
    ishift+=nc*2;
    b1[i].set_su2bond(nc,dir,cval,bdim);
    delete []cval;
    delete []bdim;
  }
  u1.set_su2struct(nb,lc,b1);
  nten=u1.get_nten();
  iflag=new int[nten];
  MPI::COMM_WORLD.Recv(iflag,nten,MPI_INT,src,4);
  delete []b1;
  tptr->set_tensor_su2(u1,ptr1,ptr2,iflag);
  delete []iflag;
  delete []isend;
  delete []ptr1;
  delete []ptr2;
}

//------------------------------------------------------------------------------
void dmrg_su2::synchronize_site(int i, int lr, int src, int des){
//------------------------------------------------------------------------------
  int j,k,tag;
  int n1,n2;
  if(lr==0){
    n1=0;
    n2=i+1;
  }
  else if(lr==1){
    n1=i;
    n2=ns;
  }
  for(k=0;k<2+exci;k++){
    if(comm_rank==src){
      if(k==0)
	send_tensor_su2(uu[i],des);
      else if(k==1){
	if(hh[i].is_null()) tag=1;
	else tag=0;
	MPI::COMM_WORLD.Send(&tag,1,MPI_INT,des,0);
	if(tag==0) send_tensor_su2(hh[i],des);
      }
      else if(k>1)
	send_tensor_su2(ovlp[k-2][i],des);
    }
    else if(comm_rank==des){
      if(k==0)
	recv_tensor_su2(uu[i],src);
      else if(k==1){
	MPI::COMM_WORLD.Recv(&tag,1,MPI_INT,src,0);
	if(tag==0) recv_tensor_su2(hh[i],src);
      }
      else if(k>1)
	recv_tensor_su2(ovlp[k-2][i],src);
    }
  }
  for(k=n1;k<n2;k++)
    if(k==i||(lr==0&&fll[i][k]!=0)||(lr==1&&frr[i][k]!=0))
      if(comm_rank==src){
	MPI::COMM_WORLD.Barrier();
	send_tensor_su2(opr[i][k],des);
      }
      else if(comm_rank==des){
	MPI::COMM_WORLD.Barrier();
	recv_tensor_su2(opr[i][k],src);
      }
  for(k=0;k<ns;k++){
    if(plqflg[i][k]==0)continue;
    if(src==1&&i>plqpos[k][2])continue;
    if(src==0&&i<plqpos[k][1])continue;
    if(comm_rank==src){
      MPI::COMM_WORLD.Barrier();
      send_tensor_su2(plqopr[i][k],des);
      send_tensor_su2(plqopr[i][k+ns],des);
    }
    else if(comm_rank==des){
      MPI::COMM_WORLD.Barrier();
      recv_tensor_su2(plqopr[i][k],src);
      recv_tensor_su2(plqopr[i][k+ns],src);
    }
  }
}

//------------------------------------------------------------------------------
void dmrg_su2::save_mps1(){
//------------------------------------------------------------------------------
  int i,j,k,nele,nb,nc,n1,n2;
  double *ptr;
  char base[100],name[100],len[10],dim[10],sec[10],pos[10],exc[10],id[10],dk[10];
  ofstream fout;
  tensor_su2* tptr;
  if(memory_flag==1){
    strcpy(name,"rsync -avz ./data");
    //sprintf(id,"%d",jobid);
    sprintf(dk,"%d",comm_rank+3);
    strcat(name,dk);
    //strcat(name,"/");
    //strcat(name,id);
    strcat(name,"/uu-* ./");
    system(name);
    system(name);
    system(name);
  }
}

//------------------------------------------------------------------------------
void dmrg_su2::save_mps2(){
//------------------------------------------------------------------------------
  int i,j,k,nele,nele1,nb,nc,n1,n2,nten,*iflag;
  double *ptr,*ptr1;
  char base[100],name[100],len[10],dim[10],sec[10],pos[10],exc[10],id[10],dk[10];
  ofstream fout;
  tensor_su2* tptr;
  if(memory_flag==1){
    if(comm_rank==0){
      n1=ns/2-2;
      n2=ns/2;
    }
    else if(comm_rank==1){
      n1=ns/2;
      n2=ns/2+2;
    }
  }
  else{
    if(comm_rank==0){
      n1=0;n2=ns/2;
    }
    else if(comm_rank==1){
      n1=ns/2;n2=ns;
    }
  }
  sprintf(len,"%d",ns);
  sprintf(sec,"%d",totspin);
  sprintf(dim,"%d",max_dcut);
  sprintf(exc,"%d",exci);

  strcpy(name,"uu-");
  strcat(name,len);
  strcat(name,"-");
  strcat(name,dim);
  strcat(name,"-");
  strcat(name,sec);
  strcat(name,"-");
  strcat(name,exc);
  strcat(name,"-");
  strcpy(base,name);
  for(i=n1;i<n2;i++){
    tptr=&(uu[i]);
    sprintf(pos,"%d",i);
    strcpy(name,base);
    strcat(name,pos);
    strcat(name,".dat");
    fout.open(name,ios::out);
    if(fout.is_open()){
      fout<<tptr->get_nbond()<<"\t"<<tptr->get_locspin()<<endl;
      nb=tptr->get_nbond();
      for(j=0;j<nb;j++){
	fout<<tptr->get_bonddir(j)<<"\t"<<tptr->get_nmoment(j)<<endl;
	nc=tptr->get_nmoment(j);
	for(k=0;k<nc;k++)
	  fout<<tptr->get_angularmoment(j,k)<<"\t";
	fout<<endl;
	for(k=0;k<nc;k++)
	  fout<<tptr->get_bonddim(j,k)<<"\t";
	fout<<endl;
      }
      fout.close();
    }
    sprintf(pos,"%d",i);
    strcpy(name,base);
    strcat(name,pos);
    strcat(name,".bin");
    fout.open(name,ios::out | ios::binary);
    if(fout.is_open()){
      tptr->get_nelement(nele,nele1);
      ptr=new double[nele+nele1];
      ptr1=&(ptr[nele]);
      nten=tptr->get_nten();
      iflag=new int[nten];
      tptr->get_telement(ptr,ptr1,iflag);
      fout.write((char*)iflag,nten*sizeof(int));
      fout.write((char*)ptr,(nele+nele1)*sizeof(double));
      delete []ptr;
      delete []iflag;
      fout.close();
    }      
  }
}

//------------------------------------------------------------------------------
bool dmrg_su2::read_mps(int read1,int read2){
//------------------------------------------------------------------------------
  int i,j,k,state,nele,nele1,bd[3],nb,nc,lc,dir,*angm,*bdim,read,n1,n2,*iflag,nten;
  double *ptr,*ptr1;
  char base[100],name[100],len[10],dim[10],sec[10],pos[10],exc[10];
  ifstream fin;
  tensor_su2* tptr;
  su2bond* b1;
  su2struct* u1;
  bool pass=true,pass1=true;

  if(comm_rank==0){
    n1=0;n2=ns/2;
  }
  else if(comm_rank==1){
    n1=ns/2;n2=ns;
  }
  u1=new su2struct[ns];
  sprintf(len,"%d",ns);
  sprintf(sec,"%d",totspin);

  for(state=0;state<=exci;state++){
    if(state==exci)
      read=read1;
    else if(state<exci){
      cout<<comm_rank<<" read mps "<<read2<<endl;
      read=read2;
      //if(state==0&&ly==10)read=3000;
    }
    sprintf(dim,"%d",read);
    sprintf(exc,"%d",state);

    if(read==max_dcut){
      strcpy(name,"uu-");
    }
    else if((read<max_dcut||read>max_dcut)&&read!=0){
      strcpy(name,"../");
      strcat(name,dim);
      strcat(name,"/uu-");
    }
    else if(read==0){
      pass=false;
      continue;
    }
    strcat(name,len);
    strcat(name,"-");
    strcat(name,dim);
    strcat(name,"-");
    strcat(name,sec);
    strcat(name,"-");
    strcat(name,exc);
    strcat(name,"-");
    strcpy(base,name);
    for(i=n1;i<n2;i++){
      if(state==exci)
	tptr=&(uu[i]);
      else tptr=&(orth[state][i]);
      sprintf(pos,"%d",i);
      strcpy(name,base);
      strcat(name,pos);
      strcat(name,".dat");
      fin.open(name,ios::in);
      if(fin.is_open()){
	fin>>nb>>lc;
	b1=new su2bond[nb];
	for(j=0;j<nb;j++){
	  fin>>dir>>nc;
	  angm=new int[nc];
	  bdim=new int[nc];
	  for(k=0;k<nc;k++)
	    fin>>angm[k];
	  for(k=0;k<nc;k++)
	    fin>>bdim[k];
	  b1[j].set_su2bond(nc,dir,angm,bdim);
	  delete []angm;
	  delete []bdim;
	}
	u1[i].set_su2struct(nb,lc,b1);
	delete []b1;
	fin.close();
      }
      else{
	cout<<comm_rank<<"\tstate="<<state<<"\t dat pos="<<i<<endl;
	if(state==exci)	pass=false;
	else if(state<exci) pass1=false;
      }
      strcpy(name,base);
      strcat(name,pos);
      strcat(name,".bin");
      fin.open(name,ios::in | ios::binary);
      if(fin.is_open()){
	nten=u1[i].get_nten();
	iflag=new int[nten];
	fin.read((char*)iflag,nten*sizeof(int));
	u1[i].get_nelement(nele,nele1,iflag);
	ptr=new double[nele+nele1];
	ptr1=&(ptr[nele]);
	fin.read((char*)ptr,(nele+nele1)*sizeof(double));
	tptr->set_tensor_su2(u1[i],ptr,ptr1,iflag);
	delete []ptr;
	delete []iflag;
	fin.close();
      }
      else{
	cout<<comm_rank<<"\tstate="<<state<<"\t bin pos="<<i<<endl;
	if(state==exci)
	  pass=false;    
	else if(state<exci)
	  pass1=false;
      }
    }
  }
  delete []u1;
  if(pass1==false){
    cout<<"can not do excited state without precalculated ground state"<<endl;
    exit(0);
  }
  if(pass==false){
    cout<<"can not read uu mps"<<endl;
  }
  MPI::COMM_WORLD.Barrier();
  return pass;
}

//------------------------------------------------------------------------------
void dmrg_su2::save_enr(){
//------------------------------------------------------------------------------
  char base[100],name[100],len[10],dim[10],sec[10],pos[10],exc[10];
  ofstream fout;
  if(comm_rank==1)return;
  sprintf(len,"%d",ns);
  sprintf(sec,"%d",totspin);
  sprintf(dim,"%d",max_dcut);
  sprintf(exc,"%d",exci);

  strcpy(name,"enr-");
  strcat(name,len);
  strcat(name,"-");
  strcat(name,dim);
  strcat(name,"-");
  strcat(name,sec);
  strcat(name,"-");
  strcat(name,exc);
  strcat(name,".dat");
  fout.open(name,ios::app);
  if(fout.is_open()){
    fout<<setprecision(16);
    fout<<gs_enr[exci]<<endl;
    fout.close();
  }
}

//------------------------------------------------------------------------------
void dmrg_su2::read_enr(int read1, int read2){
//------------------------------------------------------------------------------
  int state,read;
  char base[100],name[100],len[10],dim[10],sec[10],pos[10],exc[10];
  ifstream fin;
  sprintf(len,"%d",ns);
  sprintf(sec,"%d",totspin);
  for(state=0;state<exci+1;state++){
    if(state==exci)
      read=read1;
    else if(state<exci){
      cout<<comm_rank<<" read enr enr "<<read2<<endl;
      read=read2;
      //if(state==0&&ly==10)read=3000;
    }
    sprintf(dim,"%d",read);
    sprintf(exc,"%d",state);
    if(read==max_dcut)
      strcpy(name,"enr-");
    else if((read<max_dcut||read>max_dcut)&&read!=0){
      strcpy(name,"../");
      strcat(name,dim);
      strcat(name,"/enr-");
    }
    else if(state!=exci&&read==0){
      cout<<"need ground state enr to calculate excited state"<<endl;
      exit(0);
      continue;
    }
    strcat(name,len);
    strcat(name,"-");
    strcat(name,dim);
    strcat(name,"-");
    strcat(name,sec);
    strcat(name,"-");
    strcat(name,exc);
    strcat(name,".dat");
    fin.open(name,ios::in);
    if(fin.is_open()){
      do{
	fin>>gs_enr[state];
      }while(!fin.eof());
      cout<<comm_rank<<" state="<<state<<" gs_enr["<<state<<"]="<<gs_enr[state]<<endl;
      fin.close();
    }
  }
  cout<<comm_rank<<" done read enr"<<endl;
}

//------------------------------------------------------------------------------
void dmrg_su2::save_ww(){
//------------------------------------------------------------------------------
  char base[100],name[100],len[10],dim[10],sec[10],pos[10],exc[10],id[10],dk[10];
  ofstream fout;
  if(comm_rank==1)return;
  sprintf(len,"%d",ns);
  sprintf(sec,"%d",totspin);
  sprintf(dim,"%d",max_dcut);
  sprintf(exc,"%d",exci);
  sprintf(pos,"%d",ns/2);

  strcpy(name,"ww-");
  strcat(name,len);
  strcat(name,"-");
  strcat(name,dim);
  strcat(name,"-");
  strcat(name,sec);
  strcat(name,"-");
  strcat(name,exc);
  strcat(name,"-");
  strcat(name,pos);
  strcpy(base,name);
  strcpy(name,base);
  strcat(name,".bin");
  fout.open(name,ios::out | ios::binary);
  if(fout.is_open()){
    fout.write((char*)ww[exci],max_dcut*phdim*sizeof(double));    
    fout.close();
  }
}

//------------------------------------------------------------------------------
void dmrg_su2::read_ww(int read1,int read2){
//------------------------------------------------------------------------------
  char base[100],name[100],len[10],dim[10],sec[10],pos[10],exc[10],id[10],dk[10];
  ifstream fin;
  int state,i;
  sprintf(len,"%d",ns);
  sprintf(sec,"%d",totspin);
  sprintf(pos,"%d",ns/2);

  for(state=0;state<=exci;state++){
    if(state==exci)
      read=read1;
    else if(state<exci){
      read=read2;
      //if(state==0&&ly==10)read=3000;
    }
    sprintf(dim,"%d",read);
    sprintf(exc,"%d",state);

    strcpy(name,"ww-");
    if(read==max_dcut){
      strcpy(name,"ww-");
    }
    else if(read<max_dcut||read>max_dcut){
      strcpy(name,"../");
      strcat(name,dim);
      strcat(name,"/ww-");
    }
    strcat(name,len);
    strcat(name,"-");
    strcat(name,dim);
    strcat(name,"-");
    strcat(name,sec);
    strcat(name,"-");
    strcat(name,exc);
    strcat(name,"-");
    strcat(name,pos);
    strcpy(base,name);
    strcpy(name,base);
    strcat(name,".bin");
    fin.open(name,ios::in | ios::binary);
    if(fin.is_open()){
      fin.read((char*)ww[state],read*sizeof(double));    
      fin.close();
    }
    else{
      cout<<comm_rank<<"\tstate="<<state<<"\t ww set 1"<<endl;
      for(i=0;i<read;i++)
	ww[state][i]=1;
    }
  }
}

//------------------------------------------------------------------------------
void dmrg_su2::save_memory(int p,int lr){
//------------------------------------------------------------------------------
  int i,j,k,l,m,nele,nele1,nten,nb,nc,*iflag;
  double *ptr,*ptr1;
  char base[100],len[10],dim[10],sec[10],pos[10],exc[10],name[100],pos2[10],id[10],dk[10];
  ofstream fout;
  tensor_su2 *tptr;
  if(memory_flag==0)return;

  sprintf(len,"%d",ns);
  sprintf(sec,"%d",totspin);
  sprintf(dim,"%d",max_dcut);
  sprintf(exc,"%d",exci);
  sprintf(pos,"%d",p);
  sprintf(id,"%d",jobid);
  sprintf(dk,"%d",comm_rank+3);

  //#pragma omp parallel for default(shared) private(myrank,i,j,k,l,m,nb,nc,nele,ptr)
  for(l=0;l<3;l++){
    strcpy(base,"./data");
    strcat(base,dk);
    //strcat(base,"/");
    //strcat(base,id);
    if(l==0) strcat(base,"/opr-");
    else strcat(base,"/plqopr-");
    strcat(base,len);
    strcat(base,"-");
    strcat(base,dim);
    strcat(base,"-");
    strcat(base,sec);
    strcat(base,"-");
    strcat(base,exc);
    strcat(base,"-");
    for(i=0;i<ns;i++){
      if(l==0&&(lr==0&&fll[p][i]==1||lr==1&&frr[p][i]==1||i==p)){
	tptr=&(opr[p][i]);
	sprintf(pos2,"%d",i);
      }
      else if(l>0&&plqflg[p][i]!=0&&(lr==0&&p>=plqpos[i][1]||lr==1&&p<=plqpos[i][2])){
	tptr=&(plqopr[p][i+(l-1)*ns]);
	sprintf(pos2,"%d",i+(l-1)*ns);
      }
      else continue;
      strcpy(name,base);
      strcat(name,pos);
      strcat(name,"-");
      strcat(name,pos2);
      strcat(name,".dat");
      fout.open(name,ios::out);
      if(fout.is_open()){
	fout<<tptr->get_nbond()<<"\t"<<tptr->get_locspin()<<endl;
	nb=tptr->get_nbond();
	for(j=0;j<nb;j++){
	  fout<<tptr->get_bonddir(j)<<"\t"<<tptr->get_nmoment(j)<<endl;
	  nc=tptr->get_nmoment(j);
	  for(k=0;k<nc;k++)
	    fout<<tptr->get_angularmoment(j,k)<<"\t";
	  fout<<endl;
	  for(k=0;k<nc;k++)
	    fout<<tptr->get_bonddim(j,k)<<"\t";
	  fout<<endl;
	}
	fout.close();
      }
      strcpy(name,base);
      strcat(name,pos);
      strcat(name,"-");
      strcat(name,pos2);
      strcat(name,".bin");
      fout.open(name,ios::out | ios::binary);
      if(fout.is_open()){
	tptr->get_nelement(nele,nele1);
	ptr=new double[nele+nele1];
	ptr1=&(ptr[nele]);
	nten=tptr->get_nten();
	iflag=new int[nten];
	tptr->get_telement(ptr,ptr1,iflag);
	fout.write((char*)iflag,nten*sizeof(int));
	fout.write((char*)ptr,(nele+nele1)*sizeof(double));
	delete []ptr;
	delete []iflag;
	fout.close();
      }      
    }
  }
  //save hh,uu/uu,ovlp
  for(l=0;l<2+exci*2;l++){
    strcpy(base,"./data");
    strcat(base,dk);
    //strcat(base,"/");
    //strcat(base,id);
    if(l==0){
      strcat(base,"/hh-");
      if(hh[p].get_nbond()==0)continue;
    }
    else if(l==1){
      strcat(base,"/uu-");
    }
    else if(l>=2&&l<2+exci){
      strcat(base,"/ovlp-");
      sprintf(pos2,"%d",l-2);
    }
    else if(l>=2+exci&&l<2+exci*2){
      strcat(base,"/orth-");
      sprintf(pos2,"%d",l-2-exci);
    }
    strcat(base,len);
    strcat(base,"-");
    strcat(base,dim);
    strcat(base,"-");
    strcat(base,sec);
    strcat(base,"-");
    strcat(base,exc);
    strcat(base,"-");
    if(l==0)
      tptr=&(hh[p]);
    else if(l==1)
      tptr=&(uu[p]);
    else if(l>=2&&l<2+exci)
      tptr=&(ovlp[l-2][p]);
    else if(l>=2+exci&&l<2+exci*2)
      tptr=&(orth[l-2-exci][p]);
    strcpy(name,base);
    strcat(name,pos);
    if(l>=2){
      strcat(name,"-");
      strcat(name,pos2);
    }
    strcat(name,".dat");
    fout.open(name,ios::out);
    if(fout.is_open()){
      fout<<tptr->get_nbond()<<"\t"<<tptr->get_locspin()<<endl;
      nb=tptr->get_nbond();
      for(j=0;j<nb;j++){
	fout<<tptr->get_bonddir(j)<<"\t"<<tptr->get_nmoment(j)<<endl;
	nc=tptr->get_nmoment(j);
	for(k=0;k<nc;k++)
	  fout<<tptr->get_angularmoment(j,k)<<"\t";
	fout<<endl;
	for(k=0;k<nc;k++)
	  fout<<tptr->get_bonddim(j,k)<<"\t";
	fout<<endl;
      }
      fout.close();
    }
    strcpy(name,base);
    strcat(name,pos);
    if(l>=2){
      strcat(name,"-");
      strcat(name,pos2);
    }
    strcat(name,".bin");
    fout.open(name,ios::out | ios::binary);
    if(fout.is_open()){
      tptr->get_nelement(nele,nele1);
      ptr=new double[nele+nele1];
      ptr1=&(ptr[nele]);
      nten=tptr->get_nten();
      iflag=new int[nten];
      tptr->get_telement(ptr,ptr1,iflag);
      fout.write((char*)iflag,nten*sizeof(int));
      fout.write((char*)ptr,(nele+nele1)*sizeof(double));
      delete []ptr;
      delete []iflag;
      fout.close();
    }
  }      
}

//------------------------------------------------------------------------------
void dmrg_su2::read_orth(int p){
//------------------------------------------------------------------------------
  int i,j,k,l,m,state,nele,nele1,bd[3],nb,nc,lc,dir,*cval,*bdim,read,*iflag,nten;
  double *ptr,*ptr1;
  char base[100],len[10],dim[10],sec[10],pos[10],exc[10],name[100],pos2[10],id[10],dk[10];
  ifstream fin;
  tensor_su2 *tptr;
  su2bond* b1;
  su2struct u1;
  if(memory_flag==0)return;

  sprintf(len,"%d",ns);
  sprintf(sec,"%d",totspin);
  sprintf(dim,"%d",max_dcut);
  sprintf(exc,"%d",exci);
  sprintf(pos,"%d",p);
  sprintf(id,"%d",jobid);
  sprintf(dk,"%d",comm_rank+3);
  for(l=0;l<exci;l++){
    if(orth[l][p].get_nbond()!=0)continue;
    strcpy(base,"./data");
    strcat(base,dk);
    //strcat(base,"/");
    //strcat(base,id);
    strcat(base,"/orth-");
    sprintf(pos2,"%d",l);

    strcat(base,len);
    strcat(base,"-");
    strcat(base,dim);
    strcat(base,"-");
    strcat(base,sec);
    strcat(base,"-");
    strcat(base,exc);
    strcat(base,"-");

    tptr=&(orth[l][p]);
    strcpy(name,base);
    strcat(name,pos);
    strcat(name,"-");
    strcat(name,pos2);
    strcat(name,".dat");
    fin.open(name,ios::in);
    if(fin.is_open()){
      fin>>nb>>lc;
      b1=new su2bond[nb];
      for(j=0;j<nb;j++){
	fin>>dir>>nc;
	cval=new int[nc];
	bdim=new int[nc];
	for(k=0;k<nc;k++)
	  fin>>cval[k];
	for(k=0;k<nc;k++)
	  fin>>bdim[k];
	b1[j].set_su2bond(nc,dir,cval,bdim);
	delete []cval;
	delete []bdim;
      }
      u1.set_su2struct(nb,lc,b1);
      delete []b1;
      fin.close();
    }
    strcpy(name,base);
    strcat(name,pos);
    strcat(name,"-");
    strcat(name,pos2);
    strcat(name,".bin");
    fin.open(name,ios::in | ios::binary);
    if(fin.is_open()){
      nten=u1.get_nten();
      iflag=new int[nten];
      fin.read((char*)iflag,nten*sizeof(int));
      u1.get_nelement(nele,nele1,iflag);
      ptr=new double[nele+nele1];
      ptr1=&(ptr[nele]);
      fin.read((char*)ptr,(nele+nele1)*sizeof(double));
      tptr->set_tensor_su2(u1,ptr,ptr1,iflag);
      delete []ptr;
      delete []iflag;
      fin.close();
    }
  }
}

//------------------------------------------------------------------------------
void dmrg_su2::read_memory(int p,int lr){
//------------------------------------------------------------------------------
  int i,j,k,l,m,state,nele,nele1,bd[3],nb,nc,lc,dir,*cval,*bdim,read,*iflag,nten;
  double *ptr,*ptr1;
  char base[100],len[10],dim[10],sec[10],pos[10],exc[10],name[100],pos2[10],id[10],dk[10];
  ifstream fin;
  tensor_su2 *tptr,*tptrin,tmp;
  su2bond* b1;
  su2struct u1;
  if(memory_flag==0)return;

  sprintf(len,"%d",ns);
  sprintf(sec,"%d",totspin);
  sprintf(dim,"%d",max_dcut);
  sprintf(exc,"%d",exci);
  sprintf(pos,"%d",p);
  sprintf(id,"%d",jobid);
  sprintf(dk,"%d",comm_rank+3);

  //#pragma omp parallel for default(shared) private(myrank,i,j,k,l,m,lc,dir,nb,nc,cval,bdim,b1,u1,nele,ptr)
  tptrin=&(tmp);
  for(l=0;l<3;l++){
    strcpy(base,"./data");
    strcat(base,dk);
    //strcat(base,"/");
    //strcat(base,id);
    if(l==0) strcat(base,"/opr-");
    else strcat(base,"/plqopr-");
    strcat(base,len);
    strcat(base,"-");
    strcat(base,dim);
    strcat(base,"-");
    strcat(base,sec);
    strcat(base,"-");
    strcat(base,exc);
    strcat(base,"-");
    for(i=0;i<ns;i++){
      if(l==0&&(lr==0&&fll[p][i]==1||lr==1&&frr[p][i]==1||i==p)){
	tptr=&(opr[p][i]);
	sprintf(pos2,"%d",i);
      }
      else if(l>0&&plqflg[p][i]!=0&&(lr==0&&p>=plqpos[i][1]||lr==1&&p<=plqpos[i][2])){
	tptr=&(plqopr[p][i+(l-1)*ns]);
	sprintf(pos2,"%d",i+(l-1)*ns);
      }
      else continue;
      strcpy(name,base);
      strcat(name,pos);
      strcat(name,"-");
      strcat(name,pos2);
      strcat(name,".dat");
      fin.open(name,ios::in);
      if(fin.is_open()){
	fin>>nb>>lc;
	b1=new su2bond[nb];
	for(j=0;j<nb;j++){
	  fin>>dir>>nc;
	  cval=new int[nc];
	  bdim=new int[nc];
	  for(k=0;k<nc;k++)
	    fin>>cval[k];
	  for(k=0;k<nc;k++)
	    fin>>bdim[k];
	  b1[j].set_su2bond(nc,dir,cval,bdim);
	  delete []cval;
	  delete []bdim;
	}
	u1.set_su2struct(nb,lc,b1);
	delete []b1;
	fin.close();
      }
      strcpy(name,base);
      strcat(name,pos);
      strcat(name,"-");
      strcat(name,pos2);
      strcat(name,".bin");
      fin.open(name,ios::in | ios::binary);
      if(fin.is_open()){
	nten=u1.get_nten();
	iflag=new int[nten];
	fin.read((char*)iflag,nten*sizeof(int));
	u1.get_nelement(nele,nele1,iflag);
	ptr=new double[nele+nele1];
	ptr1=&(ptr[nele]);
	fin.read((char*)ptr,(nele+nele1)*sizeof(double));
	tptr->set_tensor_su2(u1,ptr,ptr1,iflag);
	//tptrin->set_tensor_su2(u1,ptr,ptr1,iflag);
	delete []ptr;
	delete []iflag;
	fin.close();
      }
    }
  }
  //read hh uu/uu ovlp
  //#pragma omp parallel for default(shared) private(myrank,i,j,k,l,m,lc,dir,nb,nc,cval,bdim,b1,u1,nele,ptr)
  for(l=0;l<2+exci*2;l++){
    strcpy(base,"./data");
    strcat(base,dk);
    //strcat(base,"/");
    //strcat(base,id);
    if(l==0){
      strcat(base,"/hh-");
    }
    else if(l==1){
      strcat(base,"/uu-");
    }
    else if(l>=2&&l<2+exci){
      strcat(base,"/ovlp-");
      sprintf(pos2,"%d",l-2);
    }
    else if(l>=2+exci&&l<2+exci*2){
      strcat(base,"/orth-");
      sprintf(pos2,"%d",l-2-exci);
    }
    strcat(base,len);
    strcat(base,"-");
    strcat(base,dim);
    strcat(base,"-");
    strcat(base,sec);
    strcat(base,"-");
    strcat(base,exc);
    strcat(base,"-");
    if(l==0)
      tptr=&(hh[p]);
    else if(l==1)
      tptr=&(uu[p]);
    else if(l>=2&&l<2+exci)
      tptr=&(ovlp[l-2][p]);
    else if(l>=2+exci&&l<2+exci*2)
      tptr=&(orth[l-2-exci][p]);
    strcpy(name,base);
    strcat(name,pos);
    if(l>=2){
      strcat(name,"-");
      strcat(name,pos2);
    }
    strcat(name,".dat");
    fin.open(name,ios::in);
    if(fin.is_open()){
      fin>>nb>>lc;
      b1=new su2bond[nb];
      for(j=0;j<nb;j++){
	fin>>dir>>nc;
	cval=new int[nc];
	bdim=new int[nc];
	for(k=0;k<nc;k++)
	  fin>>cval[k];
	for(k=0;k<nc;k++)
	  fin>>bdim[k];
	b1[j].set_su2bond(nc,dir,cval,bdim);
	delete []cval;
	delete []bdim;
      }
      u1.set_su2struct(nb,lc,b1);
      delete []b1;
      fin.close();
    }
    strcpy(name,base);
    strcat(name,pos);
    if(l>=2){
      strcat(name,"-");
      strcat(name,pos2);
    }
    strcat(name,".bin");
    fin.open(name,ios::in | ios::binary);
    if(fin.is_open()){
      nten=u1.get_nten();
      iflag=new int[nten];
      fin.read((char*)iflag,nten*sizeof(int));
      u1.get_nelement(nele,nele1,iflag);
      ptr=new double[nele+nele1];
      ptr1=&(ptr[nele]);
      fin.read((char*)ptr,(nele+nele1)*sizeof(double));
      tptr->set_tensor_su2(u1,ptr,ptr1,iflag);
      //tptrin->set_tensor_su2(u1,ptr,ptr1,iflag);
      delete []ptr;
      delete []iflag;
      fin.close();
    }
  }
}

//------------------------------------------------------------------------------
void dmrg_su2::release_memory_all(int p,int lr){
//------------------------------------------------------------------------------
  int i,j;
  if(memory_flag==0)return;
  for(i=0;i<ns;i++){
    if(lr==0&&fll[p][i]==1||lr==1&&frr[p][i]==1)
      opr[p][i].clean();
    if(plqflg[p][i]){
      plqopr[p][i].clean();
      plqopr[p][i+ns].clean();
    }
  }
  hh[p].clean();
  uu[p].clean();
  for(i=0;i<exci;i++){
    ovlp[i][p].clean();
    orth[i][p].clean();
  }
}

//------------------------------------------------------------------------------
void dmrg_su2::release_memory2(int p,int lr){
//------------------------------------------------------------------------------
  int i,j,l,m;
  char base[100],len[10],dim[10],sec[10],pos[10],exc[10],name[100],pos2[10],id[10],dk[10];

  for(i=0;i<ns;i++){
    if(lr==0&&fll[p][i]==1||lr==1&&frr[p][i]==1)
      opr[p][i].clean();
    if(plqflg[p][i]){
      plqopr[p][i].clean();
      plqopr[p][i+ns].clean();
    }
  }
  hh[p].clean();
  if(memory_flag==0)return;
  sprintf(len,"%d",ns);
  sprintf(sec,"%d",totspin);
  sprintf(dim,"%d",max_dcut);
  sprintf(exc,"%d",exci);
  sprintf(pos,"%d",p);
  sprintf(id,"%d",jobid);
  sprintf(dk,"%d",comm_rank+3);

  for(l=0;l<3;l++){
    strcpy(base,"rm -f ./data");
    strcat(base,dk);
    //strcat(base,"/");
    //strcat(base,id);
    if(l==0) strcat(base,"/opr-");
    else strcat(base,"/plqopr-");
    strcat(base,len);
    strcat(base,"-");
    strcat(base,dim);
    strcat(base,"-");
    strcat(base,sec);
    strcat(base,"-");
    strcat(base,exc);
    strcat(base,"-");
    for(i=0;i<ns;i++){
      if(l==0&&(lr==0&&fll[p][i]==1||lr==1&&frr[p][i]==1||i==p))
	sprintf(pos2,"%d",i);
      else if(l>0&&plqflg[p][i]!=0&&(lr==0&&p>=plqpos[i][1]||lr==1&&p<=plqpos[i][2]))
	sprintf(pos2,"%d",i+(l-1)*ns);
      else continue;
      strcpy(name,base);
      strcat(name,pos);
      strcat(name,"-");
      strcat(name,pos2);
      strcat(name,".dat");
      system(name);
      strcpy(name,base);
      strcat(name,pos);
      strcat(name,"-");
      strcat(name,pos2);
      strcat(name,".bin");
      system(name);
    }
  }
}

//------------------------------------------------------------------------------
void dmrg_su2::wavefunc_transformation(int i,int flag){
//------------------------------------------------------------------------------
  int j;
  if(flag==0){
    uu[i].right2left_vectran();
    for(j=0;j<max_exci;j++)
      orth[j][i].right2left_vectran();
  }
  else if(flag==1){
    uu[i].left2right_vectran();
    for(j=0;j<max_exci;j++)
      orth[j][i].left2right_vectran();
  }
}
