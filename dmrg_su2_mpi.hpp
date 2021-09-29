#ifndef DMRG_SU2_H
#define DMRG_SU2_H

#include "lanczos_su2.hpp"

#include "lut/LookupTable.hpp"

#include "typedefs.h"


using namespace std;

class dmrg_su2{
private:
  int ly,lx,bondd,ph,ns,read,read2,exci,max_exci,totspin,xleft,xrght,badposition,phdim,nfree,nbb,minlan,physpn;
  int **fll,**frr,**hmap,**plqpos,**plqflg,*qfll,*qfrr,**plqind;
  tensor_su2 *uu,*vv,*hh,*qq,*rr,***mm,*hh_bond,hhleft,hhrght,**opr,**plqopr,**ovlp,**orth,sigma[4],ring[10],qterm[10],*overlapvec,*plqplq;
  double gs_enr[100],enr[100],err,**ww,*wtmp,*bond_enr,*bond_enr_all,plqcorrelation[100];
  //this is to store the cgc coefficient

  void initialize_local_operators();
  void prepare_translationysquare_operators();
  void lanczos_solve_eigenvector(int,int,tensor_su2&);
  void lanczos_solve_eigenvector_idmrg(int,int,tensor_su2&);
  void transport_center_two_sites();
  void sweep_to_left();
  void sweep_to_right();
  void meet_again_at_center(int);
  void prepare_site_operator_from_left(int,int);
  void prepare_site_operator_from_right(int,int);
  void prepare_site_operator_from_left(int);
  void prepare_site_operator_from_right(int);
  void send_tensor_su2(tensor_su2&,int);
  void recv_tensor_su2(tensor_su2&,int);
public:
  dmrg_su2(int,int,int,int,int,int);
inline  dmrg_su2(){};
  ~dmrg_su2();
  double get_enr(){return gs_enr[exci];}
  void do_idmrg();
  void build_uu();
  void prepare_sweep();
  void prepare_measurement();
  void test_mpi();
  void synchronize_site(int,int,int,int);
  void sweep();
  void hamiltonian_vector_multiplication_idmrg(int,int,tensor_su2&,tensor_su2&);
  void save_mps2();
  void save_mps1();
  void save_enr();
  void read_enr(int,int);
  bool read_mps(int,int);
  bool read_mps_exci(int,int);
  bool read_mps_glide(int,int);
  void save_ww();
  void read_ww(int,int);
  void read_ww_exci(int,int);
  void read_ww_glide(int,int);
  void mea_enr();
  void mea_enr_part2();
  void mea_stgm();
  void save_memory(int,int);
  void read_memory(int,int);
  void read_orth(int);
  void release_memory(int,int);
  void release_memory_all(int,int);
  void release_memory2(int,int);
  void prepare_input_vector(int,int,tensor_su2&);
  void prepare_input_vector_initial(int,int,tensor_su2&);
  void makeup_clebsch_gordan_coefficient_tensors(int);
  void test_clebsch_gordan_coefficient();
  void check_mps_cgc();
  void wavefunc_transformation(int,int);
  void test_su2bond();
  
  // 5D sparse lookup tables
  LookupTable_5 fac_permutation_left, fac_permutation_rght;
	LookupTable_5 fac_operator_transformation_left, fac_operator_transformation_rght ;

  // 4D sparse lookup tables
  LookupTable_4 fac_operator_onsite_left ;
  LookupTable_4 fac_operator_onsite_rght ;
  LookupTable_4 fac_operator_pairup_left ;
  LookupTable_4 fac_operator_pairup_rght ;
  
  // 3D dense lookup table
  VecDbl3 fac_hamilt_vec ;

};


#endif
