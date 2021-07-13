#include "tensor.hpp"
#include "su2struct.hpp"

using namespace std;

class tensor_su2{
private:
  int nbond,nten,locspin;
  tensor *tarr,**parr;
  tensor *tcgc,**pcgc;
  su2struct cgc;
public:
  tensor_su2();
  ~tensor_su2();
  void set_tensor_su2(int,int,su2bond*,tensor**);
  void set_tensor_su2(su2struct&,double*,double*,int*);
  void set_tensor_su2(su2struct&,double*,double*);
  void set_tensor_su2(su2struct&,double*);
  inline int get_nbond(){return nbond;}
  inline int get_nten(){return nten;}
  inline int get_locspin(){return locspin;}
  void clean();
  void make_spinor_start(int);
  void make_spinor_end(int);
  void make_spinhalf_permutation_rightmove(tensor_su2&);
  void make_spinhalf_permutation_leftmove();
  void make_spinhalf_permutation_rightmove_4legs(tensor_su2&);
  void make_spinhalf_permutation_leftmove_4legs();
  void make_spinhalf_permutation_rightend();
  void make_spinhalf_permutation_rightend_v2();
  void make_spinhalf_permutation_leftend();
  void make_spinhalf_on_vec_left(tensor_su2&,tensor_su2&,int);
  void make_spinhalf_on_vec_right(tensor_su2&,tensor_su2&,int);
  void make_spinhalf_Qterm(tensor_su2&,tensor_su2&,tensor_su2&,tensor_su2&);
  void contract(tensor_su2&,int,tensor_su2&,int);
  void conjugate(int);
  void take_conjugate(int);
  void take_conjugate();
  void make_standard_cgc();
  void fuse(int,int);
  void fuse(int,int,int,int);
  void fuse(su2bond&,su2bond&);
  void fuse(su2bond&,su2bond&,int);
  void fuse_to_multiplet(su2bond&,su2bond&,int);
  void fuse_to_multiplet(int,int,int);
  void fuse_to_singlet(su2bond&,su2bond&);
  void operator_tensor_product_identity(tensor_su2&,su2bond&);
  void operator_tensor_product_identity(tensor_su2&,tensor_su2&);
  void print();
  bool is_null(int);
  bool is_null();
  bool check_angularmoments(int*);
  void shift(int,int);
  void exchangeindex(int,int);
  void left2right_vectran();
  void right2left_vectran();
  void reflection_vectran();
  int svd(tensor_su2&,double,tensor_su2&,double,double*);
  void multiply_singular_value(int,double*);
  void devide_singular_value(int,double*);
  void get_su2bond(int,su2bond&);
  int get_angularmoment_index(int,int);
  int get_tensor_index(int*);
  bool get_tensor_argument(int,int*,int*,int*);
  int get_nmoment(int);
  int get_bonddir(int);
  int get_angularmoment(int,int);
  int get_bonddim(int,int);
  int get_cgcdim(int,int);
  tensor_su2& operator = (double); 
  tensor_su2& operator = (tensor_su2&); 
  tensor_su2& operator += (tensor_su2&);
  tensor_su2& operator -= (tensor_su2&);
  tensor_su2& operator *= (double);
  tensor_su2& operator /= (double);
  bool operator == (tensor_su2&); 
  bool operator != (tensor_su2&); 
  double inner_prod(tensor_su2&);
  double take_trace();
  double ss_inner_prod(tensor_su2&);
  tensor* get_parr(int);
  tensor* get_pcgc(int);
  su2struct& get_cgc();
  void get_nelement(int&,int&);
  void get_telement(double*,double*);  
  void get_telement(double*,double*,int*);  
  void cgc_make_scalar_operator();
  void cgc_make_cgc();
  void cgc_make_identity_cgc();
  void cgc_make_wigner_eckart_theorem(int);
  double normalize_vector();
  void diagonalize();
  bool check_unitary_transformation_cgc(int);
  void makeup_input_vector();
  void initialize_input_vector();
  void direct_sum(int,tensor_su2&,tensor_su2&);
  void orthogonalize(tensor_su2*,int);
  void overlap_initial(tensor_su2&,tensor_su2&,int);
  void overlap_transformation(tensor_su2&,tensor_su2&,tensor_su2&,int);
  void operator_initial(tensor_su2&,tensor_su2&,tensor_su2&,int);
  void operator_initial(tensor_su2&,tensor_su2&,tensor_su2&,tensor_su2&,int);
  void operator_transformation(tensor_su2&,tensor_su2&,tensor_su2&,int);
  void operator_pairup(tensor_su2&,tensor_su2&,tensor_su2&,tensor_su2&,int);
  void hamiltonian_vector_multiplication(tensor_su2&,tensor_su2&,tensor_su2&);
  void permutation(tensor_su2&,tensor_su2&,tensor_su2&,tensor_su2&,tensor_su2&,int);
};
