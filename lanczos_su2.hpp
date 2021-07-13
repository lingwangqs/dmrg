#include "tensor_su2.hpp"
class lanczos_su2{
private:
  int nrep,mlanc;
  double *aal,*nnl;
  double *eig,*vec;
  tensor_su2 *ff;

  void diatridiag(int);
public:
  lanczos_su2();
  ~lanczos_su2();
  void diag_op(int,int,tensor_su2&,tensor_su2&);
  void initialize_lanczos(tensor_su2&,int);
  void lanczos1(int,int,tensor_su2&,int);
  void compute_eigenvector(int,int,tensor_su2&,int,int&,int&);
  double get_eigval(){return eig[0];}
  void check_eigenvector(int,int,tensor_su2&,double&,double&);
};
