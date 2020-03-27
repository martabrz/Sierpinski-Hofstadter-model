#include <armadillo>
using namespace std;
using namespace arma;

/* Input:
  states (complex mat) - eigenstates
  Ef (double) - the Fermi level

Output:
  bott (double) - Bott index for given occupied states */

mat BottIndex(cx_mat states, vec V, int Ef, double flux){
  complex<double> ii(0,1);
  double pi = datum::pi;
  mat xy;
  mat bott(1, 3);
  xy.load("xy_position.dat");
  vec x = xy.col(0);

  vec y = xy.col(1);
  int nstates = Ef;
  vec x_temp = unique(x);
  uvec c = hist(x, x_temp);
  mat Ny(x_temp.n_rows, 2);
  Ny.col(0) = x_temp;
  Ny.col(1) = conv_to<vec>::from(c);
  vec y_temp = unique(y);
  c = hist(y, y_temp);
  mat Nx(y_temp.n_rows, 2);
  Nx.col(0) = y_temp;
  Nx.col(1) = conv_to<vec>::from(c);
  int counter_x = 0;
  int counter_y = 0;
  cx_mat Px; Px.zeros(x.n_elem, x.n_elem);
  cx_mat Py; Py.zeros(y.n_elem, y.n_elem);
  for(int j = 0; j < (Ny.col(0)).n_elem; j++){
    for(int i = 0; i < Nx(j,1); i++){
      Px(i+counter_x, i+counter_x) = exp(x(i+counter_x)*2*pi*ii/double(Nx(j,1)));
    }
    counter_x += Nx(j,1);
    for(int i = 0; i < Ny(j,1); i++){
      Py(i+counter_y, i+counter_y) = exp(y(i+counter_y)*2*pi*ii/double(Ny(j,1)));
      }
    counter_y += Ny(j,1);
  }

  cout << "state: " << nstates << endl;
  cx_mat P_temp; P_temp.zeros(nstates, x.n_elem);

  cx_mat P; P.zeros(x.n_elem, nstates);
  P = (states.cols(0, nstates))*(states.cols(0, nstates)).t();
  cx_mat Px_proj(x.n_elem, x.n_elem);
  Px_proj = P*Px*P;
  cx_mat Py_proj(x.n_elem, x.n_elem);
  Py_proj = P*Py*P;
  cx_mat UVUV(x.n_elem, x.n_elem);
  UVUV = Px_proj*Py_proj*(Px_proj.t())*(Py_proj.t());
  cx_vec UVUV_eigval; cx_mat UVUV_eigvec;
  eig_gen(UVUV_eigval, UVUV_eigvec, UVUV);
  complex<long double> bbott;
  bbott = 0;
  for(int i = 0; i < UVUV_eigval.n_elem; i++){
    if (abs(UVUV_eigval(i)) > 1e-8){
    bbott += log(UVUV_eigval(i));
    }
  }
   bott(0) = imag(bbott)/(2*pi);
   bott(1) = V(nstates);
   bott(2) = flux;
   return bott;

}
