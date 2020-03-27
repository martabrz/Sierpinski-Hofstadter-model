#include <armadillo>
using namespace std;
using namespace arma;

mat LocalizationSpatial(cx_mat U, int state, double eval){
  complex<double> ii(0,1);
  double pi = datum::pi;
  mat xy;
  xy.load("xy_position.dat");
  cx_mat localiz(xy.n_elem/2, 4);
  for(int i = 0; i < xy.n_elem/2; i++){
    localiz(i, 0) = xy(i, 0);
    localiz(i, 1) = xy(i, 1);
    localiz(i, 2) = U(i, state)*conj(U(i, state));
    localiz(i, 3) = eval;
  }
  mat real_local = conv_to<mat>::from(localiz);
   return real_local;
}
