#include <armadillo>
using namespace std;
using namespace arma;

mat EdgeMarker(cx_mat U, double flux, vec V, string filename){
  complex<double> ii(0,1);
  double pi = datum::pi;
  mat xy;
  xy.load("xy_pos.dat");
  vec idx;
  idx.load(filename);
  cx_mat localiz(xy.n_elem/2, 3);
  cx_double temp;
  for(int i = 0; i < U.n_cols; i++){
    localiz(i, 0) = flux;
    localiz(i, 1) = V(i);
    temp = 0;
    for (int j = 0; j < idx.n_elem; j++){
      temp += U(idx(j),i)*conj(U(idx(j),i));
    }
    localiz(i, 2) = temp;
  }
  mat real_local = conv_to<mat>::from(localiz);
  return real_local;
}
