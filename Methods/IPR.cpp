#include <armadillo>
using namespace arma;

mat InvParticipationRatio(cx_mat states, vec V, double flux){
  cx_mat IPR(states.n_cols, 3);
  mat IPR_real(states.n_cols, 3);
  IPR.zeros();
  cx_double frac_upp, frac_low;

  for(uword i = 0; i < states.n_cols; i++){
    frac_upp = 0; frac_low = 0;
    for(uword j = 0; j < states.n_rows; j++){
      frac_upp += pow(abs(states(j, i)), 4.0);
      frac_low += pow(abs(states(j, i)), 2.0);
    }
    IPR(i, 0) = frac_upp/pow(frac_low, 2);
    IPR(i, 1) = V(i);
    IPR(i, 2) = flux;
  }
  IPR_real = conv_to<mat>::from(IPR);
  return IPR_real;
}
