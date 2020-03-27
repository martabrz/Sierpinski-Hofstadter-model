#include <armadillo>
using namespace std;
using namespace arma;
/* Calculates the real-space Chern number
Input:
  states (complex mat) - eigenstates
  V (double vec) - eigenvalues
  Ef (double) - Fermi level
Output:
  chern (double) - local Chern number for a given set of occupied states */

mat LocalChern(cx_mat states, vec V, double Ef){
  complex<double> ii(0,1); double pi = datum::pi;
  mat xy; xy.load("xy_positions.dat");
  vec x = xy.col(0); vec y = xy.col(1);
  uvec idx = linspace<uvec>(0, x.n_elem - 1, x.n_elem - 1);
  uvec occ = find(V <= Ef);
  int countA, countB, countC;
  double phi, tx, ty;
  double x1, y1, x2, y2, x3, y3, x4, y4;
  double x10, y10, x20, y20, x30, y30, x40, y40;
  double px, py, p21sq, p41sq;
  int counter;
  double middle_x, middle_y;
  cx_mat P; complex<double> P_temp;
  int nstates = occ.n_elem - 1;
  colvec vertices0(8);
  colvec verticesS(8);
  mat scaling = zeros(8, 8);
  vec s = linspace<vec>(0.01, 2, 20);
  mat chern(s.n_elem, 3);
  // define a rectangle as in the article
  x10 = 26-8; y10 = 62; x20 = 26-8; y20 = 45; x30 = 0+8; y30 = 45; x40 = 0+8; y40 = 62;
  middle_x = (x40+x20)/2; //max(x)/2.0;
  middle_y = (y40+y20)/2; //max(y)/2.0;
  vertices0 << x10 << endr << x20 << endr << x30 << endr << x40 << endr
            << y10 << endr << y20 << endr << y30 << endr << y40 << endr;
  for(uword z = 0; z < s.n_elem; z++){
    uvec sites(x.n_elem);
    uvec A(sites.n_elem);
    uvec B(sites.n_elem);
    uvec C(sites.n_elem);
    scaling = s(z)*eye(8, 8);
    verticesS = scaling*vertices0; // clockwise convention
    x1 = verticesS(0); x2 = verticesS(1); x3 = verticesS(2); x4 = verticesS(3);
    y1 = verticesS(4); y2 = verticesS(5); y3 = verticesS(6); y4 = verticesS(7);
    counter = 0;
    p21sq = pow(x2-x1, 2) + pow(y2-y1, 2);
    p41sq = pow(x4-x1, 2) + pow(y4-y1, 2);
    for(uword k = 0; k < x.n_elem; k++){
      px = x(k) - x1;
      py = y(k) - y1;
      if ((px*(x2-x1) + py*(y2-y1) >=0) && px*(x2-x1) + py*(y2-y1) <= p21sq){
        if ((px*(x4-x1) + py*(y4-y1) >=0) && px*(x4-x1) + py*(y4-y1) <= p41sq){
          sites(counter) = k;
          counter++;
        }
      }
     }
    sites = sites.head(counter);
    countA = 0; countB = 0; countC = 0;
    for(uword i = 0; i < sites.n_elem; i++){
      ty = y(sites(i)) - middle_y;
      tx = x(sites(i)) - middle_x;
      phi = atan2(-ty, tx);
      if(phi < 0){
        phi += 2*pi;
      }
      if ((phi >= (11.0/6.0)*pi && phi < 2*pi) || (phi >= 0*pi && phi < 0.5*pi)){
        A(countA) = sites(i);
        countA++;
      }
      if (phi >= pi*0.5 && phi < (7.0/6.0)*pi){
        B(countB) = sites(i);
        countB++;
      }
      if (phi >= (7.0/6.0)*pi && phi < (11.0/6.0)*pi){
        C(countC) = sites(i);
        countC++;
      }
    }
    A = A.head(countA); B = B.head(countB); C = C.head(countC);
    P = (states.cols(0, nstates))*(states.cols(0, nstates)).t();
    P_temp = 0;
    for(uword j = 0; j < A.n_elem; j++){
      for(uword k = 0; k < B.n_elem; k++){
        for(uword l = 0; l < C.n_elem; l++){
        P_temp += P(A(j), B(k))*P(B(k), C(l))*P(C(l), A(j)) - P(A(j), C(l))*P(C(l), B(k))*P(B(k), A(j));
        }
      }
    }
   P_temp *= 12*pi;
   chern(z, 0) = imag(P_temp);
   chern(z, 1) = sites.n_elem;
   chern(z, 2) = s(z);
  }
 return chern;
}
