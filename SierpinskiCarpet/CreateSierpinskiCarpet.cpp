/* Create Sierpinski carpet using the Hadamard matrix
Input:
 n (integer) - iteration of a fractal
 s (integer) - lattice size
Output:
 lattice (umat) - Sierpinski carpet as a lattice */

#include <armadillo>
using namespace std;
using namespace arma;

umat CreateHadamard(int n, int s){
  umat hadamard; hadamard.zeros(s, n);
  for (uword i = 0; i < s; i++){
    for(uword j = 0; j < n; j++){
      hadamard(i,j) = (( int(i/pow(3, j) ) % 3 ) % 2);
    }
  }
  return hadamard;
}

umat CreateSierpinskiCarpet(int n, int s, umat hadamard){
  umat lattice; lattice.zeros(s, s);
  int sum;
  for (uword i = 0; i < s; i++){
    for(uword j = 0; j < s; j++){
      sum = 0;
      for(uword k = 0; k < n; k++){
        sum += hadamard(i,k)*hadamard(j,k);
      }
      if(sum == 0){
        lattice(i,j) = 1;
      }
      else{
        lattice(i,j) = 0;
      }
    }
  }
  return lattice;
}
