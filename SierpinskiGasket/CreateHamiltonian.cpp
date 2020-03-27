#include <random>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

bool DoesFileExist(const string& name) {
  /* Checks if file 'name' exists */
    return f.good();
}



cx_mat CreateHamiltonianPBC(umat lattice, double t, double flux){
  complex<double> ii(0,1);
  int l_size; l_size = lattice.n_rows; int h_size; h_size = pow(l_size, 2);
  cx_mat hamiltonian; hamiltonian.zeros(h_size, h_size);
  uvec idx = linspace<uvec>(0, lattice.n_elem - 1, lattice.n_elem);
  umat pos = ind2sub(size(lattice), idx).t();
  uvec idx_zeros = find(lattice == 0);
  umat pos_zeros = ind2sub(size(lattice), idx_zeros).t();
  double r;
  for(uword i = 0; i < idx.n_elem; i++){
    for(uword j = 0; j < idx.n_elem; j++){
      r = sqrt(pow((double(pos(i, 0)) - double(pos(j, 0))), 2) +
               pow((double(pos(i, 1)) - double(pos(j, 1))), 2));
      if(r > 1.35 && r < 1.45 && pos(i, 0) == pos(j, 0) - 1 &&
         pos(i, 1) == pos(j, 1)-1 && pos(i, 0) != l_size - 1){
        hamiltonian(i, j) = t*exp(-ii*flux*(-1.0*double(pos(i, 0)) - 0.5));
        hamiltonian(j, i) = conj(hamiltonian(i, j));
      }
      if(r > 0.95 && r < 1.05 && pos(i, 1) == pos(j, 1)){
        hamiltonian(i, j) = t;
        hamiltonian(j, i) = conj(hamiltonian(i, j));
      }
      if(r > 0.95 && r < 1.05 && pos(i, 0) == pos(j, 0) &&
         pos(i, 1) < pos(j, 1)){
        hamiltonian(i, j) = t*exp(ii*double(pos(i, 0))*flux);
        hamiltonian(j, i) = conj(hamiltonian(i, j));
      }

      if(sqrt(pow((pos(i, 0) - pos(j, 0)),2) +
          pow((pos(i, 1) - pos(j, 1)), 2)) == l_size - 1){
            if(pos(i, 0) == pos(j, 0)){
              hamiltonian(i, j) = 0*t*exp(ii*double(pos(i, 0))*flux);
              hamiltonian(j, i) = conj(hamiltonian(i, j));
            }
            if(pos(i, 1) == pos(j, 1)){
              hamiltonian(i, j) = 0*t*exp(-ii*double(pos(i, 1)*flux*l_size));
              hamiltonian(j, i) = conj(hamiltonian(i, j));
            }
          }
      if(r < sqrt(pow(l_size - 1, 2) + 1) + 0.1 &&
         r > sqrt(pow(l_size - 1, 2) + 1) - 0.1 ){
            if(pos(i, 0) == l_size - 1 && pos(j, 0) == 0 && pos(i, 1) < pos(j, 1)){
              hamiltonian(i, j) = 0*t*exp(-ii*(0.5 + (double(pos(i, 1)))*l_size)*flux);
              hamiltonian(j, i) = conj(hamiltonian(i, j));
            }
            if(pos(i, 1) == l_size - 1 && pos(j, 1) == 0 && pos(i, 0) < pos(j, 0)){
              hamiltonian(i, j) = 0*t*exp(-ii*flux*(-1.0*double(pos(i, 0)) - 0.5));
              hamiltonian(j, i) = conj(hamiltonian(i, j));
            }
          }
      if(pos(i, 0) == l_size - 1 && pos(j, 0) == 0 &&
         pos(i, 1) == l_size - 1 && pos(j, 1) == 0){
         hamiltonian(i, j) = 0*t*exp(-ii*(0.5 + (double(pos(i, 1)))*l_size)*flux);
         hamiltonian(j, i) = conj(hamiltonian(i, j));
       }
     }
   }

   if(idx_zeros.n_elem > 0){

   uvec erase(idx_zeros.n_elem); uvec e;
   for(uword k = 0; k < idx_zeros.n_elem; k++){
     e = find(pos_zeros(k, 0) == pos.col(0) && pos_zeros(k, 1) == pos.col(1));
     erase(k) = e(0);
   }

   uvec l_triangle(idx_zeros.n_elem/2 - l_size/2);
   uvec u_triangle(idx_zeros.n_elem/2 - l_size/2);
   int counter = 0;
   double x1, y1, x2, y2, x3, y3;
   double dot1, dot2, dot3;
   double x, y;
   x1 = 0; y1 = 1; x2 = 0; y2 = l_size - 1; x3 = l_size - 2; y3 = l_size -1;
   for(uword k = 0; k < idx_zeros.n_elem; k++){
     x = pos(erase(k), 1);
     y = pos(erase(k), 0);
     dot1 = (y2 - y1)*(x - x1) + (-x2 + x1)*(y - y1);
     dot2 = (y3 - y2)*(x - x2) + (-x3 + x2)*(y - y2);
     dot3 = (y1 - y3)*(x - x3) + (-x1 + x3)*(y - y3);
     if(dot1 >= 0 && dot2 >= 0 && dot3 >= 0){
       l_triangle(counter) = erase(k);
       counter++;
     }
   }
   counter = 0;
   x1 = l_size - 1; y1 = 0; x2 = 1; y2 = 0; x3 = l_size - 1; y3 = l_size -2 ;
   for(uword k = 0; k < idx_zeros.n_elem; k++){
     x = pos(erase(k), 1);
     y = pos(erase(k), 0);
     dot1 = (y2 - y1)*(x - x1) + (-x2 + x1)*(y - y1);
     dot2 = (y3 - y2)*(x - x2) + (-x3 + x2)*(y - y2);
     dot3 = (y1 - y3)*(x - x3) + (-x1 + x3)*(y - y3);
     if(dot1 >= 0 && dot2 >= 0 && dot3 >= 0){
       u_triangle(counter) = erase(k);
       counter++;
     }
   }
   for(uword k = 0; k < l_triangle.n_elem; k++){
     hamiltonian(l_triangle(k), l_triangle(k) - 1) = 0;
     hamiltonian(l_triangle(k) - 1, l_triangle(k)) = 0;
     hamiltonian(l_triangle(k), l_triangle(k) + l_size) = 0;
     hamiltonian(l_triangle(k) + l_size, l_triangle(k)) = 0;
     hamiltonian(l_triangle(k) - 1 , l_triangle(k) + l_size) = 0;
     hamiltonian(l_triangle(k) + l_size, l_triangle(k) - 1) = 0;
   }

   for(uword k = 0; k < u_triangle.n_elem; k++){
     hamiltonian(u_triangle(k), u_triangle(k) + 1) = 0;
     hamiltonian(u_triangle(k) + 1, u_triangle(k)) = 0;
     hamiltonian(u_triangle(k), u_triangle(k) - l_size) = 0;
     hamiltonian(u_triangle(k) - l_size, u_triangle(k)) = 0;
     hamiltonian(u_triangle(k) + 1 , u_triangle(k) - l_size) = 0;
     hamiltonian(u_triangle(k) - l_size, u_triangle(k) + 1) = 0;
   }

   uvec rem(3*idx_zeros.n_elem);
   counter = 0;
   for(uword k = 1; k < l_triangle.n_elem; k++){
     for(uword l = 0; l < l_triangle.n_elem - 1; l++){
       if((pos(l_triangle(k), 1) == pos(l_triangle(k-1), 1)) &&
          (pos(l_triangle(k), 0) == pos(l_triangle(l), 0) + 1) &&
          (pos(l_triangle(k), 1)  == pos(l_triangle(l), 1) -1) &&
          (pos(l_triangle(k), 0)  == pos(l_triangle(l), 0) +1) &&
          (pos(l_triangle(l), 1) == pos(l_triangle(l+1), 1)) &&
          (pos(l_triangle(l), 0) == pos(l_triangle(l+1), 0)- 1) &&
          (pos(l_triangle(k), 0) == pos(l_triangle(l+1), 0))
        ){
            rem(counter) = l_triangle(k-1);
            counter++;
            rem(counter) = l_triangle(l);
            counter++;
            rem(counter) = l_triangle(l+1);
            counter++;
        }
     }
   }
   for(uword k = 0; k < u_triangle.n_elem - 1; k++){
     for(uword l = 1; l < u_triangle.n_elem; l++){
       if(pos(u_triangle(k), 1) == pos(u_triangle(k+1), 1) &&
          pos(u_triangle(k), 0) == pos(u_triangle(k+1), 0) -1 &&
          pos(u_triangle(k), 0) == pos(u_triangle(l), 0) - 1 &&
          pos(u_triangle(k), 1) == pos(u_triangle(l), 1) + 1 &&
          pos(u_triangle(l), 1) == pos(u_triangle(l-1), 1)){
            rem(counter) = u_triangle(l-1);
            counter++;
            rem(counter) = u_triangle(l);
            counter++;
            rem(counter) = u_triangle(k+1);
            counter++;
        }
     }
   }
   if(counter < 3*idx_zeros.n_elem - 1){
     rem.shed_rows(counter, 3*idx_zeros.n_elem - 1);
   }
   counter = 0;
   if (rem.n_elem > 0){
   uvec rem_sorted = sort(rem);
   for(int i = 0 ; i < rem_sorted.n_elem - 1; i++){
    if(rem_sorted(i) == rem_sorted(i+1)){
      rem_sorted.shed_row(i);
      i--;
    }
   }
   for(int i = idx.n_elem - 1; i >=0; i--){
     for(int j = 0; j < rem_sorted.n_elem; j++){
       if(idx(i) == rem_sorted(j)){
         idx.shed_row(i);
       }
     }
   }
 }
 int c = 0;
 uvec rem_upper((l_size*l_size - l_size)/2);
 for(int i = 1; i < l_size; i++){
  for(int j = 0; j < i; j++){
   rem_upper(c) = i*l_size + j;
         c++;
       }
     }

     for(int i = idx.n_elem - 1; i >=0; i--){
       for(int j = 0; j < rem_upper.n_elem; j++){
         if(idx(i) == rem_upper(j)){
           idx.shed_row(i);
         }
       }
     }
  uvec diag_join(l_size);
  for(int i = 0; i < l_size; i++){
    diag_join(i) = i*l_size + i;
  }
// if two triangles are decoupled
  for(int i = idx.n_elem - 1; i >=0; i--){
    for(int j = 0; j < diag_join.n_elem; j++){
      if(idx(i) == diag_join(j)){
        idx.shed_row(i);
      }
    }
  }
  if (!DoesFileExist("xy_pos.dat")){
    umat latt(idx.n_elem, 2);
    latt = pos.rows(idx);
    latt.swap_cols(0,1);
    latt.save("xy_pos.dat", raw_ascii);
  }

   return hamiltonian(idx, idx);
 }
 else{
  return hamiltonian;
 }
}

cx_mat OnSiteDisorder(cx_mat hamiltonian, double W){
  /* Generates random on-site potential drawn from [-W/2, W/2] based on 32-bit Mersenne Twister generator */
  double seed = chrono::high_resolution_clock::now().time_since_epoch().count();
  mt19937 generator(seed);
  uniform_real_distribution<double> dis(-W/2, W/2);

  for(int i = 0; i < hamiltonian.n_cols; i++){
      hamiltonian(i, i) = dis(generator);
    }
  return hamiltonian;
}
