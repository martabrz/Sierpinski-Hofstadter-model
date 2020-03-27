#include <armadillo>
#include <random>
#include <cmath>
using namespace std;
using namespace arma;

bool DoesFileExist(const string& name) {
  /* Checks if file 'name' exists */
    ifstream f(name.c_str());
    return f.good();
}

cx_mat CreateHamiltonian(umat lattice, double t, double flux){
  /* Builds a nearest-neighbour tight-binding Hamiltonian with OBC */
  complex<double> ii(0,1);
  int l_size; l_size = lattice.n_rows; int h_size; h_size = pow(l_size, 2);
  cx_mat hamiltonian; hamiltonian.zeros(h_size, h_size);
  uvec idx_ones = find(lattice == 1); // position of all lattice sites
  umat pos_ones = ind2sub(size(lattice), idx_ones).t();
  uvec idx = linspace<uvec>(0, lattice.n_elem - 1, lattice.n_elem);
  umat pos = ind2sub(size(lattice), idx).t();
  for(uword i = 0; i < idx.n_elem; i++){
    for(uword j = 0; j < idx.n_elem; j++){
      if(sqrt(pow((pos(i, 0) - pos(j, 0)),2) +
          pow((pos(i, 1) - pos(j, 1)), 2)) == 1){
            if(pos(i, 0) == pos(j, 0)){  // x direction
             hamiltonian(i, j) = t*exp(-ii*double(pos(i, 0)*flux));
             hamiltonian(j, i) = conj(hamiltonian(i, j));
            }
           if(pos(i, 1) == pos(j, 1)){  // y direction
             hamiltonian(i, j) = t;
             hamiltonian(j, i) = conj(hamiltonian(i, j));
           }
      }
    }
  }
  if (!DoesFileExist("xy_positions.dat")){
    pos_ones.swap_cols(0,1);
    pos_ones.save("xy_positions.dat", raw_ascii);
  }
  return hamiltonian.submat(idx_ones, idx_ones);
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
