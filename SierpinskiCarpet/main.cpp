#include <iostream>
#include <fstream>
#include <armadillo>
#include "CreateHamiltonian.h"
#include "CreateSierpinskiCarpet.h"
#include "LocalChern.h"
using std::cout;
using std::ofstream;
using namespace arma;

/* main code for calculating the spectrum and the Chern number */

int main(int argc, char* argv[]) {
    wall_clock timer;
    timer.tic();
    double pi = datum::pi;
    ofstream chern("chern.dat"); chern.precision(10); chern.setf(ios::fixed);
    int n = 4; int s = pow(3, n); // iteration
    umat lat; lat = CreateSierpinskiCarpet(n, s, CreateHadamard(n, s));
    double flux = 2*pi*0.25;
    vec V; cx_mat U; cx_mat H;
    H = CreateHamiltonian(lat, -1, flux);
    cout << "Creating Hamiltonian in (sec): " << timer.toc() << endl;
    timer.tic();
    eig_sym(V, U, H, "dc").raw_print;
    cout << "Diagonalizing Hamiltonian in (sec): " << timer.toc() << endl;
    LocalChern(U, V, -1.5).raw_print(chern);
    timer.tic();
    cout << "Computing local Chern number in (sec): " << timer.toc() << endl;
    return 0;
  }
