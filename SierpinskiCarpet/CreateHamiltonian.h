arma::cx_mat CreateHamiltonian(arma::umat lattice, double t, double flux);
arma::cx_mat OnSiteDisorder(arma::cx_mat hamiltonian, double W);
bool DoesFileExist(const std::string& name);
