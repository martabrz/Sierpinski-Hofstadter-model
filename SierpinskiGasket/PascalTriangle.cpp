#include <armadillo>

using namespace arma;

int cal(int pow, int val, int MOD){
  if(pow == 0)
    return 1;
  int v = cal(pow/2, val, MOD);
    if(pow % 2 == 0)
      return (v*v) % MOD;
    else
      return (((v*val) % MOD)*v) % MOD;
}

umat CreateTriangle(int n, int m){
  int r = pow(2, n);
  // int r = 81;
  umat lattice; lattice.zeros(r, r);
  umat triangle; triangle.zeros(r, r);

  // create normal triangle
  for(int j = 0; j < r; j++){
    for(int i = 0; i <= j; i++){
        triangle(j, i) = 1;
    }
  }

  // left lower corner
  lattice(r - 1, 0) = 1;
  lattice(r - 2, 0) = 1; lattice(r - 1, 1) = 1;
  if(m > 1){
    // Pascal's triangle
    for(int i = r - 1; i >= 0; i--){
      for(int j = 0; j <= i; j++){
        if((i == r - 1) && j != 0 && j != 1){
          lattice(i, j) = lattice(i, j - 1);
        }
        if(j == 0 && (i != r - 1) && (i != r - 2)){
          lattice(i, j) = lattice(i + 1, j);
        }
        if((i != r - 1) && j != 0)
        {
          lattice(i, j) = lattice(i, j - 1) + lattice(i + 1, j);
        }
      }
    }
    // create lattice modulus m
    for(int j = 0; j < r; j++){
      for(int i = 0; i < r; i++){
        lattice(j, i) %= m;
        // cout << lattice(j,i) << endl;
        if(lattice(j, i) > 1){
          lattice(j, i) = 1;
        }
      }
    }
    return lattice;
  }

  if(m == 1){
    throw std::invalid_argument("m = 1: empty lattice");
  }
  if(m == 0){
    std::cout << "m = 0: creating triangular lattice" << endl;
    return triangle;
  }
  if(m < 0){
    throw std::invalid_argument("m < 0: invalid argument");
  }
}

umat CreateTrianglePBC(int n, int m){
  int r = pow(2, n) + 1;
  umat lattice; lattice.zeros(r, r);
  umat triangle; triangle.zeros(r, r);
  umat lattice_all; lattice_all.zeros(r+1, r+1);
  // create normal triangle
  for(int j = 0; j < r; j++){
    for(int i = 0; i < r; i++){
        triangle(j, i) = 1;
    }
  }

  // left lower corner
  lattice(r - 1, 0) = 1;
  lattice(r - 2, 0) = 1; lattice(r - 1, 1) = 1;
  if(m > 1){
    // Pascal's triangle
    for(int i = r - 1; i >= 0; i--){
      for(int j = 0; j <= i; j++){
        if((i == r - 1) && j != 0 && j != 1){
          lattice(i, j) = lattice(i, j - 1);
        }
        if(j == 0 && (i != r - 1) && (i != r - 2)){
          lattice(i, j) = lattice(i + 1, j);
        }
        if((i != r - 1) && j != 0)
        {
          lattice(i, j) = lattice(i, j - 1) + lattice(i + 1, j);
        }
      }
    }



    // create lattice modulus m
    for(int j = 0; j < r - 1; j++){
      for(int i = 0; i <= j; i++){
        lattice(j, i) %= m;
        if(lattice(j, i) > 1){
          lattice(j, i) = 1;
        }
      }
    }
    // lattice_all = lattice;
    // lattice_all = lattice + flipud(fliplr(lattice));
    // for(int i = 0; i < r; i++){
    //   lattice_all(i, i) = 1;
    // }
    for(int i = 1; i <= r; i++){
        for(int j = 0; j < r; j++){
        lattice_all(i,j) = lattice(i-1,j);
      }
    }
    // mirror copy
    for(int i = 0; i < r; i++){
        for(int j = r; j > i ; j--){
        lattice_all(i,j) = lattice(r-1-i, r - j);
        // lattice_all(i,j) = 0;
      }
    }
    return lattice_all;
  }

  if(m == 1){
    throw std::invalid_argument("m = 1: empty lattice");
  }
  if(m == 0){
    std::cout << "m = 0: creating triangular lattice" << endl;
    return triangle;
  }
  if(m < 0){
    throw std::invalid_argument("m < 0: invalid argument");
  }
}
