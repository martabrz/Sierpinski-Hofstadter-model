//Level distribution for given eigenvalues
#include <armadillo>

using std::cout;
using namespace arma;

mat LevelStatistics(vec V, vec epsilon, int k){
  mat levels((2*k + 1)*epsilon.n_elem, 2);
  vec temp_s(2*k + 1);
  double spacingsquared;
  double meanspacing;
  int idx;
  int a;
  int counter = 0;

  for(uword i = 0; i < epsilon.n_elem; i++){
    idx = index_min(abs(V - epsilon(i)));
    // cout << idx << endl;
      if(V(idx) > epsilon(i)){
          for(int j = -k; j <= k; j++){
            if(idx + j - 1 < 0 || idx + j > V.n_elem - 1){
              temp_s(j + k) = 100;
            }

            else{
              temp_s(j + k) = V(idx + j) - V(idx + j - 1);
            }
            levels(counter, 0)  = epsilon(i);
            counter++;

          }

          meanspacing = sum(temp_s)/double(temp_s.n_elem);
          temp_s /= meanspacing;

          a = 0;
          for(int m = i*(2*k +1); m < (i+1)*temp_s.n_elem; m++){
            levels(m, 1) = temp_s(a);
            a++;
          }
        }
// ////////////////////////////////////////////////////////////////
      else{
        for(int j = -k; j <= k; j++){
          if(idx + j < 0 || idx + j + 1 > V.n_elem - 1){
            temp_s(j + k) = 100;
          }

          else{
            temp_s(j + k) = V(idx + j + 1) - V(idx + j);
          }
          levels(counter, 0)  = epsilon(i);
          counter++;
        }

        meanspacing = sum(temp_s)/double(temp_s.n_elem);
        temp_s /= meanspacing;

        a = 0;
        for(int m = i*(2*k + 1); m < (i + 1)*temp_s.n_elem; m++){
          levels(m, 1) = temp_s(a);
          a++;
        }
      }
  }
  return levels;
}


mat VarianceDisorder(mat levels, int dis_aver, int k, int eps_elem, double flux){
   levels = levels.rows(stable_sort_index(levels.col(0)));
   int range = (2*k + 1)*dis_aver;
   mat var(eps_elem, 3);
   double mean_s2;
   double mean_s;

   for(int i = 0; i < eps_elem; i++){
     mean_s2 = 0;
     mean_s = 0;
     var(i, 0) = levels(i*range, 0);
     var(i, 2) = flux;
     for(int j = i*range; j < (i+1)*range; j++)
     {
       mean_s2 += pow(levels(j, 1), 2);
       mean_s += levels(j, 1);
     }
     mean_s2 /= range;
     mean_s /= range;
     var(i, 1) = mean_s2 - pow(mean_s, 2);
   }
  return var;
}
