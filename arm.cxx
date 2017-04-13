// sudo apt-get install libarmadillo-dev libarmadillo6 libarmadillo6-dbgsym
// g++ arm.cxx -std=c++11 -larmadillo --coverage -Wall -o arm && ./arm

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main() {

  mat A = randu<mat>(100,2);

  for (int i = 0; i < 100; i++) {
    A(i,1) = 0.5*A(i,0) + 0.5*A(i,1);
  }

  // cout << A << endl;

  mat coeff;
  mat score;
  vec latent;
  vec tsquared;

  princomp(coeff, score, latent, tsquared, A);

  cout << coeff  << endl;
  cout << latent << endl;
  //cout << score  << endl;

  return 0;
}
