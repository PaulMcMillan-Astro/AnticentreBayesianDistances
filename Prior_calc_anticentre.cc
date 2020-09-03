/*-----------------------------------------------------------------------------
*
* Distance_calc_anticentre.cc
* Author: PJM
* email: paul@astro.lu.se
* Github: https://github.com/PaulMcMillan-Astro/
*
* Description: Code which calculates distances based on a Bayesian scheme, a la
* McMillan 2018
*
*-----------------------------------------------------------------------------*/



#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>


#include "Pi.h"
#include "Prior.h"
#include "Gauss-Legendre-header.h"

using namespace std;

// Value of a Gaussian, given value, mean, sigma
template <class T>
   T gaussian(T x, T x0, T sig) {
   T isig = 1./sig, tmp = (x-x0)*isig, tmp2 = -0.5*tmp*tmp;
   const T isrttpi = 1./ 2.506628274631000502415765284811;
   return isrttpi*isig*exp(tmp2);
 }

// Value of a Gaussian, given value, mean, sigma, unnormalised
// Only use when making comparisons!
 template <class T>
   T gaussian_unnorm(T x, T x0, T sig) {
   T isig = 1./sig, tmp = (x-x0)*isig, tmp2 = -0.5*tmp*tmp;
   return exp(tmp2);
 }


const double deg2rad  = Pi/180.;


/*

Main programme

*/
int main(int argc, char *argv[] ) {
// Prior chosen - can alter by hand
  Basic_prior Prior;
  double l_min_deg,l_max_deg;
  ofstream to;
  if(argc==4) {
    to.open(argv[1]);
    l_min_deg = atof(argv[2]);
    l_max_deg = atof(argv[3]);
  } else {
    cerr << "Usage: "<<argv[0]<<"outputname l_min_deg l_max_deg\n";
    return 1;
  }
  if (l_min_deg>=l_max_deg) {
    cerr << "l_min must be smaller than l_max";
    return 1;
  }
  // input
  double l_min=l_min_deg*deg2rad,l_max=l_max_deg*deg2rad,b_max=10*deg2rad,r_min=100, r_max=9900;
  int nr=50, nl=50, nb=1000;
  double integral[nr];  

  to << "distance,integrated_prior\n";

  for(int i=0;i<nr;i++) {
      integral[i] = 0.;
  }
      
    
  for(int j=0;j<nl;j++){
    double weight_l = 1.;
    if(j==0 || j== nl-1) weight_l = 0.5;
    double l = l_min+j*(l_max-l_min)/(nl-1);
    for(int k=0;k<nb;k++) {
        double weight_b = 1.;
        if(k==0 || k== nb-1) weight_b = 0.5;
        double b = k*(b_max)/(nb-1);
        Prior.prepare_sightline(l, b);
        for(int i=0;i<nr;i++) {
            double r = r_min+i*(r_max-r_min)/(nr-1);
            integral[i] += weight_b*weight_l*Prior.prior(r)*Prior.cosb;
        }
    }
  }
    for(int i=0;i<nr;i++)   {
        to << 1e-3*(r_min+i*(r_max-r_min)/(nr-1)) << ',' << integral[i]<<'\n';
    }
    return 0;
}
