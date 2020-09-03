/*-----------------------------------------------------------------------------
*
* Prior.cc
* Author: PJM
* email: paul@astro.lu.se
* Github: https://github.com/PaulMcMillan-Astro/
*
* Description: Priors used in Bayesian distance determination.
*    Specifically here: For distances derived from parallaxes given in Gaia's
*    DR2
*-----------------------------------------------------------------------------*/

#include <cmath>
#include <iostream>
#include "Prior.h"



/*-----------------------------------------------------------------------------
*
* Class: Basic_prior
*
* See Prior.h
*
*-----------------------------------------------------------------------------*/

/*---------  Basic_prior::prepare_sightline: Set l,b etc ------------*/

void Basic_prior::prepare_sightline(double L,double B) {
  l=L; b=B;
  cosl=cos(l); sinl=sin(l);
  cosb=cos(b); sinb=sin(b);
}

/* ---------- setup: set Galactic position - used by prior()  -------*/
void Basic_prior::setup(double s) {
  // the below are slighty wrong because the Galactic coordinate
  // system isn't aligned with the plane of the disc.
  R =  sqrt(R0*R0 + s*s*cosb*cosb - 2.*R0*s*cosb*cosl);
  z = s*sinb+z0;
  r = sqrt(R*R+z*z);
  offgrid = ( s<0.)? true : false;
}

/* ---------- prior(s) : returns prior given distance (in pc) -----*/
double Basic_prior::prior(double s) {
  setup(s);
  if(offgrid) return 0.;  // if s<0
  return (prior_thin()+prior_thick()+
	  prior_halo()+prior_bulge())* //position
    s*s;                         //r^2 (volume between s & s+ds)
}


/*----------------     Constructor  -------------------*/
Basic_prior::Basic_prior() :
  R0(8330.), z0(0.), iRd_thin(1./2600.), izd_thin(1./300.),
 norm_thick(0.0596),
  iRd_thick(1./3600.), izd_thick(1./900.),
  norm_halo(3.774e9),
  power_halo(-3.39), // see Mcmillan et al 2018 (or Burnett & Binney 2010)
  rcut_halo(3000.),  // cutoff otherwise dominant
  qb(0.5), gammab(0.),betab(1.8),
  r0b(0.075*1000.), rcutb( 2.1*1000.) // Bulge - see McMillan 2017
{

  // Setup bulge central density such that ratio to disc mass is correct
  // Some 'magic' number from know properties of bulge model. Done
  // arithmetically to stop me screwing it up.
  double Md = 8*3.141592*((1./izd_thin)*powf(1./iRd_thin,2)*1. +
                          (1./izd_thin)*powf(1./iRd_thin,2)*norm_thick);
  double MBonMS = 9.23/54.3;
  double Mb = MBonMS*Md/(1.-MBonMS);
  rho0b = 1e-9*97.4*(0.51/0.5)/8.96*Mb;
}


/*---- Thin disc prior -----*/
double Basic_prior::prior_thin() {
  //  if(age>10.) return 0.; // age out of range
  const double isrttpi = 1./ 2.506628274631000502415765284811;
  double exponent = - R*iRd_thin-fabs(z)*izd_thin ;
  double out = isrttpi*exp(exponent);
  return out;
}


/*---- Thick disc prior -----*/
double Basic_prior::prior_thick() {
  //  if(age<8. || age>12.) return 0.; // age out of range
  const double isrttpi = 1./ 2.506628274631000502415765284811;
  double exponent = - R*iRd_thick-fabs(z)*izd_thick;       // position
  double out = norm_thick*isrttpi*exp(exponent);
  return out;
}


/*---- Halo prior -----*/
double Basic_prior::prior_halo() {
  if(r<rcut_halo) return norm_halo*pow(rcut_halo,power_halo);
  // power law diverges at low r - Carollo model applied where it shouldn't
  double out = norm_halo*pow(r,power_halo);
  return out;
}
/*----- Bulge prior ------*/
double Basic_prior::prior_bulge() {
  double rho = rho0b, m = sqrt(R*R+(z/qb*z/qb)), m0 = m/r0b;
  if(gammab != 0.) rho *= powf(m0,-gammab);
  m0+=1.;
  rho *= powf(m0,-betab);
  //cerr << R << ' ' << rho << '\n';
  if(rcutb) rho *= exp(-powf(m/rcutb,2));
  //cerr << R << ' ' << rho << '\n';
  return rho;
}






/*-----------------------------------------------------------------------------
*
* Class: GRVS_prior
*
* See Prior.h
*
*-----------------------------------------------------------------------------*/

/*- setup_flat_DR2() sets prior assuming flat metallicity & age distribution -*/
/*                  Values found in external code */
void GRVS_prior::setup_flat_DR2() {
  // Max/min magnitude (N.B. at larger magnitudes assumed constant prob)
  Mmin = -12;
  Mmax = 7;
  // divisions in piecewise linear fit
  div0_photprior  =  -2.34936054702 ;
  div1_photprior  =  1.04283066069 ;
  div2_photprior  =  3.6743145648 ;
  // Gradients in each piece
  m0_photprior  =  0.636749185704 ;
  m1_photprior  =  0.1865447312 ;
  m2_photprior  =  0.599726656521 ;
  m3_photprior  =  0.144988560077 ;
  // Red clump parameters
  RC_cen  =  -0.435918445572 ;
  RC_sig  =  0.23199035803 ;
  RC_norm  =  8185.80862358 ;
  // Normalisation (so that red clump is the right size)
  c1_photprior   =  3.40503108665 ;

  // ensure continuity
  c2_photprior = c1_photprior+(m1_photprior- m2_photprior)*div1_photprior;
  c3_photprior = c2_photprior+(m2_photprior - m3_photprior)*div2_photprior;
  c0_photprior = c1_photprior+(m1_photprior- m0_photprior)*div0_photprior;
}



/*- setup_flat_DR2() sets prior thin-disc-like distribution -*/
/*                  Values found in external code            */
void GRVS_prior::setup_thin_DR2() {

    // Max/min magnitude (N.B. at larger magnitudes assumed constant prob)
  Mmin = -12;
  Mmax = 7;



  // divisions in piecewise linear fit
div0_photprior  =  -2.02022468025 ;
div1_photprior  =  1.27873640026 ;
div2_photprior  =  3.5343533639 ;
// Gradients in each piece
m0_photprior  =  0.678238934106 ;
m1_photprior  =  0.213637402684 ;
m2_photprior  =  0.670052484961 ;
m3_photprior  =  0.141224575915 ;
// Red clump parameters
RC_cen  =  -0.349737628723 ;
RC_sig  =  0.150875363233 ;
RC_norm  =  9919.71867715 ;
// Normalisation (so that red clump is the right size)
c1_photprior   =  3.3209869255 ;
// ensure continuity
  c2_photprior = c1_photprior+(m1_photprior- m2_photprior)*div1_photprior;
  c3_photprior = c2_photprior+(m2_photprior - m3_photprior)*div2_photprior;
  c0_photprior = c1_photprior+(m1_photprior- m0_photprior)*div0_photprior;

}




GRVS_prior::GRVS_prior() {
// Anything I need to do to set up the GRVS prior
  setup_thin_DR2();
}

GRVS_prior::GRVS_prior(string type) {
// Anything I need to do to set up the GRVS prior
  if(type=="thin_DR2")      setup_thin_DR2();
  else if(type=="flat_DR2") setup_flat_DR2();
  else {
    cerr << "No prior of that type known - exiting\n";
    exit(1);
  }
}

/* ----- Setup sightline ------- */
void GRVS_prior::prepare_sightline(double l, double b) {
    Non_mag_prior.prepare_sightline(l,b);
}

/* ------ Prior -------- */
double GRVS_prior::prior(double s_in, double GRVS_in) {
  M_GRVS = GRVS_in - 5*log10(s_in/10.);
  return prior_mag() * Non_mag_prior.prior(s_in);
}

/* ------ Prior on G_RVS -------- */
double GRVS_prior::prior_mag() {
  double log_pm;
  if(M_GRVS<Mmin) return 0.;
  else if(M_GRVS<div0_photprior) {
    log_pm = c0_photprior+m0_photprior*M_GRVS;
  }
  else if(M_GRVS<div1_photprior) {
    log_pm = c1_photprior+m1_photprior*M_GRVS;
  } else if(M_GRVS<div2_photprior) {
    log_pm = c2_photprior+m2_photprior*M_GRVS;
  } else if(M_GRVS<Mmax) {
    log_pm = c3_photprior+m3_photprior*M_GRVS;
  } else if(M_GRVS>=Mmax) {
    log_pm = c3_photprior+m3_photprior*Mmax;
  }
  double pm = pow(10.,log_pm);
  pm += RC_norm * exp(-pow(M_GRVS-RC_cen,2)/(2.*RC_sig*RC_sig));

  return pm;

}



Anticentre_prior::Anticentre_prior() {
  setup_basic();
}

void Anticentre_prior::setup_basic(){
  // initial guess based on fit to p/ep>5 stars
  //logscale = -1.1635621154541131; // data initial

  logscale = -1.9746760143252673; // gogum initial
  logscale = -1.325897540666271; // gogum iteration 1
  logscale = -1.2087596366581825;
}

double Anticentre_prior::prior(double s_in)      // Return prior given s
{
  return Non_mag_prior.prior(s_in) * prior_select(s_in);
}


double Anticentre_prior::prior_select(double s_in)      // Return prior given s
{
  return exp(logscale*0.001*s_in);
}


/* ----- Setup sightline ------- */
void Anticentre_prior::prepare_sightline(double l, double b) {
    Non_mag_prior.prepare_sightline(l,b);
}



