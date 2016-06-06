#ifndef DMTPC_THEORY_H
#define DMTPC_THEORY_H

#include "TROOT.h"

/** Functions for calculating dark matter theory distributions from Lewin & Smith \author Asher Kaboth

*/

namespace DmtpcTheory{

   /**
      Conversts a cross section to a rate
      \param xsect cross section
      \param v_0 rms of dark matter velocity
      \param A A of target nucleus
      \param M_dark mass of dark matter particle
      \param rho_D density of dark matter
    */
   double ConvertXsectToRate( const Double_t xsect , const Double_t v_0, 
			      const Double_t A, const Double_t M_dark, 
			      const Double_t rho_D );
   /** 2D distribution with no escape velocity; incorrectly normalized
       \par x 2-D array containing energy, cos(theta_cygnus)
       \par parameters: 0-dark matter mass, 1-A of target, 2-mass of target 3-v_0, 4-v_earth, 5-cross section, 6-density of dark matter

   */
   double dRdERbetweenVearthandInfinity2DOverR0(Double_t *x, Double_t *par);
   /** 2D distribution with no escape velocity; correctly normalized
       \par x 2-D array containing energy, cos(theta_cygnus)
       \par parameters: 0-dark matter mass, 1-A of target, 2-mass of target 3-v_0, 4-v_earth, 5-cross section, 6-density of dark matter

   */
   double dRdERbetweenVearthandInfinity2D(Double_t *x, Double_t *par);
   /** 1D distribution with no escape velocity; incorrectly normalized
       \par x 1-D array containing energy
       \par parameters: 0-dark matter mass, 1-A of target, 2-mass of target 3-v_0, 4-v_earth, 5-cross section, 6-density of dark matter

   */
   double dRdERbetweenVearthandInfinity1DOverR0(Double_t *x, Double_t *par);
   /** 1D distribution with no escape velocity; correctly normalized
       \par x 1-D array containing energy
       \par parameters: 0-dark matter mass, 1-A of target, 2-mass of target 3-v_0, 4-v_earth, 5-cross section, 6-density of dark matter

   */
   double dRdERbetweenVearthandInfinity1D(Double_t *x, Double_t *par);
   /** 1D distribution with escape velocity; incorrectly normalized
       \par x 1-D array containing energy
       \par parameters: 0-dark matter mass, 1-A of target, 2-mass of target 3-v_0, 4-v_earth, 5-cross section, 6-density of dark matter, 7-v_escape

   */
   double dRdERbetweenVearthandVesc1DOverR0(Double_t *x, Double_t *par);
   /** 1D distribution with escape velocity; correctly normalized
       \par x 1-D array containing energy
       \par parameters: 0-dark matter mass, 1-A of target, 2-mass of target 3-v_0, 4-v_earth, 5-cross section, 6-density of dark matter, 7-v_escape

   */
   double dRdERbetweenVearthandVesc1D(Double_t *x, Double_t *par);
      /** 1D distribution with no escape velocity; incorrectly normalized
       \par x 1-D array containing cos(theta_cygnus)
       \par parameters: 0-dark matter mass, 1-A of target, 2-mass of target 3-v_0, 4-v_earth, 5-cross section, 6-density of dark matter, 7-energy

   */
   double dRdThOverR0(Double_t *x, Double_t *par);
      /** 1D distribution with no escape velocity; correctly normalized
       \par x 1-D array containing cos(theta_cygnus)
       \par parameters: 0-dark matter mass, 1-A of target, 2-mass of target 3-v_0, 4-v_earth, 5-cross section, 6-density of dark matter, 7-energy

   */
   double dRdTh(Double_t *x, Double_t *par);
      

}

#endif
