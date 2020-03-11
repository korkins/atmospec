/****************************************************************************************************/
/*   
	Continuum model:  Water Vapor & CO2:
    { revision:  $Revision: 1.00 $
      created:   $Date: 2003/02/12 18:50:04 $

     This  version of the water vapor continuum, mt_ckd_1.00, is a completely new 
     continuum model based on a collision induced  component and a sub-Lorentzian 
     line wing component.  Continua related to other species are the same as for 
     ckd_2.4.2.
	 }

	This program calculates absorption OT of H2O, CO2 and NO2 continua for spectral range 
	  NuBeg - NuEnd with step d_NU. The continuum absorption is calculated with step DV_H2O=10 cm-1.
	  At arbitrary frequency, it is found by linear interpolation.

  The tau_absorption profile is computed in ib (SHARM) points (tau_abs[0]=0).
	Integration over altitude is performed using Simpson's quadrature.
	Spectral range of H2O continuum: 0 - 20000 cm-1, step 10 cm-1.
	Spectral range of CO2 continuum: 0 - 10000 cm-1, step 10 cm-1.
	Spectral range of NO2 continuum (12593.73 - 25000 with variable step 3-7 cm-1).
	
	Ozone absorption is calculated separately after calculation of average vertical profiles.
*/

#include <math.h>
#include <stdio.h>

/* FUNCTION RADFN CALCULATES THE RADIATION TERM FOR THE LINE SHAPE */

double RADFN(double vi, double xkt)
{
/* Korkin: see e.g. Clough Kneizys Davies 1989: Atm.Res. 23, pp229-241, Eq.(2)
	abs_coef = k(nu) ~ R(nu)=nu*(1 - exp(x))/(1 + exp(x)),
	where x = -(hc/kT)*nu = -RADCN2*nu/T
*/

  double f, var, expon;
	
	var = vi/xkt; // Korkin: xkt = t0[iatm][k]/RADCN2; vi = nu => var = nu*RADCN2/T
	if(var <= 0.01)
		f = 0.5 * vi * var; // Korkin: Taylor
	else if(var > 0.01 && var <= 10.)
	{	expon = exp(-var);
		f = vi*(1.- expon)/(1.+ expon); // Korin: General expression, Eq.(2)
	}
	else
		f = vi; // exp(-BIG) = 0 => R(nu) = nu
	return f;
}
