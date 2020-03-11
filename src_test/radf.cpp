#include <math.h> /* exp */
//
double radf(double const nu, double const Tkelv) {
/*--------------------------------------------------------------------------------------------------
PURPOSE:
	To compute the radiation field term for continuum
IN:
	nu      d   wavenumber in cm-1
	Tkelv   d   temperatute in K
OUT:
	radf    d   radiation field, [radf]=[nu]=cm-1
NOTE:
    See [1: Eqs.3-5] for equations and units.
REFS:
	1. Clough SA, Kneizys FX, and Davies RW, 1989: Atmos.Res. 23, p.229
PROTOTYPE:
	double simpson(double const *, int const, double const);
--------------------------------------------------------------------------------------------------*/
	double const
		c2 = 1.4387770; // 2nd radiation constant, c2=hc/k (cm.K), in
	double
		e;
//--------------------------------------------------------------------------------------------------
	e = exp(-c2*nu/Tkelv);
	return nu*(1.0 - e)/(1.0 + e);
} // radf
/*--------------------------------------------------------------------------------------------------
2020-02-13:
	First created and tested vs original subroutine (below) nu = 1:1:50000, T=200:1:400, including
	nu0=20000, To=296K printed out explictely. untime in debug S=0.68s, A=0.40; in release S=A=0.00!
----------------------------------------------------------------------------------------------------
Original subroutine from Alexei's continuum.cpp

// FUNCTION RADFN CALCULATES THE RADIATION TERM FOR THE LINE SHAPE //

double RADFN(double vi, double xkt)
{
// Korkin: see e.g. Clough Kneizys Davies 1989: Atm.Res. 23, pp229-241, Eq.(2)
//	abs_coef = k(nu) ~ R(nu)=nu*(1 - exp(x))/(1 + exp(x)),
//	where x = -(hc/kT)*nu = -RADCN2*nu/T

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
*/