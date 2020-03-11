#include <math.h> /* exp, log */
//
double intexp(double const *f, double const *x, int const nx) {
/*--------------------------------------------------------------------------------------------------
PURPOSE:
	To integrate f(x) analyticaly using approximation: f(x) = f0*exp(g*(x-x0)), x0 = x[0].
INPUT:
	f   d[nx]   f(x)
	x   d[nx]   x-points
	nx  i       total number of points, nx > 1
OUTPUT:
	ie   d       result of integration
COMMENT:
	log(f/f0) = g(x-x0) leads to a 1-parameter linear least squares problem [1]:
		y = gx -> xTy = gxTx -> g = xTy/xTx
	WARNING: g = 0 is NOT processed separately
REFERENCESS:
	1.https://en.wikipedia.org/wiki/Linear_least_squares
PROTOTYPE:
	double intexp(double const *, double const *, int const);
--------------------------------------------------------------------------------------------------*/
	int
		ix;
	double
		f0, x0, xTx, xTy, g, ie;
//--------------------------------------------------------------------------------------------------
	f0 = f[0];
	x0 = x[0];
	xTx = 0.0;
	xTy = 0.0;
	for (ix = 1; ix < nx; ix++) { // skip ix=0: x-x0=0, ln(f[0]/f0)=0
		xTx += (x[ix] - x[0])*(x[ix] - x[0]);
		xTy += (x[ix] - x[0])*log(f[ix]/f0);
	} // for ix = 1:nx-1
	g = xTy/xTx;
//  THINKME: g=0 - f(x) == constant; if g < tiny*(x_max + x_min)/2 then integrate const
	ie = (f0/g)*(exp(g*(x[nx-1] - x0)) - 1.0);
	return ie;
} // intexp
/*--------------------------------------------------------------------------------------------------
2020-01-30:
	First created and tested
	A) analytical integration for f(x) = c0*exp(g(x-x0))
		1) c0=+1, g=-2, x=0:1:10, x0=0
		2) c0=-1, g=+2, x=0:1:10, x0=0
		3) c0=+9.0, g=-0.5, x=0:1:10, x0=0
		4) c0=+9.0, g=-0.5, x=0:0.001:1, x0=0
		5) c0=+3.0, g=-0.125, x=5:0.001:6, x0=5
	B) TO BE TESTED: Interpolation of the MODTRAN profiles
	C) TO BE TESTED: use exp-interpolation in Excel, get coeficients
--------------------------------------------------------------------------------------------------*/