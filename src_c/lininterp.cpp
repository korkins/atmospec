void lininterp(int const nx, double const *x, double const *y,                                //  in
	           int const np, double const *xp,
			   double *yp) {															      // out
/*--------------------------------------------------------------------------------------------------
PURPOSE:
	To compute piecewise linear interpolation.
INPUT:
	nx   i       number of grid points
	x    d[nx]   abscissae: x[ix] < x[ix+1]
	y    d[nx]   y(x)
	np   i       number of interpolation points
	xp   d[np]   interpolation points: xp[ip] < xp[ip+1]
OUTPUT:
	yp   d[np]   y(xp)
COMMENT:
	'yp' is defined outside.
	No extrapolation, unlike e.g. in [1]
REFERENCESS:
	1.https://docs.scipy.org/doc/numpy/reference/generated/numpy.interp.html
PROTOTYPE:
	void lininterp(int const, double const *, double const *, int const, double const *, double *);
--------------------------------------------------------------------------------------------------*/
	int
		ip, ix;
	double
		x_min, xp_min, xp_max, x0, x1, y0, y1, k;
//--------------------------------------------------------------------------------------------------
//
	x_min = x[0];
	xp_min = xp[0];
	xp_max = xp[np-1];
	x0 = x[0];
	y0 = y[0];
	x1 = x[1];
	y1 = y[1];
	ix = 1;
	ip = 0;
	if (x_min <= xp_min)
	{
		while (x1 < xp_min && ix < nx) // move x-bins right
		{
			x0 = x1;
			y0 = y1;
			ix += 1;
			x1 = x[ix];
			y1 = y[ix];
		} // while (x1 < xp_min && ix < nx)
	}// if x_min < xp_min
	else
		while(xp[ip] < x_min) // move xp-points right
			ip += 1;
//
	while (x0 < xp_max && ix < nx)
	{
		k = (y1 - y0)/(x1 - x0);
		while (xp[ip] <= x1 && ip < np)
		{
			yp[ip] = y0 + k*(xp[ip] - x0);
			ip += 1;
		} // while (xp[ip] <= x1 && ip < np)
		x0 = x1;
		y0 = y1;
		ix += 1;
		x1 = x[ix];
		y1 = y[ix];
	} // while (x1 <= xp_max && ix < nx)
} // lininterp
/*--------------------------------------------------------------------------------------------------
2020-01-17:
	First created & tested: x = [-10:1:10] for ALL tests
	A) xp within x :: x_min < xp_min < xp_max < x_max:
	   y = x
	   xp = [-9.5 -5.0 -4.9 -4.5 -4.1 0.0 1.0 1.5 2.0 2.2 2.8 3.0 5.5 9.9]
	   yp = xp - ok.
	B) x within xp :: xp_min < x_min < x_max < xp_max:
	   x = [-10:1:10]
	   y = x
	   xp = [-12 -11 -9.5 -5.0 -4.9 -4.5 -4.1 0.0 1.0 1.5 2.0 2.2 2.8 3.0 5.5 9.9 11 12]
	   yp = [0 0 xp 0 0] - ok
	C) General test:
	   y = sin(x)
	   xp = xp(case B)
	   xp, yp =  -12.0   0.00000000e+00 <- defined outside
                 -11.0   0.00000000e+00 <- defined outside
                  -9.5   6.59513128e-02
                  -5.0   9.58924275e-01
                  -4.9   9.38712097e-01
                  -4.5   8.57863385e-01
                  -4.1   7.77014673e-01
                   0.0   0.00000000e+00
                   1.0   8.41470985e-01
                   1.5   8.75384206e-01
                   2.0   9.09297427e-01
                   2.2   7.55661943e-01
                   2.8   2.94755492e-01
                   3.0   1.41120008e-01
                   5.5  -6.19169886e-01
                   9.9  -4.48407151e-01
                  11.0   0.00000000e+00 <- defined outside
                  12.0   0.00000000e+00 <- defined outside
		Tested vs numpy.interp - ok: abs dif ~e-10 (for 9 digits), except for
		2 first and 2 last points where python extrapolates values.
*/