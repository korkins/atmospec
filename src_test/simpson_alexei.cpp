/* Program calculates the integral using n-point Simpson quadrature.
	n - number of equidistant points. Must BE ODD (n>=3). 
	rab[0 : n-1] - integrand function
	dx - distance between points (delta x)
	*/
double simpson_alexei(double *rab, int n, double dx)
{
  int j;
  double x, z;

	x = rab[0] + 4*rab[1] + rab[n-1];

	for(j=2; j<n-1; j++)
	{
		if(j%2 == 0)	z = 2;	// even number
		else			z = 4;	// odd number
		x += z*rab[j];
	}
	return x*dx/3;
}