#include <stdio.h>  /* printf, fprintf */
#include <iostream> /* system("pause") */
#include <math.h>
//
double intexp(double const *, double const *, int const);
//
int main()
{
	int const nx = 1001;
	int ix;
	double c0, x0, g, dx, x[nx], f[nx], val1, val2;
//
	dx = 0.001;
	for(ix = 0; ix < nx; ix++)
		x[ix] = 5.0 + ix*dx;
	x0 = x[0];
	g = -0.125;
	c0 = +3.0;
	for(ix = 0; ix < nx; ix++)
		f[ix] = c0*exp(g*(x[ix]-x0));
	val1 = (c0/g)*( exp(g*(x[nx-1]-x0)) - 1.0 );
	val2 = intexp(f, x, nx);
	printf("\n exact=%12.4e, ratio=int/exct=%10.6f", val1, val2/val1);
	printf("\nDone!\n");
	system("pause");
	return 0;
}