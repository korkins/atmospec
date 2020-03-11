#include <stdio.h>  /* printf, fprintf */
#include <iostream> /* system("pause") */
#include <math.h>   /* sin */
//
void lininterp(int const, double const *, double const *, int const, double const *, double *);
//
int main()
{
	int const nx = 21, np = 14+4;
	int ix, ip;
	double
		x[nx] = {-10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0,
		          0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0},
		//y[nx] = {-10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0},
		//xp[np] = {-9.5, -5.0, -4.9, -4.5, -4.1, 0.0, 1.0, 1.5, 2.0, 2.2, 2.8, 3.0, 5.5, 9.9},
		xp[np] = {-12.0, -11.0, -9.5, -5.0, -4.9, -4.5, -4.1, 0.0, 1.0, 1.5, 2.0, 2.2, 2.8, 3.0, 5.5, 9.9, 11.0, 12.0},
		*y,
		*yp;

	yp = new double[np];
	for (ip = 0; ip < np; ip++)
		yp[ip] = 0.0;
	y = new double[nx];
	for (ix = 0; ix < nx; ix++)
		y[ix] = sin(x[ix]);
	lininterp(nx, x, y, np, xp, yp);
	printf("\n xp, yp:");
	for (ip = 0; ip < np; ip++)
		printf("\n %3.1f %16.8e", xp[ip], yp[ip]);
	printf("\nDone!\n");
	system("pause");
	return 0;
}