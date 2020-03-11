#include <stdio.h>  /* printf, fprintf */
#include <iostream> /* system("pause") */
//
double simpson(double const *, int const, double const);
double simpson_alexei(double *, int , double);
//
int main()
{
	int const nx = 101;
	int ix;
	double x, dx, f[nx], sim_S, sim_A;
//
	int c24, sign;

	c24 = 4;
	sign = 1;
	for(ix = 2; ix < 9; ix++) {
		sign = -sign;
		c24 = c24 + sign*2;
		printf("\n c24=%i", c24);
	}
//
	dx = 0.1;
	for (ix = 0; ix < nx; ix++) {
		x = -5.0 + ix*dx;
		//f[ix] = x*x*x*x;
		f[ix] = fabs(sin(x));
	}
	sim_S = simpson(f, nx, dx);
	sim_A = simpson_alexei(f, nx, dx);
	printf("\n z=%12.6e,    dif=%6.1e", sim_S, sim_S - sim_A);

	printf("\nDone!\n");
	system("pause");
	return 0;
}