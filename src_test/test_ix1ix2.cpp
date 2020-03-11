#include <stdio.h>  /* printf, fprintf */
#include <iostream> /* system("pause") */
//
int ix1ix2(double const, double const, double const *, int const, int &, int&);
//
int main()
{
	int const nx = 101;
	int ix, ix1, ix2, code;
	double x0, dx, x[nx];
//
	int c24, sign;
//
	x0 = 13.0;
	dx = 0.0;
//
	for (ix = 0; ix < nx; ix++)
		x[ix] = 1.0*ix;
//
	code = ix1ix2(x0, dx, x, nx, ix1, ix2);
	if (code > 0)
		printf("\n ix1=%i, ix2=%i", ix1, ix2);
	else
		printf("\n code=%i", code);

	c24 = 4;
	sign = 1;
	for(ix = 2; ix < 9; ix++) {
		sign = -sign;
		c24 = c24 + sign*2;
		//if (ix % 2 == 0)
		//	c24 = 2;
		//else
		//	c24 = 4;
		printf("\n c24=%i", c24);
	}


	printf("\nDone!\n");
	system("pause");
	return 0;
}