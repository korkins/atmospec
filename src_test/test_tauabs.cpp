#include <stdio.h>  /* printf, fprintf */
#include <iostream> /* system("pause") */
#include <math.h>
//
//void tauabs25(double const *ext_km, double const *ztau, int const ntau, double *tau);
void testarr(double *a, int const na);
//
int main()
{
	int const
		ntau = 10, nz_mod = 50;
	int
		iz;
	double
		ztau[ntau] = {25.0, 24.9, 24.1, 23.5, 22.5, 21.5, 13.3, 2.9, 0.5, 0.0};
	    //{0.0, 0.5, 2.9, 13.3, 21.5, 22.5, 23.5, 24.1, 24.9, 25.0};
	double
		*ext_km, *tau, *atau;
//
	testarr(atau, ntau);
//
	ext_km = new double [nz_mod];
	for (iz = 0; iz < nz_mod; iz++)
		ext_km[iz] = 1.0;
//	tau = new double [ntau];
//	tauabs25(ext_km, ztau, ntau, tau);


    //printf("\n iparab = %10.6f", ipar);
	printf("\nDone!\n");
	system("pause");
	return 0;
}
//
void testarr(double *a, int const na) {
	int ia;
	a = new double [na];
	for (ia = 0; ia < na; ia++) a[ia] = -999.0;
}