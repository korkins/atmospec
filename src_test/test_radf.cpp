#include <stdio.h>  /* printf, fprintf */
#include <iostream> /* system("pause") */
#include <math.h>
#include <time.h>
//
double radf(double const, double const);
double RADFN(double vi, double xkt);
//
int main()
{
	int const
		nnu = 50000, nT = 200;
	double const
		RADCN2 = 1.4387770,
		nu0 = 20000.0,
		T0 = 296.0;
	int
		inu, iT;
	double
		s, a, dif, t1_sec, t2_sec;
	double
		nu[nnu], T[nT];
	s = 0.0;
	a = 0.0;
	dif = 0.0;
	for (inu = 0; inu < nnu; inu++)
	{
		nu[inu] = 1.0 + inu;
		for (iT = 0; iT < nT; iT++)
		{
			T[iT] = 200.0 + iT;
			dif += fabs( radf(nu[inu], T[iT]) - RADFN(nu[inu], T[iT]/RADCN2) );
		}
	}
	printf("\n dif=%6.1e", dif/nnu/nT);
	s = radf(nu0, T0);
	a = RADFN(nu0, T0/RADCN2);
	printf("\n 2 values: %12.6e, %12.6e, %12.6e", s, a, s - a);
//
	t1_sec = (double)clock() /(double)CLOCKS_PER_SEC;
	s = 0.0;
	for (inu = 0; inu < nnu; inu++)
		for (iT = 0; iT < nT; iT++)
			s += radf(nu[inu], T[iT]);
	t2_sec = (double)clock() /(double)CLOCKS_PER_SEC;
	printf("\n runtime for S: %6.2fs\n", t2_sec - t1_sec);
//
	t1_sec = (double)clock() /(double)CLOCKS_PER_SEC;
	a = 0.0;
	for (inu = 0; inu < nnu; inu++)
		for (iT = 0; iT < nT; iT++)
			a += RADFN(nu[inu], T[iT]/RADCN2);
	t2_sec = (double)clock() /(double)CLOCKS_PER_SEC;
	printf("\n runtime for A: %6.2fs\n", t2_sec - t1_sec);
//
	dif = s - a;
	printf("\nDone!\n");
	system("pause");
	return 0;
}