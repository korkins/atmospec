#include <cstring>       /* strcpy, strcat */
#include <stdio.h>       /* FILE */
#include "const_param.h" /* char path_tips[], char fname_hitdb[], int path_len_max, int fname_len_max */
//
int isotops(int const molec_id, double const T_kelvin,
	        double *Qratio, double *mmass_iso, double *Ia_iso) {
/*--------------------------------------------------------------------------------------------------
PURPOSE:
	To compute the total internal partition sum (TIPS) ratio, Q(296K)/Q(T), for a given temperature
	T[K], and return for each isotop the number of isotops, molar mass, and abundance.
INPUT:
	molec_id   i   as in HITRAN
	T_kelvin   d   temperature, [K]
OUTPUT:
	niso        i             number of isotops
	Qratio      d[niso_max]   Q(296K)/Q(T): *** WARNING: reference/arbitrary - not vice versa! ***
	mmass_iso   d[niso_max]   molar mass of isotops [g/mol]
	Ia_iso      d[niso_max]   natural terrestrial isotopic abundances
COMMENT:
    For output arrays,  only first 'niso' values are used.
	Q(T) is linear intepolated.
	Q(296K) comes from [1], not from [2] (less digits).
	See [3,4] for the definition of TIPS, and [5,6] for a Fortran code, TIPS.for.
REFERENCESS:
	1. https://hitran.org/docs/iso-meta/
	2. https://hitran.org/media/molparam.txt
	3. https://hitran.org/docs/definitions-and-units/
	4. Gamache RR, et al., 2017: JQSRT 203, 70-87 (doi: 10.1016/j.jqsrt.2017.03.045)
	5. https://hitran.org/supplementary/
	6. https://hitran.org/suppl/TIPS/
PROTOTYPE:
	int isotops(int const, double const, double *, double *, *double);
--------------------------------------------------------------------------------------------------*/
	FILE *fin;
	char
		fpath[path_len_max], fname_iso[niso_max][fname_len_max];
	int
		iso, niso;
	double
		Q, Q1, Q2, T1, T2;
	double
		Qref[niso_max];
//--------------------------------------------------------------------------------------------------
//
	niso = -999;
	T1 = 0.0;
	Q1 = 0.0;
//
	switch (molec_id)
	{
		case 6: // ch4
			niso = 4;
			for (iso = 0; iso < niso; iso++)
			{
				Qref[iso] = Qref_ch4[iso];
				mmass_iso[iso] = molar_mass_ch4[iso];
				Ia_iso[iso] = Ia_iso_ch4[iso];
				strcpy(fname_iso[iso], fname_iso_ch4[iso]);
			} // for iso
			break;
		default:
			printf("\n(in isotops.cpp) default is undefined");
	} // switch (molec_id)
//
	for (iso = 0; iso < niso; iso++)
	{
		strcpy(fpath, path_TIPS);
		strcat(fpath, fname_iso[iso]);
		fin = fopen(fpath, "r");
		fscanf(fin, "%lf %lf", &T2, &Q2);
		while (T2 < T_kelvin)
		{
//          THINKME: potentially infinite loop if T_kelvin > T_grid_max
//          THINKME: skip T from 1 to 190 as in Alexei's code? Cut lut from top
			T1 = T2;
			Q1 = Q2;
			fscanf(fin, "%lf %lf", &T2, &Q2);
		} // while T2 < T_kelvin
		Q = Q1 + (T_kelvin - T1)*(Q2 - Q1)/(T2 - T1);
		Qratio[iso] = Qref[iso]/Q;
		fclose(fin);
		printf("\n(in isotops.cpp) isotop file processed: %s", fpath);
	} // for iso
//
	return(niso);
} // isotops
/*--------------------------------------------------------------------------------------------------
2020-01-09:
	Renamed: molar_mass_iso -> mmass_iso
	Added:  Ia_iso
	Minor:  printf removed/modified
	Tested: compiled - ok.
2020-01-07:
    First created and tested for molec_id=6(ch4), T=Tref=296K -> Qratio[:]=1.0
	T=Tref+0.5K -> Q=0.5*(Q1+Q2) - ok
	niso(out) = 4 - ok
	mmass_iso = {16.0313, 17.034655, 17.037475, 18.04083} - ok
*/