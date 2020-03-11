#include <stdio.h>  /* printf, fprintf */
#include <iostream> /* system("pause") */
#include <math.h>   /* pow */
#include <time.h>
#include "const_param.h"
//
int isotops(int const, double const, double *, double *, double *);
void count_lines(int const, double const, double const, int &, int &);
void read_hitran160(int const, int const, int const, int *, double *, double *, double *,
	                double *, double *, double *, double *);
int ix1ix2(double const, double const, double const *, int const, int &, int &);
/*
	humlicek(x, y):
		- location: C:\CODES\Alexei\Common_C_Release
		- requires: cmplx.cpp, Cmplx.h, cmplx.h - different from Cmplx.h? location of cmplx.h is unknown! ??????????????
		  other sources located in C:\CODES\Alexei\Common_C_Release.
		So far, i dont know the details of humlicek(x, y)
*/
double humlicek(double, double);
//
int main()
{
	FILE *pFile;
	int code, molec_id, iline, iline0, nlines, niso, iso_ix, nnu, inu, inu1, inu2; // inu_lo, inu_hi;
	double time_start, time_end,
		   nu_hit_min, nu_hit_max, T_kelvin, p_atm, l_cm, nu_usr_min, nu_usr_max, dnu,
		   e11, e21, e12, e22, SijT, nuij_pshift, alf_doppler, gam_lorentz, x, y, n_column;
	int *isotop_id;
	double *nuij, *Sij, *gamma_air, *gamma_self, *Epp,
		   *n_air, *delta_air, *nu, *k_abs, *tau_abs;
	double Qratio[niso_max], mmass_iso[niso_max], Ia_iso[niso_max];
	char fname_out[fname_len_max] = "tau_abs_gcell.txt";
/*------------------------------------------------------------------------------------------------*/
//
//  Input:
//	gas cell paramters
	molec_id = 6;
	T_kelvin = 296.0; // temperature, [K]
	p_atm = 0.986923;      // pressure, [atm]
	l_cm = 8.0;      // cell length, [cm]
//  spectral range
	nu_usr_min = 4081.901000;
	nu_usr_max = 4505.699001;
	dnu = 0.002;
/*------------------------------------------------------------------------------------------------*/
//
	time_start = (double)clock() /(double)CLOCKS_PER_SEC;
	nu_hit_min = nu_usr_min - delta_nu;
	nu_hit_max = nu_usr_max + delta_nu;
//  Column (areal) *number* density: *number* of particles per unit area, [cm-2]
//  [k.T] = [erg]; 1erg = 1 g.cm2/s2: https://en.wikipedia.org/wiki/Erg 
	n_column = l_cm*atm_to_cm_g_s*p_atm/(k_boltzman*T_kelvin); // n_column = l_cm * n_volume
//  Echo some parameters
	printf("PROJECT=tau_abs_gas_cell:");
	printf("\n-delta_nu=%6.2f [cm-1]", delta_nu);
//
//  User-defined spectral grid, nu, and spectral absorption coefficient, k_abs
	nnu = ceil((nu_usr_max - nu_usr_min)/dnu);
	nu = new double [nnu];
	k_abs = new double [nnu];
	tau_abs = new double [nnu];
	for (inu = 0; inu < nnu; inu++) {
		nu[inu] = nu_usr_min + inu*dnu;
		k_abs[inu] = 0.0;
		tau_abs[inu] = 0.0;
	}
//
//  Read HITRAN data, compute TIPS_ratio
	count_lines(molec_id, nu_hit_min, nu_hit_max, iline0, nlines);
	printf("\niline0 = %i, nlines = %d", iline0, nlines);
//
	isotop_id = new int[nlines];
	nuij = new double[nlines];
	Sij = new double[nlines];
	gamma_air = new double[nlines];  // not used for gas cell (one gas)
	gamma_self = new double[nlines];
	Epp = new double[nlines];
	n_air = new double[nlines];
	delta_air = new double[nlines];
//
	read_hitran160(molec_id, iline0, nlines, isotop_id, nuij, Sij,
				       gamma_air, gamma_self, Epp, n_air, delta_air);
//
	niso = isotops(molec_id, T_kelvin, Qratio, mmass_iso, Ia_iso);
//
	for (iline = 0; iline < nlines; iline++) {
		if (ix1ix2(nuij[iline], delta_nu, nu, nnu, inu1, inu2) > 0) {
			iso_ix = isotop_id[iline]-1;
			nuij_pshift = nuij[iline] + delta_air[iline]*p_atm; // Alexei checks for non-zero pressure shift delta_air > tiny
			alf_doppler = nuij_pshift*sqrt(2.0*n_avogadro*k_boltzman*T_kelvin*ln2/mmass_iso[iso_ix])/c_light; // save many operations
			gam_lorentz = pow(T_ref/T_kelvin, n_air[iline])*gamma_self[iline]*p_atm; // 1 gas: p_atm == p_self, gamma_air is not used
			y = sqrt_ln2*gam_lorentz/alf_doppler;
//          Temperature-corrected line intensity
//          THINKME: combine e11 & e21 and save one exp call
	        e11 = exp(-c2_rad*Epp[iline]/T_kelvin);
	        e21 = exp(-c2_rad*Epp[iline]/T_ref);
	        e12 = 1.0 - exp(-c2_rad*nuij_pshift/T_kelvin); // nuij_pressure_shift ??????????????????
	        e22 = 1.0 - exp(-c2_rad*nuij_pshift/T_ref);    // nuij_pressure_pshift ?????????????????
//          THINK ME: Alexei if-checks for 1-exp(-big_number) = 1
	        SijT = Ia_iso[iso_ix]*Sij[iline]*Qratio[iso_ix]*e11*e12/e21/e22;
			for (inu = inu1; inu < inu2+1; inu++) { // note: inu2 is included
				x = sqrt_ln2*fabs(nu[inu] - nuij_pshift)/alf_doppler; // save operation sqrt(ln2)/gam_D
				k_abs[inu] += SijT*(sqrt_ln2/sqrt_pi)*humlicek(x, y)/alf_doppler; // in hum(x, y) make sure variables are in the right order !!!!!!!!!!!!!!!!!!!!!
			} // inu = inu1 : inu2
		} // if iline contributes to 'nu'
		if (!(iline%1000)) printf("\n-processed: iline=%d", iline);
	} // for iline
//
	pFile = fopen(fname_out, "w");
	fprintf(pFile, "nlines=%i\n", nlines); // using tau_abs is not necessary
	for (inu = 0; inu < nnu; inu++)
	{
		tau_abs[inu] = k_abs[inu]*n_column;
		fprintf(pFile, "%9.3f  %12.6e\n", nu[inu], tau_abs[inu]); // using tau_abs is not necessary
	}
	fclose(pFile);
//
	delete[] tau_abs;
	delete[] k_abs;
	delete[] nu;
	delete[] isotop_id;
	delete[] nuij;
	delete[] Sij;
	delete[] gamma_air;
	delete[] gamma_self;
	delete[] Epp;
	delete[] n_air;
	delete[] delta_air;
//
	time_end = (double)clock() /(double)CLOCKS_PER_SEC;
	printf("\ncpu total time %6.2fs\n", time_end - time_start);
	system("\npause");
}
/*------------------------------------------------------------------------------
yy/mm/dd - First created
------------------------------------------------------------------------------*/