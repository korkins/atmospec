#include <stdio.h>  /* printf, fprintf */
#include <iostream> /* system("pause") */
#include <cstring>  /* strcpy, strcat */
#include <time.h>
#include "const_param.h"
#include "hprofiles.h"
//
int hisotops(int const, double const *, int const, double *, double *, double *);
double simpson(double const *, int const, double const);
double intparab(double const, double const, double const, double const, double const,
	                                        double const, double const, double const);
void count_lines(int const, double const, double const, int &, int &);
void read_hitran160(int const, int const, int const, int *, double *, double *, double *,
	                double *, double *, double *, double *);
void kabs(int const,
	      double const *, int const,
	      int const, double const *, double const *, double const *,
	      int const, int *, double *, double *, double *,
	                        double *, double *, double *, double *,
	      double *);
void cabs(int const, double const *, int const,
	      double const *, double const *, double const *,
	      int const, double *);
void tauabs25(double const *, double const *, double const *, int const, double *);
//
/*
http://modtran.spectral.com/modtran_faq
*/
int main()
{
	FILE *pFile;
	int
		molec_id, iatm, iz, iline0, nlines, nnu, inu, ik, ntau;
	double
		time_start, time_end, hint0_24km, hint24_25km, hint25_50km, hint50_120km,
		nmolec, vol_m3, vol_cm3, scol_cm2, hcol_cm,
		nu_usr_min, nu_usr_max, dnu, nmolec_all, column_amount_usr, column_amount_mod, scalef,
		nu_hit_min, nu_hit_max, WV_mass_density, hcol_water_cm, c_stp,
		t1_sec, t2_sec;
	int
		*isotop_id;
	double
		*Patm, *Tkelv, *conc_cm3, *ppmv, *Pgas;
	double
		*nuij, *Sij, *gamma_air, *gamma_self, *Epp, *n_air, *delta_air, *nu,
		Qratio[niso_max*nz_mod], mmass_iso[niso_max], Ia_iso[niso_max], *knu, *gas_ratio, *ext_km,
		*tau_abs, *ztau, *conc_all, *kcont, *atmcm_km;
	char
		txt_ext[5] = ".txt",
		fname_out[fname_len_max] = "tau_abs_hprofile.txt";
/*------------------------------------------------------------------------------------------------*/
//
//  Input:
	molec_id = 7;
	iatm = 6;
	column_amount_usr = -1.9; // for standard amount use column_amount < 0
	nu_usr_min = 13000.0; 4200.0; floor(10000.0/2.5);
	nu_usr_max = 13050.0; 4450.1; ceil(10000.0/1.5);
	dnu = 1.0/10000; // dnu=0.1: t=200s.; dnu=0.01: t=2,000s
	ntau = 5;
	tau_abs = new double [ntau];
	ztau = new double [ntau];
	ztau[0] = 0.0; ztau[1] = 2.0; ztau[2] = 5.0; ztau[3] = 10.0; ztau[4] = 20.0; 
/*------------------------------------------------------------------------------------------------*/
//
	time_start = (double)clock() /(double)CLOCKS_PER_SEC;
	printf(" PROJECT: tau_abs_hprofile");
	printf("\n atmosphere: %i - %s", iatm, atmos_name[iatm-1]);
	printf("\n molecule: %i - %s", molec_id, molec_name[molec_id-1]);
	printf("\n number of altitudes: %i", ntau);
	printf("\n altitudes, z(km):");
	for (iz = 0; iz < ntau; iz++) printf("%8.3f", ztau[iz]);
	printf("\n spectral range: %9.2f - %9.2f (cm-1)", nu_usr_min, nu_usr_max);
	printf("\n output spectral resolution, dnu: %6.4f (cm-1)", dnu);
	printf("\n line wings range, delta_nu (in const_param.h): %5.2f (cm-1)", delta_nu);
	strcpy(fname_out, molec_name[molec_id-1]);
	strcat(fname_out, ".txt");
//
//  User-defined spectral grid, nu, and spectral absorption coefficient, k_abs
	nnu = ceil((nu_usr_max - nu_usr_min)/dnu);
	nu = new double [nnu];
	for (inu = 0; inu < nnu; inu++)
		nu[inu] = nu_usr_min + inu*dnu;
	knu = new double [nnu*nz_mod];
	for (ik = 0; ik < nnu*nz_mod; ik++)
		knu[ik] = 0.0;
//
	gas_ratio = new double [nz_mod];
	switch(molec_id) {
		case 1:
			for (iz = 0; iz < nz_mod; iz++)
				gas_ratio[iz] = H2O_ppmv[iatm-1][iz]/1.0e6;
			break;
		case 2:
			for (iz = 0; iz < nz_mod; iz++)
				gas_ratio[iz] = CO2_ppmv[iz]/1.0e6;
			break;
		case 3:
			for (iz = 0; iz < nz_mod; iz++)
				gas_ratio[iz] = O3_ppmv[iatm-1][iz]/1.0e6;
			break;
		case 4:
			for (iz = 0; iz < nz_mod; iz++)
				gas_ratio[iz] = N2O_ppmv[iatm-1][iz]/1.0e6;
			break;
		case 5:
			for (iz = 0; iz < nz_mod; iz++)
				gas_ratio[iz] = CO_ppmv[iatm-1][iz]/1.0e6;
			break;
		case 6:
			for (iz = 0; iz < nz_mod; iz++)
				gas_ratio[iz] = CH4_ppmv[iatm-1][iz]/1.0e6;
			break;
		case 7:
			for (iz = 0; iz < nz_mod; iz++)
				gas_ratio[iz] = O2_ppmv[iz]/1.0e6;
			break;
		case 10:
			for (iz = 0; iz < nz_mod; iz++)
				gas_ratio[iz] = NO2_ppmv[iz]/1.0e6;
			break;
	} // switch(molec_id)
//
	Patm = new double[nz_mod];
	Tkelv = new double[nz_mod];
	conc_cm3 = new double[nz_mod];
	conc_all = new double[nz_mod];
	Pgas = new double [nz_mod];
	for (iz = 0; iz < nz_mod; iz++)
	{
		Tkelv[iz] = Tkelv_mod[iatm-1][iz];
		Patm[iz] = Pmbar_mod[iatm-1][iz]*mbar_to_atm;
		Pgas[iz] = Pmbar_mod[iatm-1][iz]*mbar_to_atm*gas_ratio[iz];
		conc_cm3[iz] = Dcm3_mod[iatm-1][iz]*gas_ratio[iz];
		conc_all[iz] = Dcm3_mod[iatm-1][iz];
	}
//
	hint0_24km = simpson(conc_cm3, 25, 1.0);
//	hint24_25km = 0.5*(conc_cm3[24] + conc_cm3[25]);
	hint24_25km = intparab(zkm_mod[24], zkm_mod[25], zkm_mod[23], zkm_mod[24], zkm_mod[25],
			                                         conc_cm3[23], conc_cm3[24], conc_cm3[25]);
	hint25_50km = simpson(&conc_cm3[25], 11, 2.5);
	hint50_120km = simpson(&conc_cm3[35], 15, 5.0);
	nmolec = (hint0_24km + hint24_25km + hint25_50km + hint50_120km)*1.0e5; // 1.0e5 acounts for km->cm;
	vol_m3 = nmolec*k_boltzman*To_stp/Po_stp * 1.0e-7; //1.0e-7 for k_B in SI
	vol_cm3 = vol_m3*1.0e6;
	scol_cm2 = 1.0;
	hcol_cm = vol_cm3/scol_cm2;
	printf("\n column height: %10.4e (atm-cm)", hcol_cm);
//
//======================================================================================================
//  Exercise (not used in the code):
//      Convert density profile (MODTRAN) to atm-cm/km profile (Alexei):
//          atm-cm/km = atm-cm/(1e-5.cm) = 1e5.atm-cm/cm
//      where atm-cm (numerator) is the gas column height at the zero level with
//      z0, To_stp, Po_stp and cm(denimonator) is the one at z, Tz, Pz.
//      Using the known equation for ideal gas (n is concentration in molec/m3)
//              p = nkT
//      one gets ratio of 2 gas column heights, h(z0)/h(z):
//              atm-cm/cm = [To_stp/Po_stp] * [Pz/Tz]
//      where
//              Pz/Tz = nk
//      Note, in SI [n] = molec/m3 = molec/(1e-6.cm3) = 1e6 cm-3
//      Combining together, one gets
//              atm-cm/km = 1e5.k.To_stp/Po_stp 1e6.conc_cm3
//      where the Boltzmann constant, k, must be in SI and both sides depend on z.
	c_stp = 1.0e5 * (k_boltzman * 1.0e-7) * To_stp / Po_stp * 1.0e6;
	atmcm_km = new double [nz_mod];
    for (iz = 0; iz < nz_mod; iz++) atmcm_km[iz] = c_stp*conc_cm3[iz]; // compare vs Alexei's profiles
//======================================================================================================
//
	if (molec_id == 1) // H2O in wv_cm
	{
		WV_mass_density = water_molar_mass*n_loschmidt/n_avogadro;
		printf("\n WV mass density: %10.4e (g/cm3)", WV_mass_density);
		hcol_water_cm = hcol_cm*WV_mass_density/water_mass_density;
		printf("\n MODTRAN precipitated water column amount: %6.2f (cm)", hcol_water_cm);
		if (column_amount_usr > 0.0)
		{
			printf("\n user precipitated water column amount: %6.2f (cm)", column_amount_usr);
			scalef = column_amount_usr/hcol_water_cm;
			for (iz = 0; iz < nz_mod; iz++) conc_cm3[iz] *= scalef;
//
			hint0_24km = simpson(conc_cm3, 25, 1.0);
//			hint24_25km = 0.5*(conc_cm3[24] + conc_cm3[25]);
			hint24_25km = intparab(zkm_mod[24], zkm_mod[25], zkm_mod[23], zkm_mod[24], zkm_mod[25],
			                                         conc_cm3[23], conc_cm3[24], conc_cm3[25]);
			hint25_50km = simpson(&conc_cm3[25], 11, 2.5);
			hint50_120km = simpson(&conc_cm3[35], 15, 5.0);
			nmolec = (hint0_24km + hint24_25km + hint25_50km + hint50_120km)*1.0e5; // 1.0e5 acounts for km->cm;
			vol_m3 = nmolec*k_boltzman*To_stp/Po_stp * 1.0e-7; //1.0e-7 for k_B in SI
			vol_cm3 = vol_m3*1.0e6;
			scol_cm2 = 1.0;
			hcol_cm = vol_cm3/scol_cm2;
			hcol_water_cm = hcol_cm*WV_mass_density/water_mass_density;
			printf("\n check column amount: %6.2f (cm)", hcol_water_cm);
		} // if (column_amount_usr > 0.0)
	} // if (molec_id == 1)
	else // species other than H2O
	{
		hint0_24km = simpson(conc_all, 25, 1.0);
//		hint24_25km = 0.5*(conc_all[24] + conc_all[25]);
		hint24_25km = intparab(zkm_mod[24], zkm_mod[25], zkm_mod[23], zkm_mod[24], zkm_mod[25],
			                                         conc_all[23], conc_all[24], conc_all[25]);
		hint25_50km = simpson(&conc_all[25], 11, 2.5);
		hint50_120km = simpson(&conc_all[35], 15, 5.0);
		nmolec_all = (hint0_24km + hint24_25km + hint25_50km + hint50_120km)*1.0e5; // 1.0e-5 acounts for km->cm;
//
		column_amount_mod = nmolec*1.0e6/nmolec_all;
		printf("\n MODTRAN column amount: %8.2e (ppmv)", column_amount_mod);
		if (column_amount_usr > 0.0)
		{
			printf("\n user column amount: %8.2e (ppmv)", column_amount_usr);
			scalef = column_amount_usr/column_amount_mod;
			for (iz = 0; iz < nz_mod; iz++) conc_cm3[iz] *= scalef;
			hint0_24km = simpson(conc_cm3, 25, 1.0);
//			hint24_25km = 0.5*(conc_cm3[24] + conc_cm3[25]);
			hint24_25km = intparab(zkm_mod[24], zkm_mod[25], zkm_mod[23], zkm_mod[24], zkm_mod[25],
			                                         conc_cm3[23], conc_cm3[24], conc_cm3[25]);
			hint25_50km = simpson(&conc_cm3[25], 11, 2.5);
			hint50_120km = simpson(&conc_cm3[35], 15, 5.0);
			nmolec = (hint0_24km + hint24_25km + hint25_50km + hint50_120km)*1.0e5; // 1.0e-5 for km->cm;
			column_amount_mod = nmolec*1.0e6/nmolec_all; // 1e6 for ppmv
			printf("\n check column amount: %8.2e (ppmv)", column_amount_mod);
		} // if (column_amount_usr > 0.0)
	} // else if molec_id == 1
//
//  read hitran
	nu_hit_min = nu_usr_min - delta_nu;
	nu_hit_max = nu_usr_max + delta_nu;
	printf("\n hitran spectral range: %9.2f - %9.2f (cm-1)", nu_hit_min, nu_hit_max);
//
//  run subroutine that takes hitran on input and spits k on output
	t1_sec = (double)clock() /(double)CLOCKS_PER_SEC;
	count_lines(molec_id, nu_hit_min, nu_hit_max, iline0, nlines);
	printf("\n iline0 = %i, nlines = %d", iline0, nlines);
	t2_sec = (double)clock() /(double)CLOCKS_PER_SEC;
	printf("\n runtime for 'count_lines': %6.2fs", t2_sec - t1_sec);
//
	isotop_id = new int[nlines];
	nuij = new double[nlines];
	Sij = new double[nlines];
	gamma_air = new double[nlines];
	gamma_self = new double[nlines];
	Epp = new double[nlines];
	n_air = new double[nlines];
	delta_air = new double[nlines];
//
	t1_sec = (double)clock() /(double)CLOCKS_PER_SEC;
	read_hitran160(molec_id, iline0, nlines, isotop_id, nuij, Sij,
				       gamma_air, gamma_self, Epp, n_air, delta_air);
	t2_sec = (double)clock() /(double)CLOCKS_PER_SEC;
	printf("\n runtime for 'read_hitran160': %6.2fs", t2_sec - t1_sec);
//
	t1_sec = (double)clock() /(double)CLOCKS_PER_SEC;
	kabs(molec_id, nu, nnu, nz_mod, Tkelv, Patm, Pgas, nlines, isotop_id, nuij, Sij, gamma_air,
	          gamma_self, Epp, n_air, delta_air, knu);
	t2_sec = (double)clock() /(double)CLOCKS_PER_SEC;
	printf("\n runtime for 'kabs': %6.2fs", t2_sec - t1_sec);
//
//  add continuum for H2O(2) or CO2(2) or NO2(10)
	if (molec_id == 1 || molec_id == 2 || molec_id == 10)
	{
		kcont = new double [nnu*nz_mod];
		for (ik = 0; ik < nnu*nz_mod; ik++) kcont[ik] = 0.0;
	    t1_sec = (double)clock() /(double)CLOCKS_PER_SEC;
	    cabs(molec_id, nu, nnu, Tkelv, Patm, conc_cm3, nz_mod, kcont);
	    t2_sec = (double)clock() /(double)CLOCKS_PER_SEC;
	    printf("\n runtime for 'cabs': %6.2fs", t2_sec - t1_sec);
		for (ik = 0; ik < nnu*nz_mod; ik++) knu[ik] = kcont[ik]; // ONLY continuum
		//for (ik = 0; ik < nnu*nz_mod; ik++) knu[ik] += kcont[ik]; // add continuum => check units!
	} // if molec_id with continuum
//
	t1_sec = (double)clock() /(double)CLOCKS_PER_SEC;
	pFile = fopen(fname_out, "w");
	fprintf(pFile, "# nu, zkm = ");
	for (iz = 0; iz < ntau; iz++) fprintf(pFile, "%8.3f", ztau[iz]);
	ext_km = new double [nz_mod];
	for (inu = 0; inu < nnu; inu++)
	{
		fprintf(pFile, "\n%8.2f", nu[inu]);
		for (iz = 0; iz < nz_mod; iz++)
			ext_km[iz] = knu[inu*nz_mod+iz]*conc_cm3[iz]*1.0e5; //cm-1 = 1e5.km-1
		tauabs25(ext_km, zkm_mod, ztau, ntau, tau_abs);
		for (iz = 0; iz < ntau; iz++) fprintf(pFile, "%16.6e", tau_abs[iz]);
	}
    fclose(pFile);
	t2_sec = (double)clock() /(double)CLOCKS_PER_SEC;
	printf("\n runtime for 'tauabs25' + write to file: %6.2fs", t2_sec - t1_sec);
//
	time_end = (double)clock() /(double)CLOCKS_PER_SEC;
	printf("\n cpu total time %6.2fs", time_end - time_start);
	printf("\n file = %s\n", fname_out);
	system("pause");
}
/*------------------------------------------------------------------------------
yy/mm/dd - First created
------------------------------------------------------------------------------*/