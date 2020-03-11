#include <stdio.h>  /* printf, fprintf */
#include <iostream> /* system("pause") */
#include <math.h>   /* pow */
#include <time.h>
#include "const_param.h"
//
double simpson(double const *, int const, double const);
void lininterp(int const, double const *, double const *, int const, double const *, double *);
//
int main()
{
	FILE *fin;
	char
		fpath[path_len_max];
	int
		nnu, inu, iline, nlines, nnu_src, nnu_bpf;
	double
		   nu1, nu2, s1, s2, f1, f2, t1, t2,
		   time_start, time_end,
		   nu_src_min, nu_src_max, dnu_src,
		   nu_bpf_min, nu_bpf_max, dnu_bpf,
		   nu_min, nu_max, dnu,
		   norm, nu_center, trans;
	double
		*nu, *s0, *bpf, *tau, *nu_src, *s0_src, *nu_bpf, *f_bpf, *s0_bpf, *s0_bpf_nu, *s0_bpf_et;
/*------------------------------------------------------------------------------------------------*/
//
//  Input:
	dnu = 0.1;
/*------------------------------------------------------------------------------------------------*/
//
	time_start = (double)clock() /(double)CLOCKS_PER_SEC;
//  Echo some parameters
	printf("PROJECT = prj_transmit:");
	printf("\n-delta_nu=%6.2f [cm-1]", delta_nu);
//
//  BPF
	strcpy(fpath, path_bpf);
	strcat(fpath, fname_bpf);
	fin = fopen(fpath, "r");
	printf("\n(in main) opened: %s", fpath);
	fscanf(fin, "%i %lf %lf", &nnu_bpf, &nu_bpf_min, &dnu_bpf);
	nnu_bpf += 1;
	nu_bpf = new double [nnu_bpf];
	f_bpf = new double [nnu_bpf];
	for (inu = 0; inu < nnu_bpf; inu++)
	{
		nu_bpf[inu] = nu_bpf_min + inu*dnu_bpf;
		f_bpf[inu] = 0.0;
	}
	nlines = nnu_bpf/5;
	inu = 0;
	for (iline = 0; iline < nlines; iline++)
	{
		fscanf(fin, "%lf %lf %lf %lf %lf", &f_bpf[inu], &f_bpf[inu+1], &f_bpf[inu+2], &f_bpf[inu+3], &f_bpf[inu+4]);
		inu += 5;
	}
	fclose(fin);
	printf("\n(in main) closed: %s", fpath);
//
//  USER Grid
	nu_min = nu_bpf[0];
	nu_max = nu_bpf[nnu_bpf-1];
	nnu = ceil((nu_max - nu_min)/dnu); // must be ODD for Simpson's rule
	nu = new double[nnu];
	s0 = new double[nnu];
	bpf = new double[nnu];
	for (inu = 0; inu < nnu; inu++)
	{
		nu[inu] = nu_min + inu*dnu;
		s0[inu] = 0.0;
		bpf[inu] = 0.0;
	}
	lininterp(nnu_bpf, nu_bpf, f_bpf, nnu, nu, bpf);
//
//  S0:
	strcpy(fpath, path_source);
	strcat(fpath, fname_source);
	fin = fopen(fpath, "r");
	printf("\n(in main) opened: %s", fpath);
	fscanf(fin, "%i", &nlines);
	nu_src = new double[nlines];
	s0_src = new double[nlines];
	fscanf(fin, "%lf %lf", &nu1, &s1);
	iline = 0;
	fscanf(fin, "%lf %lf", &nu2, &s2);
	iline = 1;
	while (nu2 < nu_min && iline < nlines)
	{
		nu1 = nu2;
		s1 = s2;
		fscanf(fin, "%lf %lf", &nu2, &s2);
		iline += 1;
	}
	nu_src[0] = nu1;
	s0_src[0] = s1;
	nu_src[1] = nu2;
	s0_src[1] = s2;
	inu = 1;
	while (nu_src[inu] < nu_max && iline < nlines)
	{
		inu += 1;
		fscanf(fin, "%lf %lf", &nu_src[inu], &s0_src[inu]);
		iline += 1;
	}
	iline += 1;
	if (iline <= nlines)
	{
		inu += 1;
		fscanf(fin, "%lf %lf", &nu_src[inu], &s0_src[inu]);
	}
	nnu_src = inu+1;
	lininterp(nnu_src, nu_src, s0_src, nnu, nu, s0);
	fclose(fin);
	printf("\n(in main) closed: %s", fpath);
//
//  tau_abs:
	strcpy(fpath, path_tau);
	strcat(fpath, fname_tau);
	fin = fopen(fpath, "r");
	printf("\n(in main) opened: %s", fpath);
	fscanf(fin, "%i", &nlines);
	tau = new double[nnu];
	fscanf(fin, "%lf %lf", &nu1, &t1);
	fscanf(fin, "%lf %lf", &nu2, &t2);
	iline = 1;
	while (nu2 < nu[0] && iline < nlines) // compare as |nu2-nu0| > tiny
	{
		nu1 = nu2;
		t1 = t2;
		fscanf(fin, "%lf %lf", &nu2, &t2);
		iline += 1;
	}
	tau[0] = t1;
	tau[1] = t2;
	for (inu = 2; inu < nnu; inu++)
		fscanf(fin, "%lf %lf", &nu2, &tau[inu]);
//
	s0_bpf = new double[nnu];
	s0_bpf_nu = new double[nnu];
	s0_bpf_et = new double[nnu];
	for (inu = 0; inu < nnu; inu++)
	{
		s0_bpf[inu] = s0[inu]*bpf[inu];
	    s0_bpf_nu[inu] = s0[inu]*bpf[inu]*nu[inu];
		s0_bpf_et[inu] = s0[inu]*bpf[inu]*exp(-tau[inu]);
	}
	norm = simpson(s0_bpf, nnu, dnu);
	nu_center = simpson(s0_bpf_nu, nnu, dnu)/norm;
	printf("\ncenter wavenumber: %10.2f [1/cm]", nu_center);
	trans = simpson(s0_bpf_et, nnu, dnu)/norm;
	printf("\ntransmittance: %6.4f", trans);
//
	time_end = (double)clock() /(double)CLOCKS_PER_SEC;
	printf("\ncpu total time %6.2fs\n", time_end - time_start);
	system("\npause");
}
/*------------------------------------------------------------------------------
yy/mm/dd - First created
------------------------------------------------------------------------------*/