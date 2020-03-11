#include <stdio.h>        /* printf */
#include "const_param.h"  /* c2_rad */
#include "continuum.h"    /* data for continuum */
#include <math.h>         /* floor, pow */
/*
prototype:
    void cabs(int const, double const *, int const,
	          double const *, double const *, double const *, int const,
	          double *);
*/
void cabs(int const molec_id, double const *nu, int const nnu,                                 // in
	      double const *Tkelv, double const *Patm, double const *conc_cm3, int const nz,
	      double *kcont) {                                                                    // out
/*--------------------------------------------------------------------------------------------------
TASK:
	To compute the continuum absorbtion crossection per molecule.
IN:
	molec_id   i       as in HITRAN: 1-H2O, 2-CO2, 10-NO2
	nu         d[nu]   wavenumber in [cm-1]
	nnu        i       len(nu)
	Tkelv      d[nz]   temperature profile in [K]
	Patm       d[nz]   pressure profile in [atm]
	conc_cm3   d[nz]   volume number concentration [molec/cm3]
	nz         i       len(Tkelv)=len(Patm)=len(conc_cm3)
OUT:
    kcont      d[nz*nnu]   absorption coefficient [cm2/molec]
NOTE:
    On input, knu[nz*nnu] must be initilaized with 0.0.
	In knu[nz*nnu], iz=0:nz-1 runs frst (inner loop).
	See references in [1] for FASCODE and other information.
REFS:
	1. Clough SA, Kneizys FX, and Davies RW, 1989: Atm.Res. 23, p229
--------------------------------------------------------------------------------------------------*/
	bool
		new_dv_no2;
	int
		inu, inu1, inu2, iv, iz;
	double
		v_beg, v_end, dv, nui, c, vi, e, rad_term, xsect, Tz, r_no2, a_no2, b_no2, c_no2, d_no2,
		f_h2o, s_h2o, s296_h2o, s260_h2o, r_h2o, tpow, c_h2o, c_stp, atmcm_cm;
//--------------------------------------------------------------------------------------------------
//
	printf("\n(in cabs.cpp)                                                            CONTINUUM:");
	printf("\n(in cabs.cpp) user grid: nu = %10.3f : %10.3f cm-1, nnu = %i", nu[0], nu[nnu-1], nnu);
//
	if (molec_id == 1) // H2O
	{
		v_beg = VBeg_H2O;
		v_end = VEnd_H2O;
		dv = DV_H2O;
		printf("\n(in cabs.cpp) H2O continuum: v_beg = %10.3f, v_end = %10.3f, dv = %10.3f", v_beg, v_end, dv);
		if	(nu[0] > v_end || nu[nnu-1] < v_beg)
			printf("\n(in cabs.cpp) WARNING continuum & user grids do not overlap, kcont[:] = 0.0");
		else
		{
//          range of inu matching continuum grid. general case: dnu /= const
		    inu1 = 0;
		    nui = nu[0];
		    while (nui < v_beg && inu1 < nnu)
		    {
			    inu1 += 1;
			    nui = nu[inu1];
		    }
		    inu2 = nnu-1;
		    nui = nu[inu2];
		    while (nui > v_end)
		    {
			    inu2 -= 1;
			    nui = nu[inu2];
		    }
		    inu2 += 1; // in if: inu < inu2
//
		    printf("\n(in cabs.cpp) range of nu within continuum: nu[%i] = %10.3f : nu[%i] = %10.3f", inu1, nu[inu1], inu2-1, nu[inu2-1]);
//
		    for (inu = inu1; inu < inu2; inu++)
		    {
			    nui = nu[inu];
			    iv = (int) (nui - v_beg)/dv;
			    vi = v_beg + iv*dv;
				r_h2o = (nui - vi)/dv;
				f_h2o = FH2O_296[iv] + (FH2O_296[iv+1] - FH2O_296[iv])*r_h2o;
				s296_h2o = SH2O_296[iv] + (SH2O_296[iv+1] - SH2O_296[iv])*r_h2o;
				s260_h2o = SH2O_260[iv] + (SH2O_260[iv+1] - SH2O_260[iv])*r_h2o;
			    for (iz = 0; iz < nz; iz++)
			    {
					tpow = (Tkelv[iz] - T_ref)/(Th2o_lo - T_ref); 
				    s_h2o = s296_h2o*pow(s260_h2o/s260_h2o, tpow);
//                  see c_stp in main_tau_abs_hprofile.cpp for clarification
				    c_stp = (k_boltzman*1.0e-7)*To_stp/Po_stp*1.0e6;
				    atmcm_cm = c_stp*conc_cm3[iz];
				    c_h2o = s_h2o*atmcm_cm + f_h2o*(Patm[iz]*To_stp/Tkelv[iz]); // [atm] is already a relative unit
					e = exp(-c2_rad*nui/Tkelv[iz]);
				    rad_term = nui*(1.0 - e)/(1.0 + e); // [1: Eq.(2)]
					xsect = c_h2o*rad_term; // [kcont]=[cm3/molec]*[1/cm]=[xsect]=cm2/molec
				    kcont[inu*nz+iz] = xsect*1.0e-20; // 1.0e-20 accounts for scalef in BOTH F&S continuums
			    } // for iz
		    } // for inu
		} // if	(nu[0] > v_end || nu[nnu-1] < v_beg) ... else ... - check for overlap
	} // if molec_id == 1
	else if (molec_id == 2) // CO2
	{
		v_beg = VBeg_CO2;
		v_end = VEnd_CO2;
		dv = DV_CO2;
		printf("\n(in cabs.cpp) CO2 continuum: v_beg = %10.3f, v_end = %10.3f, dv = %10.3f", v_beg, v_end, dv);
		if	(nu[0] > v_end || nu[nnu-1] < v_beg)
			printf("\n(in cabs.cpp) WARNING continuum & user grids do not overlap");
		else
		{
//          range of inu matching continuum grid. general case: dnu /= const
		    inu1 = 0;
		    nui = nu[0];
		    while (nui < v_beg && inu1 < nnu)
		    {
			    inu1 += 1;
			    nui = nu[inu1];
		    }
		    inu2 = nnu-1;
		    nui = nu[inu2];
		    while (nui > v_end)
		    {
			    inu2 -= 1;
			    nui = nu[inu2];
		    }
		    inu2 += 1; // in if: inu < inu2
//
		    printf("\n(in cabs.cpp) range of nu within continuum: nu[%i] = %10.3f : nu[%i] = %10.3f", inu1, nu[inu1], inu2-1, nu[inu2-1]);
//
		    for (inu = inu1; inu < inu2; inu++)
		    {
			    nui = nu[inu];
			    iv = (int) (nui - v_beg)/dv;
			    vi = v_beg + iv*dv;
			    c = FCO2[iv] + (FCO2[iv+1] - FCO2[iv])*(nui - vi)/dv;
			    for (iz = 0; iz < nz; iz++)
			    {
				    e = exp(-c2_rad*nui/Tkelv[iz]);
				    rad_term = nui*(1.0 - e)/(1.0 + e); // [1: Eq.(2)]
				    xsect = c*rad_term; // cm2/molec = cm3/molec * cm-1
//                  [kcont] = [xsect] = cm2/molec
				    kcont[inu*nz+iz] = xsect*1.0e-20; // 1.0e-20 accounts for scalef in FCO2/1.0e-20
			    } // for iz
		    } // for inu
		} // if	(nu[0] > v_end || nu[nnu-1] < v_beg) ... else ... - check for overlap
	} // else if molec_id == 2
	else if (molec_id == 10) // NO2
	{
		v_beg = VBeg_NO2;
		v_end = VEnd_NO2;
		printf("\n(in cabs.cpp) NO2 continuum: v_beg = %10.3f, v_end = %10.3f", v_beg, v_end);
		if	(nu[0] > v_end || nu[nnu-1] < v_beg)
			printf("\n(in cabs.cpp) WARNING continuum & user grid do not overlap");
		else
		{
//          range of inu matching continuum grid. general case: dnu /= const
		    inu1 = 0;
		    nui = nu[0];
		    while (nui < v_beg && inu1 < nnu-1)
		    {
			    inu1 += 1;
			    nui = nu[inu1];
		    }
		    inu2 = nnu-1;
		    nui = nu[inu2];
		    while (nui > v_end)
		    {
			    inu2 -= 1;
			    nui = nu[inu2];
		    }
		    inu2 += 1; // in if: inu < inu2
//
		    printf("\n(in cabs.cpp) range of nu within continuum: nu[%i] = %10.3f : nu[%i] = %10.3f", inu1, nu[inu1], inu2-1, nu[inu2-1]);
//
//          find iv - continuum grid index to the left from inu1
		    iv = 0;
		    vi = V_NO2[iv];
		    nui = nu[inu1];
		    while (vi < nui)
		    {
			    iv += 1;
			    vi = V_NO2[iv];
		    }
		    iv -= 1;
		    new_dv_no2 = true;
		    for (inu = inu1; inu < inu2; inu++)
		    {
			    if (new_dv_no2)
			    {
                    r_no2 = (nu[inu] - V_NO2[iv])/(V_NO2[iv+1] - V_NO2[iv]); 
			        a_no2 = sigNO2_A[iv] + (sigNO2_A[iv+1] - sigNO2_A[iv])*r_no2; 
			        b_no2 = sigNO2_B[iv] + (sigNO2_B[iv+1] - sigNO2_B[iv])*r_no2;
			        c_no2 = sigNO2_C[iv] + (sigNO2_C[iv+1] - sigNO2_C[iv])*r_no2;
			        d_no2 = sigNO2_D[iv] + (sigNO2_D[iv+1] - sigNO2_D[iv])*r_no2;
				    new_dv_no2 = false;
			    } // if new interval v[iv] ... v[iv+1]
//
			    for (iz = 0; iz < nz; iz++)
			    {
				    Tz = Tkelv[iz];   
				    kcont[inu*nz+iz] = a_no2 + (b_no2 + (c_no2 + d_no2*Tz)*Tz)*Tz; // [kcont] = [xsect] = cm2/molec
			    } // for iz
//
//              recompute a,b,c,d_no2 if next inu exists & does *not* belong to V_NO2[iv]:V_NO2[iv+1]
			    if (inu+1 < inu2 && nu[inu+1] > V_NO2[iv+1])
				    new_dv_no2 = true;		
		    } // for inu
		} // if	(nu[0] > v_end || nu[nnu-1] < v_beg) ... else ... - check for overlap
	} // else if molec_id == 10 (NO2)
} // cabs
//--------------------------------------------------------------------------------------------------
/*
2020-02-21:
    First created (CO2), compiled, plotted. No precise testing yet.
*/