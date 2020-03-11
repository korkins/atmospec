#include <stdio.h>  /* printf, fprintf */
#include <iostream> /* system("pause") */
#include "Profile_GasAbs.h"
//
int main()
{
	int const iatm = 6;
	int ix_atm, iz;
	double c_total, c_h2o, c_co2, c_o3, c_n2o, c_co, c_ch4, c_o2, c_no2;
//
	ix_atm = iatm-1;
	for (iz = 0; iz < N_zLevels; iz++)
	{
		c_h2o = cH2O[ix_atm][iz];
		c_co2 = cCO2[ix_atm][iz];
		c_o3 = cO3[ix_atm][iz];
		c_n2o = cN2O[ix_atm][iz];
		c_co = cCO[ix_atm][iz];
		c_ch4 = cCH4[ix_atm][iz];
		c_o2 = cO2[ix_atm][iz];
		c_no2 = cNO2[ix_atm][iz];
		c_total = c_h2o + c_co2 + c_o3 + c_n2o + c_co + c_ch4 + c_o2 + c_no2;
		printf("\n%5.1f %12.4e", z0[iz], c_ch4);
	}
	printf("\nDone!\n");
	system("pause");
	return 0;
}