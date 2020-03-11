#include <stdio.h>  /* printf, fprintf */
#include <iostream> /* system("pause") */
#include <math.h>
//
double intparab(double const, double const, double const, double const,
	                          double const, double const, double const);
//
int main()
{
	double a, b, c, x1, x2, x3, y1, y2, y3, x, ipar;
	a = 3.0;
	b = 1.5;
	c = 2.0;
	x1 = 1.0;
	x2 = 3.0;
	x3 = 4.0;
	x = (x1 + x2)/2.0;
	y1 = a*x1*x1 + b*x1 + c;
	y2 = a*x2*x2 + b*x2 + c;
	y3 = a*x3*x3 + b*x3 + c;
	ipar = intparab(x, x1, x2, x3, y1, y2, y3);
    printf("\n iparab = %10.6f", ipar);
	printf("\nDone!\n");
	system("pause");
	return 0;
}