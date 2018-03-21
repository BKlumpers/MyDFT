/*****************************************************************

MYDFT: Basic DFT-solver for simple molecules using STO-NG

Authors: B. Klumpers
         I.A.W. Filot

Published under GNU General Public License 3.0
Source code available at: https://github.com/BKlumpers/dft

Allows for SCF-computation of molecular energies for simple molecules
using both HF and DFT. Supports both closed and open shell systems.
Includes testcases for: H, He, H2, HeH+, HeH, He2, CO, and H2O.

*****************************************************************/

#include "Harmonics.h"

//compute tesseral harmonics (real spherical harmonics)
//
//{l} is the degree of the real spherical harmonic
//{m} is the order of the l-th degree real spherical harmonic
//{gridpoint} is the vector of {x,y,z} coordinates at which the real spherical harmonic is to be evaluated
//
double HarmonicReal(int l, int m, const vec3& gridpoint){
	static double radius, theta, phi;
	static unsigned int ul;
	static const double root = sqrt(2.0);
	static double out;
	radius = gridpoint.norm();
	theta = acos(gridpoint[2]/radius);
	phi = atan2(gridpoint[1],gridpoint[0]);
	ul = l; //let g++ typecasting perform conversion from int to unsigned int
	if(m == 0){
		out = boost::math::spherical_harmonic_r(ul, m, theta, phi); //take real component since imaginary part should be zero, avoid typecast complex
	}
	else if(m > 0){
		out = root*boost::math::spherical_harmonic_r(ul, m, theta, phi); //*pow(-1.0,m)
	}
	else{
		out = root*boost::math::spherical_harmonic_i(ul, -m, theta, phi); //*pow(-1.0,m)
	}
	if(isnan(out)){
		//error message in case of unpredicted behaviour
		cout << "nanerr: " << radius << "   " << theta << "   " << phi << "   " << gridpoint[1] << "   " << gridpoint[0] << "   " << gridpoint[1]/gridpoint[0] << endl;
	}
	return out;
}

//End of file
//****************************************************************
