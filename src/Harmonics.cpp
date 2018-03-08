/*****************************************************************

Basic closed-shell spin-restricted DFT-solver for simple molecules using STO-NG

Authors: B. Klumpers
		 I.A.W. Filot

Published under GNU General Public License 3.0
Source code available at: https://github.com/BKlumpers/dft

Allows for SCF-computation of molecular energies for simple molecules.
Includes testcases for: H, He, H2, HeH+, He2, CO, and H2O.

*****************************************************************/

#include "Harmonics.h"

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
		out = root*pow(-1.0,m)*boost::math::spherical_harmonic_r(ul, m, theta, phi);
	}
	else{
		out = root*pow(-1.0,m)*boost::math::spherical_harmonic_i(ul, -m, theta, phi);
	}
	if(isnan(out)){
		cout << "nanerr: " << radius << "   " << theta << "   " << phi << "   " << gridpoint[1] << "   " << gridpoint[0] << "   " << gridpoint[1]/gridpoint[0] << endl;
	}
	return out;
}

//compute real spherical harmonics up to degree 8 in Cartesian space
//
//{l} is the degree of the real spherical harmonic (max 8)
//{m} is the order of the l-th degree real spherical harmonic
//{gridpoint} is the vector of {x,y,z} coordinates at which the real spherical harmonic is to be evaluated
//
//write all coefficients as precomputed decimals, give expressions in comments
//
// double HarmonicReal(int l, int m, vec3 gridpoint){
// 	static double radius;
// 	radius = gridpoint.norm();

// 	switch(l){
// 		case 0:
// 		switch(m){
// 			case 0:
// 			return 0.28209479177; //coeff = 1/sqrt(4*pi);
// 			break;
// 			default:
// 			cout << "Encountered invalid index for real spherical harmonic function: |m| > l." << endl;
// 			cout << "For the " << l << "th degree harmonic, the requested order was " << m << "." << endl;
// 			cout << "Defaulting to 0th order;" << endl;
// 			return 0.28209479177;
// 			break;
// 		}
// 		break; //end 0th degree
// 		case 1:
// 		switch(m){
// 			case -1:
// 			return 0.48860251190 * gridpoint[1] / radius; //coeff = sqrt(3.0/(4.0*pi));
// 			break;
// 			case 0:
// 			return 0.48860251190 * gridpoint[2] / radius; //coeff = sqrt(3.0/(4.0*pi));
// 			break;
// 			case 1:
// 			return 0.48860251190 * gridpoint[0] / radius; //coeff = sqrt(3.0/(4.0*pi));
// 			break;
// 			default:
// 			cout << "Encountered invalid index for real spherical harmonic function: |m| > l." << endl;
// 			cout << "For the " << l << "th degree harmonic, the requested order was " << m << "." << endl;
// 			cout << "Defaulting to 0th order." << endl;
// 			return 0.48860251190 * gridpoint[2] / radius;
// 			break;
// 		}
// 		break; //end 1st degree
// 		case 2:
// 		switch(m){
// 			case -2:
// 			return 1.09254843059 * gridpoint[0] * gridpoint[1] / (radius*radius); //coeff = sqrt(15.0/(4.0*pi));
// 			break;
// 			case -1:
// 			return 1.09254843059 * gridpoint[1] * gridpoint[2] / (radius*radius); //coeff = sqrt(15.0/(4.0*pi));
// 			break;
// 			case 0:
// 			return 0.31539156525 * (2*gridpoint[2]*gridpoint[2] - gridpoint[0]*gridpoint[0] - gridpoint[1]*gridpoint[1]) / (radius*radius); //coeff = sqrt(5.0/pi)/4;
// 			break;
// 			case 1:
// 			return 1.09254843059 * gridpoint[0] * gridpoint[2] / (radius*radius); //coeff = sqrt(15.0/(4.0*pi));
// 			break;
// 			case 2:
// 			return 0.31539156525 * (gridpoint[0]*gridpoint[0] - gridpoint[1]*gridpoint[1]) / (radius*radius); //coeff = sqrt(5.0/pi)/4;
// 			break;
// 			default:
// 			cout << "Encountered invalid index for real spherical harmonic function: |m| > l." << endl;
// 			cout << "For the " << l << "th degree harmonic, the requested order was " << m << "." << endl;
// 			cout << "Defaulting to 0th order." << endl;
// 			return 0.31539156525 * (2*gridpoint[2]*gridpoint[2] - gridpoint[0]*gridpoint[0] - gridpoint[1]*gridpoint[1]) / (radius*radius);
// 			break;
// 		}
// 		break; //end 2nd degree
// 		case 3:
// 		switch(m){
// 			case -3:
// 			return 0.59004358993 * (3.0*gridpoint[0]*gridpoint[0] - gridpoint[1]*gridpoint[1]) * gridpoint[1] / (radius*radius*radius); //coeff = sqrt(35.0/(2.0*pi))/4;
// 			break;
// 			case -2:
// 			return 2.89061144264 * gridpoint[0] * gridpoint[1] * gridpoint[2] / (radius*radius*radius); //coeff = sqrt(105.0/pi)/2;
// 			break;
// 			case -1:
// 			return 0.45704579946 * (4.0*gridpoint[2]*gridpoint[2] - gridpoint[0]*gridpoint[0] - gridpoint[1]*gridpoint[1]) * gridpoint[1] / (radius*radius*radius); //coeff = sqrt(21.0/(2.0*pi))/4;
// 			break;
// 			case 0:
// 			return 0.37317633259 * (2.0*gridpoint[2]*gridpoint[2] - 3.0*gridpoint[0]*gridpoint[0] - 3.0*gridpoint[1]*gridpoint[1]) * gridpoint[2] / (radius*radius*radius); //coeff = sqrt(7/pi)/4;
// 			break;
// 			case 1:
// 			return 0.45704579946 * (4.0*gridpoint[2]*gridpoint[2] - gridpoint[0]*gridpoint[0] - gridpoint[1]*gridpoint[1]) * gridpoint[0] / (radius*radius*radius); //coeff = sqrt(21.0/(2.0*pi))/4;
// 			break;
// 			case 2:
// 			return 1.44530572132 * (gridpoint[0]*gridpoint[0] - gridpoint[1]*gridpoint[1]) * gridpoint[2] / (radius*radius*radius); //coeff = sqrt(105.0/pi)/2;
// 			break;
// 			case 3:
// 			return 0.59004358993 * (gridpoint[0]*gridpoint[0] - 3.0*gridpoint[1]*gridpoint[1]) * gridpoint[0] / (radius*radius*radius); //coeff = sqrt(35.0/(2.0*pi))/4;
// 			break;
// 			default:
// 			cout << "Encountered invalid index for real spherical harmonic function: |m| > l." << endl;
// 			cout << "For the " << l << "th degree harmonic, the requested order was " << m << "." << endl;
// 			cout << "Defaulting to 0th order." << endl;
// 			return 0.37317633259 * (2.0*gridpoint[2]*gridpoint[2] - 3.0*gridpoint[0]*gridpoint[0] - 3.0*gridpoint[1]*gridpoint[1]) * gridpoint[2] / (radius*radius*radius);
// 			break;
// 		}
// 		break; //end 3rd degree
// 		case 4:
// 		switch(m){
// 			case -4:
// 			return 1.09254843059 * gridpoint[0] * gridpoint[1] / (radius*radius); //coeff = 
// 			break;
// 			case -3:
// 			return 1.09254843059 * gridpoint[0] * gridpoint[1] / (radius*radius); //coeff = 
// 			break;
// 			case -2:
// 			return 1.09254843059 * gridpoint[0] * gridpoint[1] / (radius*radius); //coeff = 
// 			break;
// 			case -1:
// 			return 1.09254843059 * gridpoint[1] * gridpoint[2] / (radius*radius); //coeff = 
// 			break;
// 			case 0:
// 			return 0.31539156525 * (2*gridpoint[2]*gridpoint[2] - gridpoint[0]*gridpoint[0] - gridpoint[1]*gridpoint[1]) / radius; //coeff = 
// 			break;
// 			case 1:
// 			return 1.09254843059 * gridpoint[0] * gridpoint[2] / (radius*radius); //coeff = 
// 			break;
// 			case 2:
// 			return 0.31539156525 * (gridpoint[0]*gridpoint[0] - gridpoint[1]*gridpoint[1]) / radius; //coeff = 
// 			break;
// 			case -3:
// 			return 1.09254843059 * gridpoint[0] * gridpoint[1] / (radius*radius); //coeff = 
// 			break;
// 			case 4:
// 			return 1.09254843059 * gridpoint[0] * gridpoint[1] / (radius*radius); //coeff = 
// 			break;
// 			default:
// 			cout << "Encountered invalid index for real spherical harmonic function: |m| > l." << endl;
// 			cout << "For the " << l << "th degree harmonic, the requested order was " << m << "." << endl;
// 			cout << "Defaulting to 0th order." << endl;
// 			return 0.31539156525 * (2*gridpoint[2]*gridpoint[2] - gridpoint[0]*gridpoint[0] - gridpoint[1]*gridpoint[1]) / radius;
// 			break;
// 		}
// 		break; //end 4th degree
// 		case 5:
		
// 		break; //end 5th degree
// 		case 6:
		
// 		break; //end 6th degree
// 		case 7:
		
// 		break; //end 7th degree
// 		case 8:
		
// 		break; //end 8th degree
// 		default:
// 		if(l<0){
// 			cout << "Encountered invalid index for real spherical harmonic function: l < 0" << endl;
// 		}
// 		else{
// 			cout << "Encountered invalid index for real spherical harmonic function: l > known degrees" << endl;
// 		}
// 		cout << "Requested degree was: " << l << endl;
// 		cout << "Defaulting to 0th order." << endl;
// 		return 0.28209479177;
// 		break;
// 	}
// }

