/*****************************************************************

Basic closed-shell spin-restricted DFT-solver for simple molecules using STO-NG

Authors: B. Klumpers
         I.A.W. Filot

Published under GNU General Public License 3.0
Source code available at: https://github.com/BKlumpers/dft

Allows for SCF-computation of molecular energies for simple molecules.
Includes testcases for: H, He, H2, HeH+, He2, CO, and H2O.

*****************************************************************/

#ifndef _ORBITALS_H
#define _ORBITALS_H
//Header-only

#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "prelim_math.h"

typedef Eigen::Vector3d vec3; //define vec3 as the Eigen::Vector3d-object

using namespace std;

//Gaussian type orbital
class GTO {
private:
    double norm = 0.0;            //storage for the normalisation constant
public:
    double alpha;           //exponent factor
    int a,b,c;              //angular momentum components in x,y,z respectively
    vec3 atompos;           //atomic coordinates of the nucleus
    double getnorm()
    {
        if(norm==0.0){
            static const double pi = 3.14159265359;
            //compute normalisation constant, only if norm was not computed previously
            norm = sqrt(pow(2.0*alpha/pi, 1.5)*pow(4.0*alpha, a + b + c)/(doublefactorial_odd(a)*doublefactorial_odd(b)*doublefactorial_odd(c)));
        }
        return norm;
    }
    double getvalue(const vec3& pos)
    {
        //get amplitude of GTO-wavefunction at r = pos
        return getnorm()*pow(pos[0] - atompos[0], a)*pow(pos[1] - atompos[1], b)*pow(pos[2]-atompos[2], c)*exp(-alpha*(pos-atompos).squaredNorm());
    }
};

//Contracted Gaussian orbital
class CGF {
public:
    vector<double> coeff; //fix variable number of GTOs per CGF
    vector<GTO> GTOlist; //list of GTOs
    //no direct normalisation, assumed CGF of normalised GTO (hence CGF not inherently normalised)
    void add_gto(double gto_coeff, double gto_alpha, int gto_a, int gto_b, int gto_c, const vec3& gto_atompos)
    {
        GTO gto_new;
        gto_new.alpha = gto_alpha;
        gto_new.a = gto_a;
        gto_new.b = gto_b;
        gto_new.c = gto_c;
        gto_new.atompos << gto_atompos;
        GTOlist.push_back(gto_new); //add GTO to list
        coeff.push_back(gto_coeff); //add coefficient to list
    }
    int get_size(){
        return coeff.size();
    }
    double getvalue(const vec3& pos)
    {
        //get expectation value of CGF at r = pos
        double out = 0.0;
        for(int i=0; i<get_size(); i++)
        {
            out = out + coeff[i]*GTOlist[i].getvalue(pos);
        }
        return out;
    }
};

#endif //_ORBITALS_H