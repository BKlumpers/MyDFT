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
    //get value of the derivative of the wavefunction w.r.t. x, y and z
    vec3 getderiv(const vec3& pos){
        static vec3 out;
        if(a == 0){
            out[0] = -2.0*alpha*getnorm()*pow(pos[1] - atompos[1], b)*pow(pos[2]-atompos[2], c)*exp(-alpha*(pos-atompos).squaredNorm()); //deriv to x
        }
        else{
            out[0] = getnorm()*(-2.0*alpha * pow(pos[0] - atompos[0], a+1) + a*pow(pos[0] - atompos[0], a-1)) *pow(pos[1] - atompos[1], b)*pow(pos[2]-atompos[2], c)*exp(-alpha*(pos-atompos).squaredNorm());
        }
        if(b == 0){
            out[1] = -2.0*alpha*getnorm()*pow(pos[0] - atompos[0], a)*pow(pos[2]-atompos[2], c)*exp(-alpha*(pos-atompos).squaredNorm()); //deriv to y
        }
        else{
            out[1] = getnorm()*(-2.0*alpha * pow(pos[1] - atompos[1], b+1) + b*pow(pos[1] - atompos[1], b-1)) *pow(pos[0] - atompos[0], a)*pow(pos[2]-atompos[2], c)*exp(-alpha*(pos-atompos).squaredNorm());
        }
        if(c == 0){
            out[2] = -2.0*alpha*getnorm()*pow(pos[1] - atompos[1], b)*pow(pos[0]-atompos[0], a)*exp(-alpha*(pos-atompos).squaredNorm()); //deriv to z
        }
        else{
            out[2] = getnorm()*(-2.0*alpha * pow(pos[2] - atompos[2], c+1) + c*pow(pos[2] - atompos[2], c-1)) *pow(pos[1] - atompos[1], b)*pow(pos[0]-atompos[0], c)*exp(-alpha*(pos-atompos).squaredNorm());
        }
        return out;
    }
    //get Hessian of the wavefunction:
    //
    //  |  xx  |  yx  |  zx  |
    //   ____________________
    //  |  xy  |  yy  |  zy  |
    //   ____________________
    //  |  xz  |  yz  |  zz  |
    //
    Eigen::MatrixXd getHessian(const vec3& pos){
        Eigen::MatrixXd out = Eigen::MatrixXd::Zero(3,3);
        double expo = getnorm()*exp(-alpha*(pos-atompos).squaredNorm());

        double factor = 4.0*alpha*alpha*pow(pos[0] - atompos[0], a+2) - 2*alpha*(2*a+1)*pow(pos[0] - atompos[0], a);
        if(a > 1){
            factor += a*(a+1)*pow(pos[0] - atompos[0], a-2);
        }
        out(0,0) = pow(pos[1] - atompos[1], b)*pow(pos[2]-atompos[2], c)*expo;
        
        factor = 4.0*alpha*alpha*pow(pos[1] - atompos[1], b+2) - 2*alpha*(2*b+1)*pow(pos[1] - atompos[1], b);
        if(b > 1){
            factor += b*(b+1)*pow(pos[1] - atompos[1], b-2);
        }
        out(1,1) = factor*pow(pos[0] - atompos[0], a)*pow(pos[2]-atompos[2], c)*expo;

        factor = 4.0*alpha*alpha*pow(pos[2] - atompos[2], c+2) - 2*alpha*(2*c+1)*pow(pos[2] - atompos[2], c);
        if(c > 1){
            factor += c*(c+1)*pow(pos[2] - atompos[2], c-2);
        }
        out(2,2) = factor*pow(pos[0] - atompos[0], a)*pow(pos[1]-atompos[1], b)*expo;
        
        factor = 4.0*alpha*alpha*pow(pos[0] - atompos[0], a+1)*pow(pos[1] - atompos[1], b+1);
        if(b > 0){
            factor -= 2*alpha*b*pow(pos[0] - atompos[0], a+1)*pow(pos[1] - atompos[1], b-1);
        }
        if(a > 0){
            factor -= 2*alpha*a*pow(pos[0] - atompos[0], a-1)*pow(pos[1] - atompos[1], b+1);
        }
        if(a > 0 && b > 0){
            factor += a*b*pow(pos[0] - atompos[0], a-1)*pow(pos[1] - atompos[1], b-1);
        }
        out(1,0) = out(0,1) = factor*pow(pos[2]-atompos[2], c)*expo;

        factor = 4.0*alpha*alpha*pow(pos[2] - atompos[2], c+1)*pow(pos[1] - atompos[1], b+1);
        if(b > 0){
            factor -= 2*alpha*b*pow(pos[2] - atompos[2], c+1)*pow(pos[1] - atompos[1], b-1);
        }
        if(c > 0){
            factor -= 2*alpha*c*pow(pos[2] - atompos[2], c-1)*pow(pos[1] - atompos[1], b+1);
        }
        if(c > 0 && b > 0){
            factor += c*b*pow(pos[2] - atompos[2], c-1)*pow(pos[1] - atompos[1], b-1);
        }
        out(2,1) = out(1,2) = factor*pow(pos[0]-atompos[0], a)*expo;

        factor = 4.0*alpha*alpha*pow(pos[0] - atompos[0], a+1)*pow(pos[2] - atompos[2], c+1);
        if(c > 0){
            factor -= 2*alpha*c*pow(pos[0] - atompos[0], a+1)*pow(pos[2] - atompos[2], c-1);
        }
        if(a > 0){
            factor -= 2*alpha*a*pow(pos[0] - atompos[0], a-1)*pow(pos[2] - atompos[2], c+1);
        }
        if(a > 0 && c > 0){
            factor += a*c*pow(pos[0] - atompos[0], a-1)*pow(pos[2] - atompos[2], c-1);
        }
        out(2,0) = out(0,2) = factor*pow(pos[1]-atompos[1], b)*expo;

        return out;

        // if(a == 0){
        //     out(0,0) = getnorm()*(  4.0*alpha*alpha*pow(pos[0] - atompos[0], a+2) - 2*alpha*pow(pos[0] - atompos[0], a)  )*pow(pos[1] - atompos[1], b)*pow(pos[2]-atompos[2], c)*expo;
                        
        //     if(b == 0){
        //         out(1,1) = getnorm()*(  4.0*alpha*alpha*pow(pos[1] - atompos[1], b+2) - 2*alpha*pow(pos[1] - atompos[1], b)  )*pow(pos[0] - atompos[0], a)*pow(pos[2]-atompos[2], c)*expo;
        //         out(1,0) = out(0,1) = getnorm()*(  4.0*alpha*alpha*pow(pos[0] - atompos[0], a+1)*pow(pos[1] - atompos[1], b+1) )*pow(pos[2]-atompos[2], c)*expo;
        //         if(c == 0){
        //             out(2,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[2] - atompos[2], c+2) - 2*alpha*pow(pos[2] - atompos[2], c)  )*pow(pos[1] - atompos[1], b)*pow(pos[0]-atompos[0], a)*expo;
        //             out(2,0) = out(0,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[0] - atompos[0], a+1)*pow(pos[2] - atompos[2], c+1) )*pow(pos[1]-atompos[1], b)*expo;
        //             out(2,1) = out(1,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[1] - atompos[1], b+1)*pow(pos[2] - atompos[2], c+1) )*pow(pos[0]-atompos[0], a)*expo;
        //         }
        //         else if(c == 1){
        //             out(2,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[2] - atompos[2], c+2) - 2*alpha*(2*c+1)*pow(pos[2] - atompos[2], c)  )*pow(pos[1] - atompos[1], b)*pow(pos[0]-atompos[0], a)*expo;
        //             out(2,0) = out(0,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[0] - atompos[0], a+1)*pow(pos[2] - atompos[2], c+1) - 2*alpha*c*pow(pos[0] - atompos[0], a+1)*pow(pos[2] - atompos[2], c-1) )*pow(pos[1]-atompos[1], b)*expo;
        //             out(2,1) = out(1,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[1] - atompos[1], b+1)*pow(pos[2] - atompos[2], c+1) - 2*alpha*c*pow(pos[1] - atompos[1], b+1)*pow(pos[2] - atompos[2], c-1) )*pow(pos[0]-atompos[0], a)*expo;
        //         }
        //         else{
        //             out(2,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[2] - atompos[2], c+2) - 2*alpha*(2*c+1)*pow(pos[2] - atompos[2], c) + c*(c+1)*pow(pos[2] - atompos[2], c-2)  )*pow(pos[1] - atompos[1], b)*pow(pos[0]-atompos[0], a)*expo;
        //             out(2,0) = out(0,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[0] - atompos[0], a+1)*pow(pos[2] - atompos[2], c+1) - 2*alpha*c*pow(pos[0] - atompos[0], a+1)*pow(pos[2] - atompos[2], c-1) )*pow(pos[1]-atompos[1], b)*expo;
        //             out(2,1) = out(1,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[1] - atompos[1], b+1)*pow(pos[2] - atompos[2], c+1) - 2*alpha*c*pow(pos[1] - atompos[1], b+1)*pow(pos[2] - atompos[2], c-1) )*pow(pos[0]-atompos[0], a)*expo;
        //         }
        //     }
        //     else if(b == 1){
        //         out(1,1) = getnorm()*(  4.0*alpha*alpha*pow(pos[1] - atompos[1], b+2) - 2*alpha*(2*b+1)*pow(pos[1] - atompos[1], b)  )*pow(pos[0] - atompos[0], a)*pow(pos[2]-atompos[2], c)*expo;
        //         out(1,0) = out(0,1) = getnorm()*(  4.0*alpha*alpha*pow(pos[0] - atompos[0], a+1)*pow(pos[1] - atompos[1], b+1) - 2*alpha*b*pow(pos[0] - atompos[0], a+1)*pow(pos[1] - atompos[1], b-1) )*pow(pos[2]-atompos[2], c)*expo;
        //         if(c == 0){
        //             out(2,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[2] - atompos[2], c+2) - 2*alpha*pow(pos[2] - atompos[2], c)  )*pow(pos[1] - atompos[1], b)*pow(pos[0]-atompos[0], a)*expo;
        //             out(2,0) = out(0,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[0] - atompos[0], a+1)*pow(pos[2] - atompos[2], c+1) )*pow(pos[1]-atompos[1], b)*expo;
        //             out(2,1) = out(1,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[1] - atompos[1], b+1)*pow(pos[2] - atompos[2], c+1) )*pow(pos[0]-atompos[0], a)*expo;
        //         }
        //         else if(c == 1){
        //             out(2,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[2] - atompos[2], c+2) - 2*alpha*(2*c+1)*pow(pos[2] - atompos[2], c)  )*pow(pos[1] - atompos[1], b)*pow(pos[0]-atompos[0], a)*expo;
        //             out(2,0) = out(0,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[0] - atompos[0], a+1)*pow(pos[2] - atompos[2], c+1) - 2*alpha*c*pow(pos[0] - atompos[0], a+1)*pow(pos[2] - atompos[2], c-1) )*pow(pos[1]-atompos[1], b)*expo;
        //             out(2,1) = out(1,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[1] - atompos[1], b+1)*pow(pos[2] - atompos[2], c+1) - 2*alpha*c*pow(pos[1] - atompos[1], b+1)*pow(pos[2] - atompos[2], c-1) )*pow(pos[0]-atompos[0], a)*expo;
        //         }
        //         else{
        //             out(2,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[2] - atompos[2], c+2) - 2*alpha*(2*c+1)*pow(pos[2] - atompos[2], c) + c*(c+1)*pow(pos[2] - atompos[2], c-2)  )*pow(pos[1] - atompos[1], b)*pow(pos[0]-atompos[0], a)*expo;
        //             out(2,0) = out(0,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[0] - atompos[0], a+1)*pow(pos[2] - atompos[2], c+1) - 2*alpha*c*pow(pos[0] - atompos[0], a+1)*pow(pos[2] - atompos[2], c-1) )*pow(pos[1]-atompos[1], b)*expo;
        //             out(2,1) = out(1,2) = getnorm()*(  4.0*alpha*alpha*pow(pos[1] - atompos[1], b+1)*pow(pos[2] - atompos[2], c+1) - 2*alpha*c*pow(pos[1] - atompos[1], b+1)*pow(pos[2] - atompos[2], c-1) )*pow(pos[0]-atompos[0], a)*expo;
        //         }
        //     }
        //     else{
        //         out(1,1) = getnorm()*(  4.0*alpha*alpha*pow(pos[1] - atompos[1], b+2) - 2*alpha*(2*b+1)*pow(pos[1] - atompos[1], b) + b*(b+1)*pow(pos[1] - atompos[1], b-2)  )*pow(pos[0] - atompos[0], a)*pow(pos[2]-atompos[2], c)*expo;
        //         out(1,0) = out(0,1) = getnorm()*(  4.0*alpha*alpha*pow(pos[0] - atompos[0], a+1)*pow(pos[1] - atompos[1], b+1) - 2*alpha*b*pow(pos[0] - atompos[0], a+1)*pow(pos[1] - atompos[1], b-1) )*pow(pos[2]-atompos[2], c)*expo;
        //     }
        // }
        // else if(a == 1){
        //     out(0,0) = getnorm()*(  4.0*alpha*alpha*pow(pos[0] - atompos[0], a+2) - 2*alpha*(2*a+1)*pow(pos[0] - atompos[0], a)  )*pow(pos[1] - atompos[1], b)*pow(pos[2]-atompos[2], c)*expo;
            
        //     if(b == 0){
        //         out(1,1) = getnorm()*(  4.0*alpha*alpha*pow(pos[1] - atompos[1], b+2) - 2*alpha*pow(pos[1] - atompos[1], b)  )*pow(pos[0] - atompos[0], a)*pow(pos[2]-atompos[2], c)*expo;
        //         out(1,0) = out(0,1) = getnorm()*(  4.0*alpha*alpha*pow(pos[0] - atompos[0], a+1)*pow(pos[1] - atompos[1], b+1) - 2*alpha*a*pow(pos[0] - atompos[0], a-1)*pow(pos[1] - atompos[1], b+1) )*pow(pos[2]-atompos[2], c)*expo;
        //         if(c == 0){
        //             out(2,2) = 0.0;
        //         }
        //         else if(c == 1){
        //             out(2,2) = 0.0;
        //         }
        //         else{
        //             out(2,2) = 0.0;
        //         }
        //     }
        //     else if(b == 1){
        //         out(1,1) = getnorm()*(  4.0*alpha*alpha*pow(pos[1] - atompos[1], b+2) - 2*alpha*(2*b+1)*pow(pos[1] - atompos[1], b)  )*pow(pos[0] - atompos[0], a)*pow(pos[2]-atompos[2], c)*expo;
        //         out(1,0) = out(0,1) = getnorm()*(  4.0*alpha*alpha*pow(pos[0] - atompos[0], a+1)*pow(pos[1] - atompos[1], b+1) - 2*alpha*b*pow(pos[0] - atompos[0], a+1)*pow(pos[1] - atompos[1], b-1) )*pow(pos[2]-atompos[2], c)*expo;
        //     }
        //     else{
        //         //
        //     }
        // }
        // else{
        //     out(0,0) = getnorm()*(  4.0*alpha*alpha*pow(pos[0] - atompos[0], a+2) - 2*alpha*(2*a+1)*pow(pos[0] - atompos[0], a) + a*(a+1)*pow(pos[0] - atompos[0], a-2)  )*pow(pos[1] - atompos[1], b)*pow(pos[2]-atompos[2], c)*expo;
        // }
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
    //get value of the derivative of the wavefunction w.r.t. x, y and z
    vec3 getderiv(const vec3& point){
        vec3 out;
        out << 0.0,0.0,0.0;
        for(int i=0; i<get_size(); i++){
            out += coeff[i]*GTOlist[i].getderiv(point);
        }
        return out;
    }
    //get Hessian of the wavefunction:
    Eigen::MatrixXd getHessian(const vec3& point){
        Eigen::MatrixXd out = Eigen::MatrixXd::Zero(3,3);
        for(int i=0; i<get_size(); i++){
            out += coeff[i]*GTOlist[i].getHessian(point);
        }
        return out;
    }
};

#endif //_ORBITALS_H