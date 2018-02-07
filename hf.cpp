/*****************************************************************

Basic HF-solver for single atoms using STO-NG

Author: B. Klumpers (bartkl@live.nl)

Published under GNU General Public License 3.0

dev-notes:
*no error handling implemented, only printing error messages
*TODO: do not calculate norm everytime for get_norm

//norm was found to be correct in case of zero angular momentum (compared with analytical solutions)
//factorial, doublefactorial_odd and binomialCoeff confirmed to yield expected results
//binomial_composite agrees with hfcxx
//overlap integral agrees with hfcxx and calculations by hand
//fmch confirmed to yield same results as numerical integration procedure

//overlapCGF and kineticGTO unconfirmed if working (WIP)
//nuclearGTO not yet implemented (empty function call)

*****************************************************************/

#include <iostream>
#include <cmath>
#include <vector>
#include <boost/math/special_functions/gamma.hpp>
#include <Eigen/Core>

typedef Eigen::Vector3d vec3; //define vec3 as the Eigen::Vector3d-object

using namespace std;

//define function calls (part1)
int factorial(int n);
int doublefactorial_odd(int k);

//Gaussian type orbital
class GTO {
public:
    double alpha;           //exponent factor
    int a,b,c;              //angular momentum components in x,y,z respectively
    vec3 atompos;           //atomic coordinates of the nucleus
    double getnorm()
    {
        static const double pi = 3.14159265359;
        //compute normalisation constant
        return sqrt(pow(2*alpha/pi, 1.5)*pow(4*alpha, a + b + c)/(doublefactorial_odd(a)*doublefactorial_odd(b)*doublefactorial_odd(c)));
    }
};

class CGF {
public:
    vector<double> coeff; //fix variable number of GTOs per CGF
    vector<GTO> GTOlist; //list of GTOs
    //no direct normalisation, assumed CGF of normalised GTO (hence CGF not inherently normalised)
    void add_gto(double gto_coeff, double gto_alpha, int gto_a, int gto_b, int gto_c, vec3 gto_atompos)
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
};

//define function calls (part2)
int binomialCoeff(int n, int k);
double overlapGTO(GTO GTO1, GTO GTO2);
double binomialComposite(int k, int a1, int a2, double r1, double r2);
double kineticGTO(GTO GTO1, GTO GTO2);
double nuclearGTO(GTO GTO1, GTO GTO2, vec3 nucleus);
double fmch(int nu, double x);
double overlapCGF(CGF CGF1, CGF CGF2);

//****************************************************************
//begin main codeblock

int main() {
    std::cout << "Initialising" << std::endl;
    double N = 3; //number of Gaussians per Slater-type orbital
    double Z = 1; //core charge

    //testblock
    GTO testGTO;
    testGTO.alpha = 1;
    testGTO.a = 0;
    testGTO.b = 0;
    testGTO.c = 0;
    testGTO.atompos << 0,0,0;
    double norm = testGTO.getnorm();

    GTO testGTO2;
    testGTO2.alpha = 2;
    testGTO2.a = 0;
    testGTO2.b = 0;
    testGTO2.c = 0;
    testGTO2.atompos << 0,0,0;
    double norm2 = testGTO2.getnorm();

    CGF cgftest;
    vec3 testpos;
    testpos << 0,0,0;
    cgftest.add_gto(0.5,1,2,3,4,testpos);
    cgftest.add_gto(0.8,5,6,7,8,testpos);
    double tempy;
    tempy = cgftest.GTOlist[1].alpha;
    cout << "cgf: " << tempy << endl;

    // for(int i = 0; i<3 ; i++){
    //     cout << "atompos: [" << i << "] " << testGTO.atompos[i] << endl;
    // }
    cout << "norm: " << norm << " also2: " << norm2 << endl;
    double overlapTest = overlapGTO(testGTO,testGTO2);
    cout << "overlap: " << overlapTest << endl;
    cout << "boosty: " << boost::math::tgamma(0.5) << " should be equal to: " << sqrt(3.14159265359) << endl; //test proper call to boost
    double Fnu_test = fmch(0,1);
    cout << "Fnu_test: " << Fnu_test << " should be: 0.74682413 (numerical)" << endl;
    double Fnu_test2 = fmch(1,1);
    cout << "Fnu_test: " << Fnu_test2 << " should be: 0.18947235 (numerical)" << endl;
    // vec3 vectest(1,1,0);
    // for(int i = 0; i<3 ; i++){
    //     cout << "vectest: [" << i << "] " << vectest[i] << endl;
    // }

    //testblock end

    std::cout << "Program has ended" << std::endl;
    return 0;
}

//end main codeblock
//****************************************************************

//compute the factorial of an integer n!
int factorial(int n)
{
    if(n < 0) {
        cout << endl << "Error: Non-valid angular momentum supplied to factorial." << endl << endl;
        return 1;
    }
    if(n > 1){
        int out = 1;
        for(int count = n; count > 0; count--){
            out = out*count;
            //cout << n << " gives " << out << endl; //debug output
        }
        return out;
    }
    else{
        return 1; //0!=1
    }
}

//compute double factorial (2k-1)!!
int doublefactorial_odd(int k)
{
    if(k < 1){
        if(k!=0){
            cout << endl << "Error: Non-valid angular momentum supplied to double factorial." << endl << endl;
        }
        return 1;
    }
    else if(k == 1){
        //cout << "feedback" << endl; //debug output
        return 1;
    }
    else{
        return factorial(2*k)/(std::pow(2,k)*factorial(k));
    }
}

//compute binomial coefficient (n  k)
int binomialCoeff(int n, int k)
{
    if((n > k) && (k >= 0)){
        return factorial(n)/(factorial(k)*factorial(n-k));
    }
    else{
        if((n > k) && (n != k)){
            cout << "Error: Negative values supplied to binomialCoeff." << endl;
        }
        else{
            if((n!=0) && (k!=0) && (n != k)){
                cout << "Error: Subset exceeds set (k > n) in binomialCoeff" << endl;
            }
        }
        return 1;
    }
}

//compute overlap integral of 2 gaussian type orbitals < GTO1 | GTO2 >
double overlapGTO(GTO GTO1, GTO GTO2)
{
    static const double pi = 3.14159265359;
//    double overlap = 0.0;
    double prefactor = exp(-GTO1.alpha*GTO2.alpha*((GTO1.atompos-GTO2.atompos).squaredNorm())/(GTO1.alpha+GTO2.alpha))*pow(pi/(GTO1.alpha+GTO2.alpha),1.5);
    vec3 rc = (GTO1.alpha*GTO1.atompos + GTO2.alpha*GTO2.atompos)/(GTO1.alpha+GTO2.alpha);
    double Sx = 0;
    for(int w=0; w<(1+0.5*(GTO1.a+GTO2.a)); w++){ //check limits
        Sx = Sx + binomialComposite(w, GTO1.a, GTO2.a, rc[0]-GTO1.atompos[0], rc[0]-GTO2.atompos[0])*doublefactorial_odd(w)/(pow(2*(GTO1.alpha+GTO2.alpha),w));
    }
    //cout << "Sx: " << Sx << " pre: " << prefactor << endl;
    double Sy = 0;
    for(int w=0; w<(1+0.5*(GTO1.b+GTO2.b)); w++){ //check limits
        Sy = Sy + binomialComposite(w, GTO1.b, GTO2.b, rc[1]-GTO1.atompos[1], rc[1]-GTO2.atompos[1])*doublefactorial_odd(w)/(pow(2*(GTO1.alpha+GTO2.alpha),w));
    }
    double Sz = 0;
    for(int w=0; w<(1+0.5*(GTO1.c+GTO2.c)); w++){ //check limits
        Sz = Sz + binomialComposite(w, GTO1.c, GTO2.c, rc[2]-GTO1.atompos[2], rc[2]-GTO2.atompos[2])*doublefactorial_odd(w)/(pow(2*(GTO1.alpha+GTO2.alpha),w));
    }
    //double overlap = GTO1.getnorm()*GTO2.getnorm()*prefactor*Sx*Sy*Sz;
    //cout << "Sy: " << Sy << " Sz: " << Sz << endl;
    //cout << "norm1: " << GTO1.getnorm() << " norm2: " << GTO2.getnorm() << endl;
    return GTO1.getnorm()*GTO2.getnorm()*prefactor*Sx*Sy*Sz;
}

//compute the multi-binomial coefficient:
double binomialComposite(int k, int a1, int a2, double r1, double r2)
{
    double out = 0;
    for(int i=max(0,2*k-a2); i<(1+min(2*k,a1)); i++){ //check limits
        out = out + binomialCoeff(a1, i)*binomialCoeff(a2, 2*k-i)*pow(r1,a1-i)*pow(r2,a2+i-2*k);
    }
    return out;
}

//compute kinetic integral of GTO
double kineticGTO(GTO GTO1, GTO GTO2)
{
    GTO GTOki = GTO2; //local manipulatable GTO
    double kinetic = GTO2.alpha*(2*GTO2.a+2*GTO2.b+2*GTO2.c+3)*overlapGTO(GTO1, GTO2); //terms of same order

    GTOki.a = GTO2.a + 2;
    kinetic = kinetic - 2*GTO2.alpha*GTO2.alpha*overlapGTO(GTO1,GTOki); //terms 2 orders higher in a
    GTOki.a = GTO2.a;
    GTOki.b = GTO2.b + 2;
    kinetic = kinetic - 2*GTO2.alpha*GTO2.alpha*overlapGTO(GTO1,GTOki); //terms 2 orders higher in b
    GTOki.b = GTO2.b;
    GTOki.c = GTO2.c + 2;
    kinetic = kinetic - 2*GTO2.alpha*GTO2.alpha*overlapGTO(GTO1,GTOki); //terms 2 orders higher in c
    GTOki.c = GTO2.c;

    GTOki.a = GTO2.a - 2;
    kinetic = kinetic - 0.5*GTO2.a*(GTO2.a-1)*overlapGTO(GTO1,GTOki); //terms 2 orders lower in a
    GTOki.a = GTO2.a;
    GTOki.b = GTO2.b - 2;
    kinetic = kinetic - 0.5*GTO2.b*(GTO2.b-1)*overlapGTO(GTO1,GTOki); //terms 2 orders lower in b
    GTOki.b = GTO2.b;
    GTOki.c = GTO2.c - 2;
    kinetic = kinetic - 0.5*GTO2.c*(GTO2.c-1)*overlapGTO(GTO1,GTOki); //terms 2 orders lower in c

    return kinetic;
}



//compute nuclear integral of GTO
double nuclearGTO(GTO GTO1, GTO GTO2, vec3 nucleus)
{
    double nuclear;
    nuclear = 1; //temp for compilation
    //compute nuclear integral of <gto1|(r-nucleus)^-1|gto2>
    return nuclear;
}

//compute the incomplete gaussian integral used in evaluation of the nuclear and 2e-integrals
//F_nu(x) = integral from 0 to 1 of: t^2nu * exp(-x*t^2)dt with respect to dt
double fmch(int nu, double x)
{
    double sum = 0;
    double m = double(nu); //cast nu to double for evaluation
    double add;
    double fmchTol = 1e-8; //convergence tolerance for the series
    if(x <= 10){
        //lower expansion
        for(int i=0; i<51; i++){
            add = pow(x,i)*boost::math::tgamma(m+0.5)/boost::math::tgamma(m+i+1.5);
            //cout << add << endl;
            if(add>fmchTol){
                sum = sum + add;
                //cout << "add" << endl;
            }
            else{
                //cout << "break" << endl;
                break;
            }
            if(i==50){
                cout << "Warning: fmch(lower) most likely not converged" << endl;
            }
        }
        return 0.5*exp(-x)*sum;
    }
    else{
        //upper expansion
        for(int i=0; i<51; i++){
            add = pow(x,-i)*boost::math::tgamma(m+0.5)/boost::math::tgamma(m-i+1.5);
            if(add>fmchTol){
                sum = sum + add;
            }
            else{
                break;
            }
            if(i==50){
                cout << "Warning: fmch(upper) most likely not converged" << endl;
            }
        }
        double limit = boost::math::tgamma(m+0.5)/(2*pow(x,m+0.5));
        return limit-0.5*exp(-x)*sum;
    }
}

//compute 2e-integral of GTO



//CGF{STO-3G} expansion:

//overlap integral for STO-NG CGF
double overlapCGF(CGF CGF1, CGF CGF2)
{
    double overlap = 0;
    for(int w=0; w<(1+CGF1.get_size()); w++){
        for(int v=0; v<(1+CGF2.get_size()); v++){
            overlap = overlap + CGF1.coeff[w]*CGF2.coeff[v]*overlapGTO(CGF1.GTOlist[w],CGF1.GTOlist[v]);
        }
    }
    return overlap;
}


//compute kinetic integral of CGF

//compute nuclear integral of CGF

//compute 2e-integral of CGF




//end of file

