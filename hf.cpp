/*****************************************************************

Basic HF-solver for single atoms using STO-NG

Author: B. Klumpers (bartkl@live.nl)

Published under GNU General Public License 3.0

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
double nuclearGTO(GTO GTO1, GTO GTO2, double charge, vec3 nucleus);
double nuclearComponent(int l, int r, int i, int a1, int a2, double PA, double PB, double PC, double alpha_c);
double fmch(int nu, double x);
double two_electronGTO(GTO GTO1, GTO GTO2, GTO GTO3, GTO GTO4);
double electricComponent(int l1, int l2, int r1, int r2, int i, int a1, int a2, int a3, int a4, double AX, double BX, double CX, double DX, double PX, double QX, double g12, double g34, double delta);
double electricFourier(int l, int a1, int a2, double PA, double PB, double r1, double gamma);
double overlapCGF(CGF CGF1, CGF CGF2);
double kineticCGF(CGF CGF1, CGF CGF2);
double nuclearCGF(CGF CGF1, CGF CGF2, double charge, vec3 pos);
double two_electronCGF(CGF CGF1, CGF CGF2, CGF CGF3, CGF CGF4);

//****************************************************************
//begin main codeblock

int main() {
    std::cout << "Initialising" << std::endl;

    //input for single hydrogen atom with STO-3G
    double Z = 1; //core charge for Hydrogen
    vec3 pos;
    pos << 0,0,0; //position nucleus at origin

    //create contracted gaussian function STO-3G
    CGF cgfH;
    cgfH.add_gto(0.154,3.425,0,0,0,pos);
    cgfH.add_gto(0.535,0.624,0,0,0,pos);
    cgfH.add_gto(0.445,0.169,0,0,0,pos);
    //Example: compute overlap, kinetic and nuclear integrals
    double overlap = overlapCGF(cgfH,cgfH);
    double kinetic = kineticCGF(cgfH,cgfH);
    double nuclear = nuclearCGF(cgfH,cgfH,1,pos);
    cout << "overlap: " << overlap << endl;
    cout << "kinetic: " << kinetic << endl;
    cout << "nuclear: " << nuclear << endl;

    /*reserved for energy minimisation
    *
    *
    *********************************/

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
    double prefactor = exp(-GTO1.alpha*GTO2.alpha*((GTO1.atompos-GTO2.atompos).squaredNorm())/(GTO1.alpha+GTO2.alpha))*pow(pi/(GTO1.alpha+GTO2.alpha),1.5);
    vec3 rc = (GTO1.alpha*GTO1.atompos + GTO2.alpha*GTO2.atompos)/(GTO1.alpha+GTO2.alpha);
    double Sx = 0;
    for(int w=0; w<(1+0.5*(GTO1.a+GTO2.a)); w++){
        Sx = Sx + binomialComposite(w, GTO1.a, GTO2.a, rc[0]-GTO1.atompos[0], rc[0]-GTO2.atompos[0])*doublefactorial_odd(w)/(pow(2*(GTO1.alpha+GTO2.alpha),w));
    }
    double Sy = 0;
    for(int w=0; w<(1+0.5*(GTO1.b+GTO2.b)); w++){
        Sy = Sy + binomialComposite(w, GTO1.b, GTO2.b, rc[1]-GTO1.atompos[1], rc[1]-GTO2.atompos[1])*doublefactorial_odd(w)/(pow(2*(GTO1.alpha+GTO2.alpha),w));
    }
    double Sz = 0;
    for(int w=0; w<(1+0.5*(GTO1.c+GTO2.c)); w++){
        Sz = Sz + binomialComposite(w, GTO1.c, GTO2.c, rc[2]-GTO1.atompos[2], rc[2]-GTO2.atompos[2])*doublefactorial_odd(w)/(pow(2*(GTO1.alpha+GTO2.alpha),w));
    }
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
    kinetic = kinetic - 2*GTO2.alpha*GTO2.alpha*overlapGTO(GTO1,GTOki)*GTO2.getnorm()/GTOki.getnorm(); //terms 2 orders higher in a
    GTOki.a = GTO2.a;
    GTOki.b = GTO2.b + 2;
    kinetic = kinetic - 2*GTO2.alpha*GTO2.alpha*overlapGTO(GTO1,GTOki)*GTO2.getnorm()/GTOki.getnorm(); //terms 2 orders higher in b
    GTOki.b = GTO2.b;
    GTOki.c = GTO2.c + 2;
    kinetic = kinetic - 2*GTO2.alpha*GTO2.alpha*overlapGTO(GTO1,GTOki)*GTO2.getnorm()/GTOki.getnorm(); //terms 2 orders higher in c
    GTOki.c = GTO2.c;

    GTOki.a = GTO2.a - 2;
    kinetic = kinetic - 0.5*GTO2.a*(GTO2.a-1)*overlapGTO(GTO1,GTOki)*GTO2.getnorm()/GTOki.getnorm(); //terms 2 orders lower in a
    GTOki.a = GTO2.a;
    GTOki.b = GTO2.b - 2;
    kinetic = kinetic - 0.5*GTO2.b*(GTO2.b-1)*overlapGTO(GTO1,GTOki)*GTO2.getnorm()/GTOki.getnorm(); //terms 2 orders lower in b
    GTOki.b = GTO2.b;
    GTOki.c = GTO2.c - 2;
    kinetic = kinetic - 0.5*GTO2.c*(GTO2.c-1)*overlapGTO(GTO1,GTOki)*GTO2.getnorm()/GTOki.getnorm(); //terms 2 orders lower in c

    return kinetic; //factors for getnorm added to correct redundant normalisation
}

//compute nuclear integral of GTO
//compute nuclear integral of <gto1|(r-nucleus)^-1|gto2>
double nuclearGTO(GTO GTO1, GTO GTO2, double charge, vec3 nucleus)
{
    static const double pi = 3.14159265359;
    double nuclear = 0;
    double prefactor = (2*pi/(GTO1.alpha + GTO2.alpha))*exp(-GTO1.alpha*GTO2.alpha*((GTO1.atompos-GTO2.atompos).squaredNorm())/(GTO1.alpha+GTO2.alpha));
    vec3 rc = (GTO1.alpha*GTO1.atompos + GTO2.alpha*GTO2.atompos)/(GTO1.alpha+GTO2.alpha);
    //x,y,z-components of the integral
    double Ax;
    double Ay;
    double Az;
    //9-fold summation to yield V(nuclear)
    for(int l=0; l<(1 + GTO1.a + GTO2.a); l++){
        for(int m=0; m<(1 + GTO1.b + GTO2.b); m++){
            for(int n=0; n<(1 + GTO1.c + GTO2.c); n++){
                for(int r=0; r<(1+0.5*l); r++){
                    for(int s=0; s<(1+0.5*m); s++){
                        for(int t=0; t<(1+0.5*n); t++){
                            for(int i=0; i<(1+0.5*l-r); i++){
                                Ax = nuclearComponent(l,r,i,GTO1.a,GTO2.a,rc[0]-GTO1.atompos[0],rc[0]-GTO2.atompos[0],rc[0]-nucleus[0],GTO1.alpha+GTO2.alpha);
                                for(int j=0; j<(1+0.5*m-s); j++){
                                    Ay = nuclearComponent(m,s,j,GTO1.b,GTO2.b,rc[1]-GTO1.atompos[1],rc[1]-GTO2.atompos[1],rc[1]-nucleus[1],GTO1.alpha+GTO2.alpha);
                                    for(int k=0; k<(1+0.5*n-t); k++){
                                        Az = nuclearComponent(n,t,k,GTO1.c,GTO2.c,rc[2]-GTO1.atompos[2],rc[2]-GTO2.atompos[2],rc[2]-nucleus[2],GTO1.alpha+GTO2.alpha);
                                        nuclear = nuclear + Ax*Ay*Az*fmch(l+m+n-2*r-2*s-2*t-i-j-k,((rc-nucleus).squaredNorm())*(GTO1.alpha+GTO2.alpha));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return -charge*prefactor*nuclear*GTO1.getnorm()*GTO2.getnorm();
}

//function for x,y,z-dependent components of the nuclear integral
double nuclearComponent(int l, int r, int i, int a1, int a2, double PA, double PB, double PC, double alpha_c)
{
    return pow(-1,l+i)*binomialComposite(l,a1,a2,PA,PB)*factorial(l)*pow(PC,l-2*r-2*i)/(factorial(i)*factorial(r)*factorial(l-2*r-2*i)*pow(4*alpha_c,r+i));
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
            if(add>fmchTol){
                sum = sum + add;
            }
            else{
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
double two_electronGTO(GTO GTO1, GTO GTO2, GTO GTO3, GTO GTO4)
{
    static const double pi = 3.14159265359;
    double electric = 0;
    double gamma12 = GTO1.alpha + GTO2.alpha;
    double gamma34 = GTO3.alpha + GTO4.alpha;
    double delta = (0.25/gamma12) + (0.25/gamma34);
    double Rab_2 = (GTO1.atompos-GTO2.atompos).squaredNorm();
    double Rcd_2 = (GTO3.atompos-GTO4.atompos).squaredNorm();
    double prefactor = 2*pi*pi*sqrt(pi)/(gamma12*gamma34*sqrt(gamma12+gamma34))*exp(-GTO1.alpha*GTO2.alpha*Rab_2/gamma12)*exp(-GTO3.alpha*GTO4.alpha*Rcd_2/gamma34);
    vec3 rp = (GTO1.alpha*GTO1.atompos + GTO2.alpha*GTO2.atompos)/gamma12;
    vec3 rq = (GTO3.alpha*GTO3.atompos + GTO4.alpha*GTO4.atompos)/gamma34;
    double Rpq_2 = (rp-rq).squaredNorm();
    //x,y,z-components of the 2e-integral:
    double Bx;
    double By;
    double Bz;
    //15-fold summation for the 2e-integral:
    for(int l=0;l<(1+GTO1.a+GTO2.a); l++){
        for(int r=0; r<(1+0.5*l); r++){
            for(int i=0; i<(0.5*l-r); i++){
                for(int l2=0; l2<(1+GTO3.a+GTO4.a); l2++){
                    for(int r2=0; r2<(1+0.5*l2); r2++){
                        Bx=electricComponent(l,l2,r,r2,i,GTO1.a,GTO2.a,GTO3.a,GTO4.a,GTO1.atompos[0],GTO2.atompos[0],GTO3.atompos[0],GTO4.atompos[0],rp[0],rq[0],gamma12,gamma34,delta);
                        for(int m=0; m<(1+GTO1.b+GTO2.b); m++){
                            for(int s=0; s<(1+0.5*m); s++){
                                for(int j=0; j<(0.5*m-s); j++){
                                    for(int m2=0; m2<(1+GTO3.b+GTO4.b); m2++){
                                        for(int s2=0; s2<(1+0.5*m2); s2++){
                                            By=electricComponent(m,m2,s,s2,j,GTO1.b,GTO2.b,GTO3.b,GTO4.b,GTO1.atompos[1],GTO2.atompos[1],GTO3.atompos[1],GTO4.atompos[1],rp[1],rq[1],gamma12,gamma34,delta);
                                            for(int n=0; n<(1+GTO1.c+GTO2.c); n++){
                                                for(int t=0; t<(1+0.5*n); t++){
                                                    for(int k=0; k<(0.5*n-t); k++){
                                                        for(int n2=0; n2<(1+GTO3.c+GTO4.c); n2++){
                                                            for(int t2=0; t2<(1+0.5*n2); t2++){
                                                                Bz=electricComponent(n,n2,t,t2,k,GTO1.c,GTO2.c,GTO3.c,GTO4.c,GTO1.atompos[2],GTO2.atompos[2],GTO3.atompos[2],GTO4.atompos[2],rp[2],rq[2],gamma12,gamma34,delta);
                                                                electric = electric + Bx*By*Bz*fmch(l+l2+m+m2+n+n2-2*(r+r2+s+s2+t+t2)-i-j-k,Rpq_2/(4*delta));
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return electric*prefactor*GTO1.getnorm()*GTO2.getnorm();
}

//function for x,y,z-dependent components of the 2e-integral
double electricComponent(int l1, int l2, int r1, int r2, int i, int a1, int a2, int a3, int a4, double AX, double BX, double CX, double DX, double PX, double QX, double g12, double g34, double delta)
{
    return pow(-1,l2+i)*electricFourier(l1,a1,a2,PX-AX,PX-BX,r1,g12)*electricFourier(l2,a3,a4,PX-CX,PX-DX,r2,g34)*pow(2*delta,2*(r1+r2))*factorial(l1+l2-2*r1-2*r2)*pow(delta,i)*pow(PX-QX,l1+l2-2*(r1+r2+i))/(pow(4*delta,l1+l2)*factorial(i)*factorial(l1+l2-2*(r1+r2+i)));
}

//function for Fourier-components of electricComponent
double electricFourier(int l, int a1, int a2, double PA, double PB, double r1, double gamma)
{
    return binomialComposite(l,a1,a2,PA,PB)*factorial(l)*pow(gamma,r1-l)/(factorial(r1)*factorial(l-2*r1));
}



//CGF{STO-3G} expansion:

//overlap integral for STO-NG CGF
double overlapCGF(CGF CGF1, CGF CGF2)
{
    double overlap = 0;
    for(int w=0; w<CGF1.get_size(); w++){
        for(int v=0; v<CGF2.get_size(); v++){
            overlap = overlap + CGF1.coeff[w]*CGF2.coeff[v]*overlapGTO(CGF1.GTOlist[w],CGF2.GTOlist[v]);
        }
    }
    return overlap;
}

//compute kinetic integral of CGF
double kineticCGF(CGF CGF1, CGF CGF2)
{
    double kinetic = 0;
    for(int w=0; w<CGF1.get_size(); w++){
        for(int v=0; v<CGF2.get_size(); v++){
            kinetic = kinetic + CGF1.coeff[w]*CGF2.coeff[v]*kineticGTO(CGF1.GTOlist[w],CGF2.GTOlist[v]);
        }
    }
    return kinetic;
}

//compute nuclear integral of CGF
double nuclearCGF(CGF CGF1, CGF CGF2, double charge, vec3 pos)
{
    double nuclear = 0;
    for(int w=0; w<CGF1.get_size(); w++){
        for(int v=0; v<CGF2.get_size(); v++){
            nuclear = nuclear + CGF1.coeff[w]*CGF2.coeff[v]*nuclearGTO(CGF1.GTOlist[w],CGF2.GTOlist[v], charge, pos);
        }
    }
    return nuclear;
}

//compute 2e-integral of CGF
double two_electronCGF(CGF CGF1, CGF CGF2, CGF CGF3, CGF CGF4)
{
    double electric = 0;
    for(int w=0; w<CGF1.get_size(); w++){
        for(int v=0; v<CGF2.get_size(); v++){
            for(int w2=0; w2<CGF3.get_size(); w2++){
                for(int v2=0; w2<CGF4.get_size(); v2++){
                    electric = electric + CGF1.coeff[w]*CGF2.coeff[v]*CGF3.coeff[w2]*CGF4.coeff[v2]*two_electronGTO(CGF1.GTOlist[w],CGF1.GTOlist[v],CGF3.GTOlist[w2],CGF4.GTOlist[v2]);
                }
            }
        }
    }
    return electric;
}



//end of file

