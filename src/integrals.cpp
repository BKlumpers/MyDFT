/*****************************************************************

Basic closed-shell spin-restricted HF-solver for simple molecules using STO-NG

Author: B. Klumpers (bartkl@live.nl)

Published under GNU General Public License 3.0

Allows for SCF-computation of molecular energies for simple molecules.
Testcases for H, He, H2, HeH+ and He2 are included.

*****************************************************************/

#include "integrals.h"

//****************************************************************
//Evaluation of integrals for Gaussian Type Orbitals:


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
    for(int i=max(0,2*k-a2); i<(1+min(2*k,a1)); i++){
        out = out + binomialCoeff(a1, i)*binomialCoeff(a2, 2*k-i)*pow(r1,a1-i)*pow(r2,a2+i-2*k);
    }
    return out;
}

//compute kinetic integral of GTO < GTO1 | -Laplacian_Delta/2 | GTO2 >
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

//compute nuclear integral of GTO < GTO1 | (r-nucleus)^-1 | GTO2 >
double nuclearGTO(GTO GTO1, GTO GTO2, double charge, const vec3& nucleus)
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
//
//Algorithm adapted from:   David B. Cook, Handbook of Computational Quantum Chemistry
//                          Dover (2005), ISBN: 978-0-486-44307-2
//
double fmch(int nu, double x)
{
    double sum = 0;
    double m = double(nu); //cast nu to double for evaluation
    double add;
    static const double fmchTol = 1e-8; //convergence tolerance for the series
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
            for(int i=0; i<(1+0.5*l-r); i++){
                for(int l2=0; l2<(1+GTO3.a+GTO4.a); l2++){
                    for(int r2=0; r2<(1+0.5*l2); r2++){
                        Bx=electricComponent(l,l2,r,r2,i,GTO1.a,GTO2.a,GTO3.a,GTO4.a,GTO1.atompos[0],GTO2.atompos[0],GTO3.atompos[0],GTO4.atompos[0],rp[0],rq[0],gamma12,gamma34,delta);
                        for(int m=0; m<(1+GTO1.b+GTO2.b); m++){
                            for(int s=0; s<(1+0.5*m); s++){
                                for(int j=0; j<(1+0.5*m-s); j++){
                                    for(int m2=0; m2<(1+GTO3.b+GTO4.b); m2++){
                                        for(int s2=0; s2<(1+0.5*m2); s2++){
                                            By=electricComponent(m,m2,s,s2,j,GTO1.b,GTO2.b,GTO3.b,GTO4.b,GTO1.atompos[1],GTO2.atompos[1],GTO3.atompos[1],GTO4.atompos[1],rp[1],rq[1],gamma12,gamma34,delta);
                                            for(int n=0; n<(1+GTO1.c+GTO2.c); n++){
                                                for(int t=0; t<(1+0.5*n); t++){
                                                    for(int k=0; k<(1+0.5*n-t); k++){
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
    return electric*prefactor*GTO1.getnorm()*GTO2.getnorm()*GTO3.getnorm()*GTO4.getnorm();
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


//****************************************************************
//Extension of the above GTO-integrals to Contracted Gaussian Functions:


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
double nuclearCGF(CGF CGF1, CGF CGF2, double charge, const vec3& pos)
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
                for(int v2=0; v2<CGF4.get_size(); v2++){
                    electric = electric + CGF1.coeff[w]*CGF2.coeff[v]*CGF3.coeff[w2]*CGF4.coeff[v2]*two_electronGTO(CGF1.GTOlist[w],CGF2.GTOlist[v],CGF3.GTOlist[w2],CGF4.GTOlist[v2]);
                }
            }
        }
    }
    return electric;
}

//sort indices of the 2e-integral using symmetry
int two_electronSort(int i, int j, int k, int l)
{
    if(i < j) {
        swap(i,j); //invariance <ij|kl> = <ji|kl>
    }
    if(k < l) {
        swap(k,l); //invariance <ij|kl> = <ij|lk>
    }

    int ij = i * (i + 1) / 2 + j; //map i,j->(ij) compound index
    int kl = k * (k + 1) / 2 + l; //map k,l->(kl) compound index

    if(ij < kl) {
        swap(ij,kl); //invariance <ij|kl> = <kl|ij>
    }

    return ij * (ij + 1) / 2 + kl; //map ij,kl->(ijkl) compound index
}

//End of file
//****************************************************************
