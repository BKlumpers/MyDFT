/*****************************************************************

Basic closed-shell spin-restricted HF-solver for simple molecules using STO-NG

Author: B. Klumpers (bartkl@live.nl)

Published under GNU General Public License 3.0

Allows for SCF-computation of molecular energies for simple molecules.
Testcases for H, He, H2, HeH+ and He2 are included.

*****************************************************************/

#include "prelim_math.h"

//****************************************************************
//preliminary mathematics:


//compute the factorial of an integer n!
int factorial(int n)
{
    if(n < 0) {
        if(EFLAG){
            cout << endl << "Error: Non-valid angular momentum supplied to factorial." << endl << endl;
        }
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
            if(EFLAG){
            cout << endl << "Error: Non-valid angular momentum supplied to double factorial." << endl << endl;
            }
        }
        return 1;
    }
    else if(k == 1){
        return 1;
    }
    else{
        return factorial(2*k)/(pow(2,k)*factorial(k));
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

//End of file
//****************************************************************
