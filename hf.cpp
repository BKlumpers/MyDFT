/*****************************************************************

Basic closed-shell spin-restricted HF-solver for simple molecules using STO-NG

Author: B. Klumpers (bartkl@live.nl)

Published under GNU General Public License 3.0

Allows for SCF-computation of molecular energies for simple molecules.
Testcases for H, He and H2 are included.

*****************************************************************/

#include <iostream>
#include <cmath>
#include <vector>
#include <boost/math/special_functions/gamma.hpp>
#include <Eigen/Dense>

typedef Eigen::Vector3d vec3; //define vec3 as the Eigen::Vector3d-object

using namespace std;

const bool EFLAG = false; //turn mathematical errors off

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

//Contracted Gaussian orbital
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

//parse output of SCF energy minimisation
class SCF_E {
public:
    double energy; //total energy of the system
    Eigen::VectorXd orbital_energies; //energies of the molecular orbitals
    Eigen::MatrixXd C_vector; //coefficient matrix
    void SCF_result(double Result_energy, Eigen::VectorXd Result_orbital_energies, Eigen::MatrixXd Result_C_vector)
    {
        energy = Result_energy;
        orbital_energies = Result_orbital_energies;
        C_vector = Result_C_vector;
    }
};

//define function calls (part2)
SCF_E SCF_HF_energy(vector<CGF> AO_list, vector<vec3> pos_list, vector<double> charge_list, vector<int> nelec_list);
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
int two_electronSort(int i, int j, int k, int l);

//****************************************************************
//begin main codeblock

int main() {
    cout << "Initialising" << endl << endl;

    //input parameters for single atoms and simple diatomics:
    double H_Z = 1; //core charge for Hydrogen
    double He_Z = 2; //core charge for Helium
    double C_Z = 6; //core charge for Carbon
    double O_Z = 8; //core charge for Oxygen
    int H_nelec = 1; //1 electron on Hydrogen
    int He_nelec = 2; //2 electrons on Helium
    int C_nelec = 6; //6 electrons on Carbon
    int O_nelec = 8; //8 electons on Oxygen
    vec3 pos;
    pos << 0,0,0; //position nucleus at origin
    vec3 pos2;
    pos2 << 0,0,1.4; //position second atom 1.4A away from the 1st atom (H-H bond length)
    vec3 pos3;
    pos3 << 0,0,2.116; //position third atom 1.4A away from the 1st atom (C-O bond length in C=O)

    //create contracted gaussian function STO-3G for Hydrogen
    CGF cgfH;
    cgfH.add_gto(0.154,3.425,0,0,0,pos);
    cgfH.add_gto(0.535,0.624,0,0,0,pos);
    cgfH.add_gto(0.445,0.169,0,0,0,pos);
    //Compute integrals for Hydrogen (only 1 electron, so no repulsion)
    double H_overlap = overlapCGF(cgfH,cgfH);
    double H_kinetic = kineticCGF(cgfH,cgfH);
    double H_nuclear = nuclearCGF(cgfH,cgfH,H_Z,pos);
    cout << "Testcase for Hydrogen: " << endl;
    cout << "H_overlap: " << H_overlap << endl;
    cout << "H_kinetic: " << H_kinetic << endl;
    cout << "H_nuclear: " << H_nuclear << endl;
    cout << "H_Energy: " << H_kinetic + H_nuclear << endl << endl;;

    //create CGF STO-3G for Helium
    CGF cgfHe;
    cgfHe.add_gto(0.154,6.362,0,0,0,pos);
    cgfHe.add_gto(0.535,1.159,0,0,0,pos);
    cgfHe.add_gto(0.447,0.314,0,0,0,pos);
    //compute integrals for Helium:
    double He_overlap = overlapCGF(cgfHe,cgfHe);
    double He_kinetic = kineticCGF(cgfHe,cgfHe);
    double He_nuclear = nuclearCGF(cgfHe,cgfHe,He_Z,pos);
    double He_repulsion = two_electronCGF(cgfHe,cgfHe,cgfHe,cgfHe);
    cout << "Testcase for Helium: " << endl;
    cout << "He_overlap: " << He_overlap << endl;
    cout << "He_kinetic: " << He_kinetic << endl;
    cout << "He_nuclear: " << He_nuclear << endl;
    cout << "He_repulsion: " << He_repulsion << endl;
    cout << "He_Energy: " << 2*He_kinetic + 2*He_nuclear + He_repulsion << endl << endl;

    //perform energy minimisation for H2:
    cout << "Testcase for H2: " << endl << endl;

    //H2-molecule
    CGF cgfH_2; //create orbitals for second hydrogen atom at pos2
    cgfH_2.add_gto(0.154,3.425,0,0,0,pos2);
    cgfH_2.add_gto(0.535,0.624,0,0,0,pos2);
    cgfH_2.add_gto(0.445,0.169,0,0,0,pos2);

    vector<CGF> AO_list; //create object containing all atomic orbitals
    AO_list.push_back(cgfH);
    AO_list.push_back(cgfH_2);

    vector<vec3> pos_list; //create list of nucleic positions
    pos_list.push_back(pos);
    pos_list.push_back(pos2);

    vector<double> charge_list; //create list of nucleic charges
    charge_list.push_back(H_Z);
    charge_list.push_back(H_Z);

    vector<int> nelec_list; //create list of electrons for each atom
    nelec_list.push_back(H_nelec);
    nelec_list.push_back(H_nelec);

    SCF_E H2_results; //create output struct for H2 SCF energy minimisation
    H2_results = SCF_HF_energy(AO_list,pos_list,charge_list, nelec_list); //run SCF energy minimisation and get results

    cout << "Program has ended" << endl;
    return 0;
}

//SCF-loop for minimisation of the total energy
SCF_E SCF_HF_energy(vector<CGF> AO_list, vector<vec3> pos_list, vector<double> charge_list, vector<int> nelec_list)
{
    //****************************************************************
    //Hartree-Fock energy minimisation:

    SCF_E out; //initialise output struct

    if(pos_list.size()!=charge_list.size()){
        cout << endl << "Error: number of nuclei not equal to supplied list of charges" << endl << endl;
        return out;
    }

    //get total number of electrons in the system
    int nelec = 0;
    for(int i=0; i<nelec_list.size(); i++) {
        nelec = nelec + nelec_list[i];
    }

    if(nelec % 2 != 0){
        cout << endl << "Error: only closed-shell systems allowed" << endl;
        cout << "nelec = " << nelec << endl << endl;
        return out;
    }

    //initialise matrices for overlap, kinetic, nuclear and repulsion
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //overlap_matrix
    Eigen::MatrixXd Kin = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //kinetic_matrix
    Eigen::MatrixXd Nucl = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //nuclear_matrix

    //calculate respective 1-integrals: overlap, kinetic and nuclear
    for(int i=0; i<AO_list.size(); i++) {
        for(int j=0; j<AO_list.size(); j++) {
            S(i,j) = overlapCGF(AO_list[i], AO_list[j]);
            Kin(i,j) = kineticCGF(AO_list[i], AO_list[j]);
            for(int n=0; n<charge_list.size(); n++){
                Nucl(i,j) = Nucl(i,j) + nuclearCGF(AO_list[i], AO_list[j], charge_list[n], pos_list[n]);
            }
        }
    }
    //get 1e-Hamiltonian
    Eigen::MatrixXd Hamil = Kin + Nucl;

    //calculate all unique 2e-integral using index sorting:
    //
    //vector of unique 2e-integrals, -1 as a reference value for non-calculated integral
    vector<double> Rep(two_electronSort(AO_list.size(),AO_list.size(),AO_list.size(),AO_list.size()), -1);
    int ij,kl,current_index;
    for(int i=0; i<AO_list.size(); i++) {
        for(int j=0; j<AO_list.size(); j++) {
            ij = i*(i+1)/2 + j;
            for(int k=0; k<AO_list.size(); k++) {
                for(int l=0; l<AO_list.size(); l++) {
                    kl = k * (k+1)/2 + l;
                    if(ij <= kl) {
                        current_index = two_electronSort(i,j,k,l);
                        //do not compute the integral if it has previously been evaluated (repulsion cannot be negative, so this is a safe default value)
                        if(Rep[current_index] != -1.0) {
                            continue;
                        }
                        Rep[current_index] = two_electronCGF(AO_list[i], AO_list[j], AO_list[k], AO_list[l]);
                    }
                }
            }
        }
    }

    //perform canonical diagonalisation:
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
    Eigen::MatrixXd Diagonal = es.eigenvalues().real().asDiagonal(); //diagonal matrix
    Eigen::MatrixXd Unitary = es.eigenvectors().real(); //unitary matrix

    //take 1/sqrt(D); element-wise since the matrix is diagonal
    for(int i=0; i<AO_list.size(); i++) {
        Diagonal(i,i) = 1.0 / sqrt(Diagonal(i,i));
    }

    //Calculate the transformation matrix
    Eigen::MatrixXd X_transform = Unitary * Diagonal;

    //matrix initialisation for energy-minimisation:
    Eigen::MatrixXd P_density = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //initialise density-matrix
    Eigen::MatrixXd P_density_new = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //initialise update-matrix for density
    Eigen::MatrixXd G_two_electron = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //Initialise two-electron Hamiltonian matrix
    Eigen::MatrixXd C_vector; //initialise eigenvectors of the Fock-matrix
    Eigen::VectorXd orbital_energies; //initialise eigenvalues of the Fock-matrix

    //optimisation parameters:
    static const double alpha = 0.5; //mixing parameter for updating density
    static const double TolEnergy = 1e-5; //tolerance for energy convergence
    static const double loop_max = 100; //maximum number of iterations

    //initialise update parameters:
    double energy_old = 0.0; //energy obtained from previous iteration for comparison
    double energy_difference = 1.0; //difference between new and old energy
    int loop_counter = 0; //initialise counter
    double energy = 0.0; //initialise energy

    cout << "Beginning electronic optimisation: " << endl;
    cout << "Iter./Energy;" << endl;

    //SCF-loop for energy optimisation:
    while(energy_difference > TolEnergy && loop_counter < loop_max) {
        loop_counter++; //increment loop counter

        //Calculate two-electron hamiltonian matrix:
        static int indexJ, indexK;
        for(int i=0; i<AO_list.size(); i++) {
            for(int j=0; j<AO_list.size(); j++) {
                G_two_electron(i,j) = 0.; //reset matrix
                for(int k=0; k<AO_list.size(); k++) {
                    for(int l=0; l<AO_list.size(); l++) {
                        indexJ = two_electronSort(i,j,l,k);
                        indexK = two_electronSort(i,k,l,j);
                        G_two_electron(i,j) = G_two_electron(i,j) + P_density(k,l) * (Rep[indexJ] - 0.5 * Rep[indexK]);
                    }
                }
            }
        }

        //Calculate Fock Matrix
        Eigen::MatrixXd Fock = Hamil + G_two_electron;

        //Transform Fock Matrix
        Eigen::MatrixXd F_transform = X_transform.transpose() * Fock * X_transform;

        //Calculate eigenvalues within the transformed basis
        es.compute(F_transform);
        Eigen::MatrixXd C_vector_transform = es.eigenvectors().real();
        Eigen::MatrixXd orbital_energies_store = es.eigenvalues().real().asDiagonal(); //store values of the orbital energies

        //Calculate total electronic energy of the system
        energy = 0.0;
        Eigen::MatrixXd E_total = Hamil + Fock; //non-weighted components of total electronic energy
        for(int i=0; i<AO_list.size(); i++) {
            for(int j=0; j<AO_list.size(); j++) {
                energy = energy + 0.5 * P_density(j,i) * E_total(i,j); //weigh matrix elements
            }
        }

        //Add nuclear repulsion to the orbital energies
        for(int i=0; i<charge_list.size(); i++){
            for(int j=0; j<charge_list.size(); j++){
                if(j>i){
                    energy = energy + charge_list[i]*charge_list[j]/(pos_list[i]-pos_list[j]).norm();
                }
            }
        }

        //Transform coefficient matrix back to original basis
        C_vector = X_transform * C_vector_transform;
        orbital_energies = orbital_energies_store.diagonal(); //transform matrix back to vector

        //generate new density matrix
        P_density_new = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size());
        for(int i=0; i<AO_list.size(); i++) {
            for(int j=0; j<AO_list.size(); j++) {
                for(int k=0; k<(nelec/2); k++) {
                    P_density_new(i,j) = P_density_new(i,j) + 2.0 * C_vector(i,k) * C_vector(j,k);
                }
            }
        }

        //update density matrix by mixing the old and new matrices
        for(int i=0; i<AO_list.size(); i++) {
            for(int j=0; j<AO_list.size(); j++) {
                P_density(i,j) = (1.0-alpha) * P_density_new(i,j) + alpha * P_density(i,j);
            }
        }

        //report energy and calculate difference with previous energy value
        energy_difference = abs(energy - energy_old);
        energy_old = energy; //store current value for the energy for comparison in the next iteration
        cout << loop_counter << ":\t" << energy << endl;
    }
    //end SCF-loop

    //output SCF results
    cout << "Stopping because energy convergence was achieved." << endl;
    cout << "System Energy = " << energy << endl << endl;

    //outputting orbital energies and orbital coefficients
    cout << "Orbital Energies:" << endl;
    for(int i=0; i<orbital_energies.size(); i++){
        cout << "e" << i+1 << " = " << orbital_energies[i] << endl;
        cout << "coeff: " << endl << C_vector.col(i) << endl << endl;
    }

    //parse results
    out.SCF_result(energy,orbital_energies,C_vector);
    return out;
}

//end main codeblock
//****************************************************************


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
                for(int v2=0; v2<CGF4.get_size(); v2++){
                    electric = electric + CGF1.coeff[w]*CGF2.coeff[v]*CGF3.coeff[w2]*CGF4.coeff[v2]*two_electronGTO(CGF1.GTOlist[w],CGF1.GTOlist[v],CGF3.GTOlist[w2],CGF4.GTOlist[v2]);
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
