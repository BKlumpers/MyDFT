/*****************************************************************

Basic closed-shell spin-restricted HF-solver for simple molecules using STO-NG

Author: B. Klumpers (bartkl@live.nl)

Published under GNU General Public License 3.0

Allows for SCF-computation of molecular energies for simple molecules.
Testcases for H, He, H2, HeH+ and He2 are included.

*****************************************************************/

#include "solvers.h"

//****************************************************************
//Hartree-Fock energy minimisation:

/*SCF-loop for minimisation of the total energy
//
//{AO_list} is a vectors of CGFs comprising the basis set
//{pos_list} is a vector of the coordinates of the nuclei
//{charge_list} is a vector of the nuclear charges
//{nelec_list} is a vector of the number of electrons associated with each atom
*/
SCF_E SCF_HF_energy(vector<CGF> AO_list, const vector<vec3>& pos_list, const vector<double>& charge_list, const vector<int>& nelec_list)
{

    SCF_E out; //initialise output struct

    if(pos_list.size()!=charge_list.size()){
        cout << endl << "Error: number of nuclei not equal to supplied list of charges" << endl << endl;
        return out;
    }

    //get total number of electrons in the system
    int nelec = 0;
    for(int i=0; i<nelec_list.size(); i++) {
        nelec += nelec_list[i];
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
                Nucl(i,j) += nuclearCGF(AO_list[i], AO_list[j], charge_list[n], pos_list[n]);
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
    if(pos_list.size() == 1){ //correct update for monoatomics
        energy_old = 1;
    }
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
                G_two_electron(i,j) = 0; //reset matrix
                for(int k=0; k<AO_list.size(); k++) {
                    for(int l=0; l<AO_list.size(); l++) {
                        indexJ = two_electronSort(i,j,l,k);
                        indexK = two_electronSort(i,k,l,j);
                        G_two_electron(i,j) += P_density(k,l) * (Rep[indexJ] - 0.5 * Rep[indexK]);
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
                energy += 0.5 * P_density(j,i) * E_total(i,j); //weigh matrix elements
            }
        }

        //Add nuclear repulsion to the orbital energies
        for(int i=0; i<charge_list.size(); i++){
            for(int j=0; j<charge_list.size(); j++){
                if(j>i){
                    energy += charge_list[i]*charge_list[j]/(pos_list[i]-pos_list[j]).norm();
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
                    P_density_new(i,j) += 2.0 * C_vector(i,k) * C_vector(j,k);
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

    //testblock for Density Functional Theory:
    //currently (incorrectly) assumes HF-SCF density matrix
    //will be moved to seperate function

    //compute number of electrons from density
    //N = sum(P*S)
    double Ndel = 0;
    for(int i=0; i<AO_list.size(); i++){
        for(int j=0; j<AO_list.size(); j++){
            Ndel += P_density(i,j)*S(i,j);
        }
    }
    cout << "Ndel: " << Ndel << endl;

    double NPdel = 0; //compute nuclear attraction from density
    for(int i=0; i<AO_list.size(); i++){
        for(int j=0; j<AO_list.size(); j++){
            NPdel += P_density(i,j)*Nucl(i,j);
        }
    }
    cout << "NPdel: " << NPdel << endl;

    double NPkin = 0; //isolate kinetic energy
    for(int i=0; i<AO_list.size(); i++){
        for(int j=0; j<AO_list.size(); j++){
            NPkin += P_density(i,j)*Kin(i,j);
        }
    }
    cout << "NPkin: " << NPkin << endl;

    int index2e,in; //Calculate two-electron matrix from density:
    double D_G_two_electron = 0;
    double K_check = 0;
    for(int i=0; i<AO_list.size(); i++) {
        for(int j=0; j<AO_list.size(); j++) {
            for(int k=0; k<AO_list.size(); k++) {
                for(int l=0; l<AO_list.size(); l++) {
                    index2e = two_electronSort(i,j,l,k);
                    in = two_electronSort(i,k,l,j);
                    D_G_two_electron += 0.5 * P_density(i,j) * P_density(k,l) * Rep[index2e];
                    K_check += -0.5 * P_density(i,j) * P_density(k,l) * 0.5 * Rep[in];
                }
            }
        }
    }
    cout << "DG: " << D_G_two_electron << endl;
    cout << "K_XC: " << K_check << endl;

    double NNrep = 0;//get nuclear repulsion
    for(int i=0; i<charge_list.size(); i++){
        for(int j=0; j<charge_list.size(); j++){
            if(j>i){
                NNrep += charge_list[i]*charge_list[j]/(pos_list[i]-pos_list[j]).norm();
            }
        }
    }
    cout << "NNrep: " << NNrep << endl;

    //really inefficient grid (Hardcoded)
    static const double Nstep = 100; //have integral run in 3d{x,y,z} with Nstep+1 discrete steps along each axis
    static const int steps = int(Nstep) + 1; //turn loop-counter into int to support openmp
    static const double SingularTol = 1e-4; //tolerance for identifying singularities
    double Density = 0; //hold value of the local electron density
    double gridsum = 0; //initialise value of the nuclear integral
    double Nsum = 0; //initialise value of the population integral
    double Jsum = 0; //initialise value of the electronic repulsion integral
    double Exchange = 0; //initialise value of the exchange energy integral
    double Correlation = 0; //initialise value of the correlation energy integral
    //define lower and upper bounds for x-axis
    double limset = 5;
    double xmin = -limset, xmax = limset;
    //define lower and upper bounds for y-axis
    double ymin = -limset, ymax = limset;
    //define lower and upper bounds for z-axis
    double zmin = -limset, zmax = limset;

    double dV = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)/pow(Nstep,3); //mesh-unit volume

    double xpos,ypos,zpos;
    vec3 origin; //mesh-point centre

    double xpos2, ypos2, zpos2;
    vec3 origin2; //mesh-point centre for second integral for electronic repulsion

    double gridDenom = 0; //nuclear integral operator
    double JDenom = 0; //two electron integral operator

    //commence discretised integration
    cout << endl << "Commencing density-based nuclear integration" << endl;
    cout << "This may take a while..." << endl;

    for(int i=0; i<steps; i++){
        xpos = xmin + i*(xmax-xmin)/Nstep;
        for(int j=0; j<steps; j++){
            ypos = ymin + j*(ymax-ymin)/Nstep;
            for(int k=0; k<steps; k++){
                zpos = zmin + k*(zmax-zmin)/Nstep;
                origin << xpos,ypos,zpos;
                Density = density(origin, P_density, AO_list);
                for(int centre=0; centre<charge_list.size(); centre++){ //loop over each nucleus
                    gridDenom = (origin-pos_list[centre]).norm();
                    if(gridDenom > SingularTol){ //avoid singularities
                        gridsum += -charge_list[centre]*dV*Density/gridDenom; //add normalised local value to nuclear integral
                    }
                }
                Nsum += dV*Density; //add normalised local value to population integral
                Exchange += dV*exchange_Dirac(Density); //Dirac exchange functional
                Correlation += dV*correlation_VWN(Density); //Vosko-Wilk-Nusair correlation functional
                //evaluate second integral in the two-electron integration
                //currenty commented since for Nstep = 20; comp. time = ca. 2hrs
                //since above loop = ca. 0.85sec for 8000 evaluations -> 8000*8000 evaluations for current loop -> 8000*0.85sec = ca. 2hrs
                // for(int i2=0; i2<steps; i2++){
                //     xpos2 = xmin + i2*(xmax-xmin)/Nstep;
                //     for(int j2=0; j2<steps; j2++){
                //         ypos2 = ymin + j2*(ymax-ymin)/Nstep;
                //         for(int k2=0; k2<steps; k2++){
                //             zpos2 = zmin + k2*(zmax-zmin)/Nstep;
                //             origin2 << xpos2,ypos2,zpos2;
                //             JDenom = (origin-origin2).norm();
                //             if(JDenom > SingularTol){ //avoid singularities
                //                 Jsum += 0.5*dV*dV*Density*density(origin2, P_density, AO_list)/JDenom;
                //             }
                //         }
                //     }
                // }
            }
        }
    }
    //output results and compare HF with DFT:
    cout << "gridsum: " << gridsum << endl;
    cout << "Nsum: " << Nsum << endl;
    //cout << "Jsum: " << Jsum << endl; //currently commented since not computed in above loop
    cout << "Exchange: " << Exchange << endl;
    cout << "Correlation: " << Correlation << endl;
    double E_DFT = NPkin + D_G_two_electron + Exchange + Correlation + NNrep + NPdel;
    double EHF = NPkin + D_G_two_electron + NPdel + NNrep + K_check;
    cout << endl << "E_DFT: " << E_DFT << " vs(HF): " << EHF << endl << endl;

    //parse results
    out.SCF_result(energy,orbital_energies,C_vector);
    return out;
}


//****************************************************************
//Density Functional Theory:


//compute value for the electron-density at a certain point {pos}
//given the density matrix {Pmatrix} and the basis set {AO_list}
double density(const vec3& pos, const Eigen::MatrixXd& Pmatrix, vector<CGF> AO_list)
{
    double density_value = 0;
    for(int i=0; i<Pmatrix.rows(); i++){
        for(int j=0; j<Pmatrix.rows(); j++){
            density_value += Pmatrix(i,j)*AO_list[i].getvalue(pos)*AO_list[j].getvalue(pos);
        }
    }
    return density_value;
}

//calculate localised exchange-part of the DFT exchange-correlation term
//according to the Dirac expression for a Homogeneous Electron Gas (LDA)
double exchange_Dirac(double Density)
{
    static const double pi = 3.14159265359;
    static const double prefactor = -0.75*pow(3/pi,1/3);
    return prefactor*pow(Density,4/3);
}

//calculate correlation-part of the DFT exchange-correlation term
//according to the Vosko-Wilk-Nusair correlation functional
double correlation_VWN(double Density)
{
    static const double pi = 3.14159265359;
    //Vosko-Wilk-Nusair V parameters
    static const double A_p = 0.0621814;
    static const double x0_p = -0.10498;
    static const double b_p = 3.72744;
    static const double c_p = 12.9352;
    double x_small = pow(3/(4*pi*Density),1/6); //x = sqrt(r0); r0 = (3/(4*pi*Density))^1/3
    double X_large_p = x_small*x_small + x_small*b_p + c_p;
    double Q_p = sqrt(4*c_p - b_p*b_p);
    double Fx_p = log(pow(x_small-x0_p,2)/X_large_p) + 2*(b_p+2*x0_p)*atan(Q_p/(b_p+2*x_small))/Q_p;
    //multiply exchange potential with Density to get local value of exchange energy, prefactor 0.5 to convert from Rydberg to Hartree
    return 0.5*A_p*(    log(x_small*x_small/X_large_p) + 2*b_p*atan(Q_p/(2*x_small+b_p))/Q_p - b_p*x0_p*Fx_p/(x0_p*x0_p + x0_p*b_p + c_p)    )*Density;
}

//End of file
//****************************************************************
