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
        cout << endl << "Error: only spin-restricted systems allowed" << endl;
        cout << "nelec = " << nelec << endl << endl;
        return out;
    }

    //initialise matrices for overlap, kinetic and nuclear Fock-matrix components
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //overlap_matrix
    Eigen::MatrixXd Kin = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //kinetic_matrix
    Eigen::MatrixXd Nucl = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //nuclear_matrix

    //calculate respective 1e-integrals: overlap, kinetic and nuclear
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
    vector<double> Rep(two_electronSort(AO_list.size(),AO_list.size(),AO_list.size(),AO_list.size()), -1.0);
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
    static const int loop_max = 100; //maximum number of iterations

    //initialise update parameters:
    double energy_old = 0.0; //energy obtained from previous iteration for comparison
    if(pos_list.size() == 1){ //correct update for monoatomics
        energy_old = 1.0;
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

    //parse results
    out.SCF_result(energy,orbital_energies,C_vector);
    return out;
}

//End of file
//****************************************************************
