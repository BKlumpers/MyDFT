/*****************************************************************

Basic closed-shell spin-restricted DFT-solver for simple molecules using STO-NG

Authors: B. Klumpers
         I.A.W. Filot

Published under GNU General Public License 3.0
Source code available at: https://github.com/BKlumpers/dft

Allows for SCF-computation of molecular energies for simple molecules.
Includes testcases for: H, He, H2, HeH+, He2, CO, and H2O.

*****************************************************************/

#include "dft.h"

//****************************************************************
//DFT energy minimisation:

/*SCF-loop for minimisation of the total energy
//
//{AO_list} is a vectors of CGFs comprising the basis set
//{pos_list} is a vector of the coordinates of the nuclei
//{charge_list} is a vector of the nuclear charges
//{nelec_list} is a vector of the number of electrons associated with each atom
*/
SCF_E SCF_DFT_energy(vector<CGF> AO_list, const vector<vec3>& pos_list, const vector<double>& charge_list, const vector<int>& nelec_list, const vector<int>& atnum_list)
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
    vector<double> Rep(two_electronSort(AO_list.size(), AO_list.size(), AO_list.size(), AO_list.size()), -1.0);
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

    //really inefficient grid (Hardcoded)
    static const double Nstep = 30.0; //have integral run in 3d{x,y,z} with Nstep+1 discrete steps along each axis
    static const int steps = int(Nstep) + 1; //turn loop-counter into int to support openmp
    static const double SingularTol = 1e-4; //tolerance for identifying singularities
    static const double DenTol = 1e-10;
    double Density = 0.0; //hold value of the local electron density
    double Local = 0.0; //hold value of local AO-density
    double weight = 0.0; //hold value of the weight of the quadratic grid
    Eigen::MatrixXd Exchange = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //initialise value of the exchange energy integral
    Eigen::MatrixXd Correlation = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //initialise value of the correlation energy integral
    double Ex = 0.0; //initialise exchange energy
    double Ec = 0.0; //initialise correlation energy
    double NELEC = 0.0; //initialise counter for electrons
    //define lower and upper bounds for x-axis
    double limset = 3.0;
    double xmin = -limset, xmax = limset;
    //define lower and upper bounds for y-axis
    double ymin = -limset, ymax = limset;
    //define lower and upper bounds for z-axis
    double zmin = -limset, zmax = limset;

    double dV = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)/pow(Nstep,3.0); //mesh-unit volume

    double xpos,ypos,zpos;
    vec3 origin; //mesh-point centre

    //Gauss-Chebyshev Lebedev grid:
    cout << "Generate GCL-grid: " << endl;
    grid GCL_grid = make_grid(pos_list, atnum_list, AO_list.size());
    cout << "Writing wavefunctions: " << AO_list.size() << endl;
    //add line for setting wavefunction amplitudes
    GCL_grid.write_amps(AO_list);
    cout << "Finished generating grid;" << endl;

    cout << "Beginning electronic optimisation: " << endl;
    cout << "Iter./Energy;" << endl;


    Eigen::MatrixXd Poisson = Exchange;
    double Jp1 = 0.0;
    double Jx = 0.0;


    //SCF-loop for energy optimisation:
    while(energy_difference > TolEnergy && loop_counter < loop_max) {
        loop_counter++; //increment loop counter

        Ex = 0.0; //reset energies
        Ec = 0.0;
        NELEC = 0.0;
        Jp1 = 0.0;
        //Calculate two-electron hamiltonian matrix:
        static int indexJ;
        for(int i=0; i<AO_list.size(); i++) {
            for(int j=0; j<AO_list.size(); j++) {
                G_two_electron(i,j) = 0.0; //reset matrices
                //compute exchange-correlation according to a Cartesian grid
                if(!Quadrature){
                    Exchange(i,j) = 0.0;
                    Correlation(i,j) = 0.0;
                    for(int xx=0; xx<steps; xx++){
                        xpos = xmin + xx*(xmax-xmin)/Nstep;
                        for(int yy=0; yy<steps; yy++){
                            ypos = ymin + yy*(ymax-ymin)/Nstep;
                            for(int zz=0; zz<steps; zz++){
                                zpos = zmin + zz*(zmax-zmin)/Nstep;
                                origin << xpos,ypos,zpos;
                                Density = density(origin, P_density, AO_list);
                                if(Density > DenTol){
                                    Local = AO_list[i].getvalue(origin)*AO_list[j].getvalue(origin);
                                    Exchange(i,j) += dV*exchange_Dirac(Density)*Local; //Dirac exchange functional
                                    Correlation(i,j) += dV*VWN_potential(Density)*Local; //Vosko-Wilk-Nusair correlation functional
                                    if(i == 0 && j == 0){
                                        //the additional nelec/2 is a scaling factor in restricted systems due to present looping scheme
                                        Ex += 0.75*dV*exchange_Dirac(Density)*Density*(nelec/2); //add factor 3/4 to convert potential to energy
                                        Ec += dV*correlation_VWN(Density)*Density*(nelec/2);
                                        NELEC += dV*Density; //check if nelec-factor is needed here as well
                                    }
                                }
                            }
                        }
                    }
                }
                for(int k=0; k<AO_list.size(); k++) {
                    for(int l=0; l<AO_list.size(); l++) {
                        indexJ = two_electronSort(i,j,l,k);
                        G_two_electron(i,j) += 0.5 * P_density(k,l) * Rep[indexJ]; //for DFT only J, not K
                    }
                }
            }
        }

        //compute exchange-correlation terms utilising numerical quadrature
        if(Quadrature){
            Exchange = Eigen::MatrixXd::Zero(Exchange.rows(), Exchange.cols()); //reset XC-matrices
            Correlation = Eigen::MatrixXd::Zero(Correlation.rows(), Correlation.cols());
            Poisson = Exchange;
            GCL_grid.set_density(P_density); //update electron density map
            for(int rpnt=0; rpnt<GCL_grid.get_rpnt(); rpnt++){ //loop across all radial points
                for(int apnt=0; apnt<GCL_grid.get_apnt(); apnt++){ //loop across all angular points
                    for(int atom=0; atom<GCL_grid.get_atomnr(); atom++){ //loop across all atomic grids
                        //get electron density at this gridpoint
                        Density = GCL_grid.get_density(atom, rpnt, apnt);//density(atom, rpnt, apnt, P_density, GCL_grid);
                        if(Density > DenTol){
                            weight = GCL_grid.get_weight(atom,rpnt,apnt);
                            Ex += 0.75*weight*exchange_Dirac(Density)*Density; //add factor 3/4 to convert potential to energy
                            Ec += weight*correlation_VWN(Density)*Density;
                            NELEC += weight*Density;
                            //Jp1 += 0.5 * weight * Density * GCL_grid.get_Poisson(atom,rpnt);
                            for(int i=0; i<AO_list.size(); i++){ //change nesting order, which limits the number of times the density must be computed
                                for(int j=0; j<AO_list.size(); j++){
                                    Local = GCL_grid.get_amps(atom,rpnt,apnt,i)*GCL_grid.get_amps(atom,rpnt,apnt,j);
                                    Exchange(i,j) += weight*exchange_Dirac(Density)*Local; //Dirac exchange functional
                                    Correlation(i,j) += weight*VWN_potential(Density)*Local; //Vosko-Wilk-Nusair correlation functional
                                    //Poisson(i,j) += 0.5*weight*Local*GCL_grid.get_Poisson(atom,rpnt);
                                }
                            }
                        }
                    }
                }
            }
        }
        // double denden=0.0;
        // for(int i=0; i<GCL_grid.get_atomnr(); i++){
        //     cout << "density: " << GCL_grid.get_point_charge(i) << endl;
        //     denden += GCL_grid.get_point_charge(i);
        // }
        // cout <<"ttl " << denden << "    " << NELEC << endl;
        //cout << "Poisson:" << endl << Poisson << endl;
        //Calculate Fock Matrix
        Eigen::MatrixXd Fock = Hamil + 2.0*G_two_electron + Exchange + Correlation;
        
        //Transform Fock Matrix
        Eigen::MatrixXd F_transform = X_transform.transpose() * Fock * X_transform;

        //Calculate eigenvalues within the transformed basis
        es.compute(F_transform);
        Eigen::MatrixXd C_vector_transform = es.eigenvectors().real();
        Eigen::MatrixXd orbital_energies_store = es.eigenvalues().real().asDiagonal(); //store values of the orbital energies

        //Calculate total electronic energy of the system
        energy = 0.0;
        // double T = 0.0;
        Jx = 0.0;
        // double Nx = 0.0;
        // double Pot = 0.0;
        //double Jp2 = 0.0;
        double Pois = 0.0;

        Eigen::MatrixXd E_total = Hamil + G_two_electron; //non-weighted components of total electronic energy
        for(int i=0; i<AO_list.size(); i++) {
            for(int j=0; j<AO_list.size(); j++) {
                energy += P_density(j,i) * E_total(i,j); //weigh matrix elements
                // T += P_density(j,i) * Kin(i,j);
                Jx += P_density(j,i) * G_two_electron(i,j);
                //Jp2 += P_density(j,i) * Poisson(i,j);
                // Pot += P_density(j,i) * (Exchange(i,j) + Correlation(i,j));
                //Pois += P_density(j,i) * Poisson(i,j);
            }
        }
        energy += (Ex + Ec); //include XC-energies

        // cout << endl << "Exc: " << (Ex + Ec)*(nelec/2) << endl;
        // cout << "#e: " << NELEC << endl;
        // cout << "Vxc: " << Pot << endl;
        // cout << "T: " << T << endl;
        //cout << "J: " << Jx << " vs " << Jp1 << " vs " << Pois << endl; // " vs " << Jp1 << endl;
        // cout << "Nucl: " << Nx << endl;
        // cout << "One: " << T+Nx << endl << endl;

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

    //debugging:
    Jp1 = 0.0;
    NELEC = 0.0;
    GCL_grid.set_density(P_density);
    GCL_grid.set_Poisson();
    for(int rpnt=0; rpnt<GCL_grid.get_rpnt(); rpnt++){ //loop across all radial points
        for(int apnt=0; apnt<GCL_grid.get_apnt(); apnt++){ //loop across all angular points
            for(int atom=0; atom<GCL_grid.get_atomnr(); atom++){ //loop across all atomic grids
                //get electron density at this gridpoint
                Density = GCL_grid.get_density(atom, rpnt, apnt);//density(atom, rpnt, apnt, P_density, GCL_grid);
                if(Density > DenTol){
                    weight = GCL_grid.get_weight(atom,rpnt,apnt);
                    NELEC += weight*Density;
                    Jp1 += 0.5 * weight * Density * GCL_grid.get_Poisson(atom,rpnt,apnt);
                    for(int i=0; i<AO_list.size(); i++){ //change nesting order, which limits the number of times the density must be computed
                        for(int j=0; j<AO_list.size(); j++){
                            Local = GCL_grid.get_amps(atom,rpnt,apnt,i)*GCL_grid.get_amps(atom,rpnt,apnt,j);
                            Poisson(i,j) += 0.5*weight*Local*GCL_grid.get_Poisson(atom,rpnt,apnt);
                        }
                    }
                }
            }
        }
    }
    cout << "Results from Poisson integration:" << endl << "WARNING: accuracy of this method is not consistent." << endl;
    cout << "Density on atom 1:" << GCL_grid.get_point_charge(0) << endl;
    //cout << "Density on atom 2:" << GCL_grid.get_point_charge(1) << endl;
    cout << "P_density: " << endl << P_density << endl;
    cout << "Poisson: " << endl << Poisson << endl;
    cout << "J: " << endl << G_two_electron << endl;
    cout << "Total: " << endl << Jp1 << " vs " << Jx << endl;
    cout << "PC check: " << endl;
    GCL_grid.PC_check(0);
    //GCL_grid.PC_check(1);
    // double Tkin = 0.0;
    // double Vnucl = 0.0;
    // double JJ = 0.0;
    // double ExP = 0.0;
    // double EcP = 0.0;
    // for(int i=0; i<AO_list.size(); i++) {
    //     for(int j=0; j<AO_list.size(); j++) {
    //         Tkin += P_density(j,i)*Kin(i,j);
    //         Vnucl += P_density(j,i)*Nucl(i,j);
    //         JJ += P_density(j,i)*G_two_electron(i,j);
    //         ExP += P_density(j,i)*Exchange(i,j);
    //         EcP += P_density(j,i)*Correlation(i,j);
    //     }
    // }
    // //Add nuclear repulsion to the orbital energies
    // double NN = 0.0;
    // for(int i=0; i<charge_list.size(); i++){
    //     for(int j=0; j<charge_list.size(); j++){
    //         if(j>i){
    //             NN += charge_list[i]*charge_list[j]/(pos_list[i]-pos_list[j]).norm();
    //         }
    //     }
    // }
    // cout << "Tkin: " << Tkin << endl;
    // cout << "Vnucl: " << Vnucl << endl;
    // cout << "JJ: " << JJ << endl;
    // cout << "Ex: " << Ex << " vs ExP: " << ExP << endl;
    // cout << "Ec: " << Ec << " vs EcP: " << EcP << endl;
    // cout << "NN: " << NN << endl;
    // cout << "TTL: " << Tkin+Vnucl+JJ+Ex+Ec+NN << endl;

    //debugging_end

    //compute number of electrons from density
    //N = sum(P*S)
    // double Ndel = 0.0;
    // for(int i=0; i<AO_list.size(); i++){
    //     for(int j=0; j<AO_list.size(); j++){
    //         Ndel += P_density(i,j)*S(i,j);
    //     }
    // }

    // vec3 posQ;
    // double Quad = 0.0;
    // //compute number of electrons using quadrature:
    // for(int atom=0; atom<pos_list.size(); atom++){
    //     for(int rpnt=0; rpnt<15; rpnt++){
    //         for(int apnt=0; apnt<110; apnt++){
    //             posQ << GCL_grid.get_gridpoint(atom,rpnt,apnt); //the 3 iteration loops go over all the gridpoints in the molecular grid
    //             Quad += GCL_grid.get_weight(atom,rpnt,apnt)*density(posQ, P_density, AO_list);
    //         }
    //     }
    // }
    
    //output SCF results
    if(loop_counter == loop_max){
        cout << "Stopping because maximum number of iterations was reached." << endl;
        cout << "Energy may not have converged." << endl;
    }
    else{
        cout << "Stopping because energy convergence was achieved." << endl;
    }
    cout << "System Energy = " << energy << endl;
    cout << "#Electrons = " << NELEC << endl;
    //cout << "#Electrons = " << Ndel << " vs SCF: " << NELEC << " vs Quad: " << Quad << endl << endl;

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
    double density_value = 0.0;
    for(int i=0; i<Pmatrix.rows(); i++){
        for(int j=0; j<Pmatrix.rows(); j++){
            density_value += Pmatrix(i,j)*AO_list[i].getvalue(pos)*AO_list[j].getvalue(pos);
        }
    }
    return density_value;
}
//**************************
double density(int atom, int rpnt, int apnt, const Eigen::MatrixXd& Pmatrix, grid Grid) //overloading to include quadrature
{
    double density_value = 0.0;
    for(int i=0; i<Pmatrix.rows(); i++){
        for(int j=0; j<Pmatrix.rows(); j++){
            density_value += Pmatrix(i,j)*Grid.get_amps(atom,rpnt,apnt,i)*Grid.get_amps(atom,rpnt,apnt,j);
        }
    }
    return density_value;
}

//calculate localised exchange-part of the DFT exchange-correlation potential
//according to the Dirac expression for a Homogeneous Electron Gas (LDA)
double exchange_Dirac(double Density)
{
    static const double pi = 3.14159265359;
    static const double prefactor = -0.75*pow(3.0/pi,1.0/3.0)*4.0/3.0; //factor 4/3 to convert energy_prefactor to potential
    double rho = pow(Density,1.0/3.0);
    return prefactor*rho;
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
    double pre = 3.0/(4.0*pi*Density);
    double x_small = pow(pre,1.0/6.0); //x = sqrt(r0); r0 = (3/(4*pi*Density))^1/3
    double X_large_p = x_small*x_small + x_small*b_p + c_p;
    double Q_p = sqrt(4.0*c_p - b_p*b_p);
    double Fx_p = log(pow(x_small-x0_p,2.0)/X_large_p) + 2.0*(b_p+2.0*x0_p)*atan(Q_p/(b_p+2.0*x_small))/Q_p;
    //multiply exchange potential with Density to get local value of exchange energy, prefactor 0.5 to convert from Rydberg to Hartree
    return 0.5*A_p*(    log(x_small*x_small/X_large_p) + 2.0*b_p*atan(Q_p/(2*x_small+b_p))/Q_p - b_p*x0_p*Fx_p/(x0_p*x0_p + x0_p*b_p + c_p)    );
}

//calculate correlation-potential as the functional derivative
//of the Vosko-Wilk-Nusair correlation functional
double VWN_potential(double Density)
{
    static const double pi = 3.14159265359;
    //Vosko-Wilk-Nusair V parameters
    static const double A_p = 0.0621814;
    static const double x0_p = -0.10498;
    static const double b_p = 3.72744;
    static const double c_p = 12.9352;
    double pre = 3.0/(4.0*pi*Density);
    double x_small = pow(pre,1.0/6.0); //x = sqrt(r0); r0 = (3/(4*pi*Density))^1/3
    double X_large_p = x_small*x_small + x_small*b_p + c_p;
    double Q_p = sqrt(4.0*c_p - b_p*b_p);
    double dF = 2.0/(x_small-x0_p)-(2.0*x_small+b_p)/X_large_p-4.0*(2.0*x0_p+b_p)/(pow(2.0*x_small+b_p,2.0)+Q_p*Q_p);

    double deriv = 0.5*A_p*(2.0/x_small - (2.0*x_small+b_p)/X_large_p - 4.0*b_p/(pow(2.0*x_small+b_p,2.0)+Q_p*Q_p) - (b_p*x0_p/(x0_p*x0_p + x0_p*b_p + c_p)) * dF);
    return correlation_VWN(Density) - x_small*deriv/6.0;
}

//End of file
//****************************************************************
