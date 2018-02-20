/*****************************************************************

Basic closed-shell spin-restricted HF-solver for simple molecules using STO-NG

Author: B. Klumpers (bartkl@live.nl)

Published under GNU General Public License 3.0

Allows for SCF-computation of molecular energies for simple molecules.
Testcases for H, He, H2, HeH+ and He2 are included.

*****************************************************************/

#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "integrals.h"
#include "orbitals.h"
#include "solvers.h"

typedef Eigen::Vector3d vec3; //define vec3 as the Eigen::Vector3d-object

using namespace std;

const bool visual = false; //process visualisation of orbitals

void plot2D(char const title[], double height, double width, double xmin, double xmax, double zmin, double zmax, vector<CGF> AO_list, const Eigen::MatrixXd& Coeff);

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
    cgfH.add_gto(0.15432897000000001,3.4252509099999999,0,0,0,pos);
    cgfH.add_gto(0.53532813999999995,0.62391373000000006,0,0,0,pos);
    cgfH.add_gto(0.44463454000000002,0.16885539999999999,0,0,0,pos);
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
    cgfHe.add_gto(0.15432897000000001,6.3624213899999997,0,0,0,pos);
    cgfHe.add_gto(0.53532813999999995,1.1589229999999999,0,0,0,pos);
    cgfHe.add_gto(0.44463454000000002,0.31364978999999998,0,0,0,pos);
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

    //check SCF for single atom
    // vector<CGF> AO_list1; //create object containing all atomic orbitals
    // AO_list1.push_back(cgfHe);
    // vector<vec3> pos_list1; //create list of nucleic positions
    // pos_list1.push_back(pos);
    // vector<double> charge_list1; //create list of nucleic charges
    // charge_list1.push_back(He_Z);
    // vector<int> nelec_list1; //create list of electrons for each atom
    // nelec_list1.push_back(He_nelec);
    // cout << "SCF test: " << endl;
    // SCF_E He1_results; //create output struct for H2 SCF energy minimisation
    // He1_results = SCF_HF_energy(AO_list1,pos_list1,charge_list1, nelec_list1); //run SCF energy minimisation and get results

    //perform energy minimisation for H2:
    cout << "Testcase for H2: " << endl << endl;

    //H2-molecule
    CGF cgfH_2; //create orbitals for second hydrogen atom at pos2
    cgfH_2.add_gto(0.15432897000000001,3.4252509099999999,0,0,0,pos2);
    cgfH_2.add_gto(0.53532813999999995,0.62391373000000006,0,0,0,pos2);
    cgfH_2.add_gto(0.44463454000000002,0.16885539999999999,0,0,0,pos2);

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

    //visualisation of H2 ground-state orbital
    
    //define the domain which will be visualised
    double width = 500; //number of gridpoints in the x-direction
    double height = 500; //number of gridpoints in the z-direction
    
    //x runs from -3,+3 (sigma axis)
    double xmin = -3, xmax = 3;
    //z runs -3/+3 apart from each atom (internuclear axis)
    double zmin = -3, zmax = 3+1.4;

    //write image of groundstate MO for H2 to file
    if(visual){
        plot2D("H2.txt",height,width,xmin,xmax,zmin,zmax,AO_list,H2_results.C_vector);
    }

    //HeH+ molecule
    cout << "Testcase for HeH+: " << endl << endl;

    //position hydrogen appropriate distance from helium:
    vec3 pos_HeHp;
    pos_HeHp << 0,0,0.772;
    CGF cgfH_He; //create orbitals for second hydrogen atom at pos2
    cgfH_He.add_gto(0.15432897000000001,3.4252509099999999,0,0,0,pos_HeHp);
    cgfH_He.add_gto(0.53532813999999995,0.62391373000000006,0,0,0,pos_HeHp);
    cgfH_He.add_gto(0.44463454000000002,0.16885539999999999,0,0,0,pos_HeHp);

    vector<CGF> AO_list_HeHp; //create object containing all atomic orbitals
    AO_list_HeHp.push_back(cgfHe);
    AO_list_HeHp.push_back(cgfH_He);

    vector<vec3> pos_list_HeHp; //create list of nucleic positions
    pos_list_HeHp.push_back(pos);
    pos_list_HeHp.push_back(pos_HeHp);

    vector<double> charge_list_HeHp; //create list of nucleic charges
    charge_list_HeHp.push_back(He_Z);
    charge_list_HeHp.push_back(H_Z);

    vector<int> nelec_list_HeHp; //create list of electrons for each atom
    nelec_list_HeHp.push_back(He_nelec-1); //He+ has one electron less
    nelec_list_HeHp.push_back(H_nelec);

    //perform energy minimisation for HeH+:
    SCF_E HeHp_results;
    HeHp_results = SCF_HF_energy(AO_list_HeHp,pos_list_HeHp,charge_list_HeHp,nelec_list_HeHp);

    if(visual){
        zmax = 3+0.772; //rescale z-axis to HeH+ bond length, then write image of ground-state MO
        plot2D("HeHp.txt",height,width,xmin,xmax,zmin,zmax,AO_list_HeHp,HeHp_results.C_vector);
    }

    //He2 molecule
    cout << "Testcase for He2: " << endl << endl;

    //position second helium appropriate distance from first helium:
    vec3 pos_He2;
    pos_He2 << 0,0,5.2;
    CGF cgfHe2; //create orbitals for second helium atom at pos2
    cgfHe2.add_gto(0.15432897000000001,6.3624213899999997,0,0,0,pos_He2);
    cgfHe2.add_gto(0.53532813999999995,1.1589229999999999,0,0,0,pos_He2);
    cgfHe2.add_gto(0.44463454000000002,0.31364978999999998,0,0,0,pos_He2);

    vector<CGF> AO_list_He2; //create object containing all atomic orbitals
    AO_list_He2.push_back(cgfHe);
    AO_list_He2.push_back(cgfHe2);

    vector<vec3> pos_list_He2; //create list of nucleic positions
    pos_list_He2.push_back(pos);
    pos_list_He2.push_back(pos_He2);

    vector<double> charge_list_He2; //create list of nucleic charges
    charge_list_He2.push_back(He_Z);
    charge_list_He2.push_back(He_Z);

    vector<int> nelec_list_He2; //create list of electrons for each atom
    nelec_list_He2.push_back(He_nelec);
    nelec_list_He2.push_back(He_nelec);

    //perform energy minimisation for He2:
    SCF_E He2_results;
    He2_results = SCF_HF_energy(AO_list_He2,pos_list_He2,charge_list_He2,nelec_list_He2);

    if(visual){
        zmax = 3+5.2; //rescale to He2 bond length, then write image of ground-state MO to file
        plot2D("He2.txt",height,width,xmin,xmax,zmin,zmax,AO_list_He2,He2_results.C_vector);
    }

    cout << "Program has ended" << endl;
    return 0;
}

//end main codeblock
//****************************************************************

//generate an image containing the amplitudes of the ground-state molcular orbital
//this image can be used to generate a colourmap of the MO in an appropriate visualisation tool
//{title} will be the filename
//{AO_list} is a vector of the atomic orbitals
//{Coeff} are the optimised orbital coefficients
//remaining parameters define the domain in which the MO is visualised
void plot2D(char const title[], double height, double width, double xmin, double xmax, double zmin, double zmax, vector<CGF> AO_list, const Eigen::MatrixXd& Coeff)
{
    //visualise the xz-plane
    //
    vec3 origin; //mesh-point centre
    double xpos,ypos,zpos; //initialise position vector components
    ypos = 0; //visualise the plane corresponding to y=0
    int orbital_vis = 0; //visualise ground-state orbital
    Eigen::MatrixXd image = Eigen::MatrixXd::Zero(width+1,height+1); //initialise matrix storing pixels

    cout << endl << "Commencing writing to file: " << title << endl;

    //compute amplitude at each point
    for(int i=0; i<(width+1); i++){
        xpos = xmin + i*(xmax-xmin)/width;
        for(int j=0; j<(height+1); j++){
            zpos = zmin + j*(zmax-zmin)/height;
            //j loops across z, i loops across x
            origin << xpos,ypos,zpos;
            for(int orbs=0; orbs<AO_list.size(); orbs++){ //loop across all atomic orbitals
                image(i,j) += AO_list[orbs].getvalue(origin)*Coeff(orbs,orbital_vis);
            }
        }
    }
    // if(string::empty(title)){ //requires #include <string>
    //     title = "2D_plot.txt"; //default filename
    // }
    std::ofstream file(title); //set file as output stream
    file << image << '\n'; //write image to file
    cout << "Finished writing." << endl << endl;
    return;
}

//End of file
//****************************************************************
