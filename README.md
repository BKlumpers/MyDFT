#origin DFT general:

/*****************************************************************

Basic closed-shell spin-restricted DFT-solver for simple molecules using STO-NG

Authors: B. Klumpers
		 I.A.W. Filot

Published under GNU General Public License 3.0
Source code available at: https://github.com/BKlumpers/dft

Allows for SCF-computation of molecular energies for simple molecules.
Includes testcases for: H, He, H2, HeH+, He2, CO, and H2O.

*****************************************************************/

#branch CGF:

/*******

Development branch; will be updated most often, but may contain missing functionality.

Any updates to this branch do however work as intended.

If you encounter any bugs, please post these as an issue on https://github.com/BKlumpers/dft/

*******/

#Dependencies:

/*****************************************************************

Eigen3

Boost 1.6.5 or later: gamma, spherical_harmonics, cubic_b_spline

(Does not work with earlier versions of boost as cubic_b_spline is not implemented)

*****************************************************************/