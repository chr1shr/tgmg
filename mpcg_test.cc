#include <cmath>
#include <cstring>
#include <limits>

#include "tgmg.hh"
#include "mpcg_fem.hh"

// Grid dimensions
const int m=65,n=65;

// Total number of gridpoints
const int mn=m*n;

// Physical dimensions of the grid
const double ax=-3,bx=3,ay=ax,by=bx;

// Grid spacings
const double dx=(bx-ax)/(m-1),dy=(by-ay)/(n-1);

int main() {

    // Total number of trials
    const int num_trials=3;

    // Total number of complete spins that the Gaussian source term is rotated
    // through. From this, define the angular step size in radians.
    const double spins=1;
    const double h=2.*M_PI*spins/num_trials;

    // Create the class that can perform the multigrid-preconditioned conjugate
    // gradient (MPCG) algorithm. It contains a TGMG solver within it.
    mpcg_fem mf(m,n,false,false,ax,dx,ay,dy);
    mf.clear_x();

    // Do a number of trials where the source term for the solve is slowly
    // moved around the domain
    for(int k=0;k<num_trials;k++) {

        // Set up the source term to be the sum of two Gaussians. Since Neumann
        // conditions are used, the source term should integrate to zero.
        mf.clear_b();
        mf.add_gaussian_source_term(1.5*cos(h*k),1.5*sin(h*k),0.3,1);
        mf.add_gaussian_source_term(-1.5*cos(h*k),-1.5*sin(h*k),0.3,-1);

        // Solve using the MPCG algorithm
        mf.solve_pre_cg(true);putchar('\n');
    }

    mf.save_fields();
}
