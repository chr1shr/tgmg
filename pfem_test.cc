#include <cmath>
#include <algorithm>
#include <cstring>
#include <limits>

#include "tgmg.hh"
#include "poisson_fem.hh"

// Grid dimensions
const int m=65,n=65;

// Total number of gridpoints
const int mn=m*n;

// Physical dimensions of the grid
const double ax=-3,bx=3,ay=ax,by=bx;

// Grid spacings
const double dx=(bx-ax)/(m-1),dy=(by-ay)/(n-1);

int main() {

    // Create the problem class
    poisson_fem pfem(m,n,false,false,true,dx,dy);

    // Creat the multigrid class
    double *b=pfem.b,*z=pfem.z;
    tgmg<poisson_fem,double,double> mg(pfem,b,z);
    mg.verbose=3;

    // Set up the multigrid hierarchy
    mg.setup();

    // Set up the source term. Since it uses Neumann boundary conditions the
    // source term must sum to zero.
    for(int ij=0;ij<mn;ij++) b[ij]=z[ij]=0;
    b[m/4+n/4*m]=1;
    b[3*m/4+3*n/4*m]=-1;

    // Solve the system using multigrid V-cycles
    tgmg_predict tp;
    mg.solve_v_cycle(tp,1,1,10);

    // Output the solutions in a format that can be read by Gnuplot using
    // the command "splot 'filename' matrix binary"
    mg.output_b("src.0",ax,dx,ay,dy);
    mg.output_z("sol.0",ax,dx,ay,dy);
}
