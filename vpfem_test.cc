#include <cmath>
#include <algorithm>
#include <cstring>
#include <limits>

#include "vpoiss_fem.hh"

// Grid dimensions
const int m=65,n=65;

// Physical dimensions of the grid
const double ax=-3,bx=3,ay=ax,by=bx;

// Grid spacings
const double dx=(bx-ax)/(m-1),dy=(by-ay)/(n-1);

int main() {
    vpoiss_fem vpf(m,n,false,false,true,dx,dy);
    int ij,mn=vpf.mn;
    double *b=vpf.b,*z=vpf.z;

    // Set up the source term
//    for(int i=0;i<mn;i++) {
//        b[i]=1;
//    }
    b[32*m+32]=1;
    b[40*m+40]=-1;

    // Set the initial guess for the solution
    for(ij=0;ij<mn;ij++) z[ij]=0.;

    // Set a block of the prefactor values to 2
    vpf.set_pre_block(0,m,0,30,20);

    // Initialize the multigrid solver, set the verbosity to maximum, and set
    // up the multigrid hierarchy
    tgmg<vpoiss_fem,double,double> mg(vpf,vpf.b,vpf.z);
    mg.setup();
    mg.verbose=3;

    // Solve the linear system using the multigrid library.
    mg.solve_v_cycle();

    // Output the solutions in a format that can be read by Gnuplot using
    // the command "splot 'filename' matrix binary"
    mg.output_b("b.0",ax,dx,ay,dy);
    mg.output_z("z.0",ax,dx,ay,dy);
}
