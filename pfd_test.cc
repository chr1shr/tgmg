#include <cmath>

#include "tgmg.hh"
#include "poisson_fd.hh"

int main() {

    // Grid dimensions
    const int m=1025,n=1025;

    // Physical dimensions of the rectangle to solve on
    const double ax=-1,bx=1,ay=ax,by=bx;

    // Create the problem class
    poisson_fd pfd(m,n,ax,bx,ay,by);

    // Create the multigrid solver
    tgmg<poisson_fd,double,double> mg(pfd,pfd.b,pfd.z);

    // Set the verbosity to maximum, and set up the multigrid hierarchy
    mg.verbose=3;
    mg.setup();

    // Set up the solution and source arrays
    mg.clear_z();
    pfd.gaussian_source_term(0.3,0.7,0.25,1);

    // Solve using multigrid V-cycles
    mg.solve_v_cycle();

    // Output the solutions in a format that can be read by Gnuplot using
    // either of the following two commands:
    //
    // 1. For color map plot
    // plot [-1:1] [-1:1] 'filename' matrix binary with image
    //
    // 2. For 3D plot
    // splot [-1:1] [-1:1] 'filename' matrix binary
    const double dx=(bx-ax)/(m-1),dy=(by-ay)/(n-1);
    mg.output_b("src.0",ax,dx,ay,dy);
    mg.output_z("sol.0",ax,dx,ay,dy);
}
