#include <cmath>
#include <ctime>

#include "tgmg.hh"
#include "poisson_fd.hh"

#ifdef _OPENMP
#include "omp.h"
inline double wtime() {return omp_get_wtime();}
#else
inline double wtime() {return 0;
}
#endif

int main() {

    // Grid dimensions and domain size
    const int m=1025,n=1025;
    const double ax=-8,bx=8,ay=ax,by=bx;

    // Total number of trials
    const int num_trials=500;

    // Total number of complete spins that the Gaussian source term is rotated
    // through. From this, define the angular step size in radians.
    const double spins=1;
    const double h=2.*M_PI*spins/num_trials;

    // Initialize the problem class, the multigrid solver, and set up the
    // multigrid hierarchy
    poisson_fd pfd(m,n,ax,bx,ay,by);
    tgmg<poisson_fd,double,double> mg(pfd,pfd.b,pfd.z);
    mg.setup();

    // Create class for predicting the number of V-cycles to perform before
    // testing for convergence
    tgmg_predict tp;

    // Initialize timers. Use the OpenMP routine for timing wall clock time,
    // and the clock() function for timing processor time
    double t0=wtime(),t1,t_setup=0,t_solve=0,dc_setup,dc_solve;
    clock_t c0=clock(),c1,c_setup=0,c_solve=0;

    // Do a number of trials where the source term for the solve is slowly
    // moved around the domain
    for(int k=0;k<num_trials;k++) {

        // Set up the source term and time it
        pfd.gaussian_source_term(6*cos(h*k),6*sin(h*k),1,1);
        t1=wtime();c1=clock();
        t_setup+=t1-t0;c_setup+=c1-c0;

        // Perform the multigrid solve and time it
        if(k%100==0) {mg.verbose=2;printf("Trial %d\n",k);}
        mg.solve_v_cycle(tp);
        mg.verbose=0;
        t0=wtime();c0=clock();
        t_solve+=t0-t1;c_solve+=c0-c1;
    }

    // Print the timing results and V-cycle statistics
    const double fac=1e3/num_trials,fac2=fac/CLOCKS_PER_SEC;
    t_setup*=fac;t_solve*=fac;
    dc_setup=double(c_setup)*fac2;
    dc_solve=double(c_solve)*fac2;
    int vc=tp.vcount;
    double vc_p_s=double(vc)/num_trials;
    printf("Source term setup:\n"
           "  Wall clock time : %.6g ms\n"
           "  Processor time  : %.6g ms [%.2f%%]\n\n"
           "Multigrid solve:\n"
           "  Wall clock time : %.6g ms\n"
           "  Processor time  : %.6g ms [%.2f%%]\n\n"
           "V-cycles:\n"
           "  Count: %d total, %g per solve\n"
           "  WC time per V-cycle    : %.6g ms\n"
           "  Proc. time per V-cycle : %.6g ms\n",
           t_setup,dc_setup,100.*dc_setup/t_setup,
           t_solve,dc_solve,100.*dc_solve/t_solve,
           vc,vc_p_s,t_solve/vc_p_s,dc_solve/vc_p_s);
}
