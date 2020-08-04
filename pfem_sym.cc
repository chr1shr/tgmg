#include <cmath>
#include <algorithm>
#include <cstring>
#include <limits>

#include "tgmg.hh"
#include "vpoiss_fem.hh"

/** Checks if two numbers are the same within a tolerance.
 * \param[in] (a,b) the two numbers.
 * \return True if they are the same, false if not. */
bool same(double a,double b) {
    const double eps=std::numeric_limits<double>::epsilon()*10;
    return fabs(b-a)<=std::max(fabs(a),fabs(b))*eps;
}

// Grid dimensions
const int m=21,n=21;

// Total number of gridpoints
const int mn=m*n;

// Physical dimensions of the grid
const double ax=-3,bx=3,ay=ax,by=bx;

// Grid spacings
const double dx=(bx-ax)/(m-1),dy=(by-ay)/(n-1);

int main() {
    int j,ij;
    vpoiss_fem pfem(m,n,false,false,true,dx,dy);
    double *b=pfem.b,*z=pfem.z;
    tgmg<vpoiss_fem,double,double> mg(pfem,b,z);
    mg.verbose=2;

    // Set up the multigrid hierarchy
    mg.setup();

    // Create a large array for storing all of the Gauss--Seidel results
    double *w=new double[mn*mn];

    // Loop over all gridpoints
    for(ij=0;ij<mn;ij++) {

        // Clear solution. Clear source term, except for a single entry set to unity.
        for(j=0;j<mn;j++) z[j]=0.;
        for(j=0;j<mn;j++) b[j]=0.;
        b[ij]=1.;

        // Perform symmetric operation
        mg.v_cycle(1,1,10,true);
//        mg.gauss_seidel();
//        mg.gauss_seidel_reverse();

        // Store solution
        memcpy(w+ij*mn,z,mn*sizeof(double));
    }

    // Check whether the operation is indeed symmetric
    for(ij=0;ij<mn;ij++) for(int ij2=ij+1;ij2<mn;ij2++) {
        if(!same(w[ij2+ij*mn],w[ij+ij2*mn])) {
            printf("%d %d %d %d %g %g\n",ij%m,ij/m,ij2%m,ij2/m,w[ij2+ij*mn],w[ij+ij2*mn]);
        }
    }

    // Free the dynamically allocated memory
    delete [] w;
}
