#include <cmath>

#include "poisson_fd.hh"
#include "tgmg.hh"

/** Creates the class representing a finite-difference discretization of the
 * Poisson equation on a rectangle for use with TGMG.
 * \param[in] (m_,n_) the number of grid points to use in the horizontal and
 *                    vertical directions.
 * \param[in] (ax_,bx_) the lower and upper x-coordinate rectangle bounds.
 * \param[in] (ay_,by_) the lower and upper y-coordinate rectangle bounds. */
poisson_fd::poisson_fd(const int m_,const int n_,const double ax_,const double bx_,const double ay_,const double by_)
        : m(m_), n(n_), mn(m_*n_), ax(ax_), ay(ay_), dx((bx_-ax)/(m-1)),
        dy((by_-ay)/(n-1)), hcc(-2/(dx*dx)-2/(dy*dy)), hcc_inv(1./hcc),
        hxx(1./(dx*dx)), hyy(1./(dy*dy)), acc(tgmg_accuracy(hcc,1e4)),
        z(new double[mn]), b(new double[mn]) {}

/** The class destructor frees the dynamically allocated memory. */
poisson_fd::~poisson_fd() {
    delete [] b;
    delete [] z;
}

/** Initializes the source term to be a Gaussian.
 * \param[in] (gx,gy) the center of the Gaussian.
 * \param[in] r the radius of the Gaussian.
 * \param[in] amp the amplitude of the Gaussian. */
void poisson_fd::gaussian_source_term(double gx,double gy,double r,double amp) {
    double lam=-0.5/(r*r);
#pragma omp parallel for
    for(int j=0;j<n;j++) {
        double y=ay-gy+dy*j,x=ax-gx;
        for(double *bp=b+j*m,*be=bp+m;bp<be;bp++,x+=dx)
            *bp=amp*exp(lam*(x*x+y*y));
    }
}

// Explicit instantiation
#include "tgmg.cc"
template class tgmg<poisson_fd,double,double>;
template class tgmg_base<tgmg_level<double,double>,double,double>;
template void tgmg_base<poisson_fd,double,double>::rat();
template void tgmg_base<poisson_fd,double,double>::clear_z();
template void tgmg_base<poisson_fd,double,double>::output(char const*,double*,double,double,double,double);
