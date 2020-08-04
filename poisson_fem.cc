#include "tgmg.hh"
#include "poisson_fem.hh"

/** An alternative constructor that independently sets up the multigrid
 * class for testing purposes.
 * \param[in] (m_,n_) the dimensions of the grid.
 * \param[in] (x_prd_,y_prd_) the periodicity of the grid.
 * \param[in] alloc_ whether to internally allocate the source and solution arrays.
 * \param[in] (dx,dy) the grid spacings in the x and y directions, respectively. */
poisson_fem::poisson_fem(int m_,int n_,bool x_prd_,bool y_prd_,bool alloc_,double dx,double dy) :
    m(m_), n(n_), mn(m*n), x_prd(x_prd_), y_prd(y_prd_), alloc(alloc_),
    dydx(dy/dx), dxdy(dx/dy), fm(4./3.*(dxdy+dydx)), fm_inv(1.0/fm),
    fey(fm_inv/3.*(-2*dxdy+dydx)), hey(0.5*fey), fex(fm_inv/3.*(-2*dydx+dxdy)),
    hex(0.5*fex), fc(-fm_inv/6.*(dxdy+dydx)), acc(tgmg_accuracy(1.,1e4)),
    z(alloc?new double[mn]:NULL), b(alloc?new double[mn]:NULL) {}

/** The class destructor frees the dynamically allocated memory. */
poisson_fem::~poisson_fem() {
    if(alloc) {
        delete [] b;
        delete [] z;
    }
}

/** Calculates the result of multiplying the matrix A by the solution vector at
 * one grid point. The diagonal term in the matrix A is omitted.
 * \param[in] i the horizontal index of the grid point.
 * \param[in] ij the overall index of the grid point. */
double poisson_fem::mul_a(int i,int ij) {
    double *w=z+ij;
    if(ij<m) {
        if(i==0) {
            return fc*w[m+1]+hey*w[m]+hex*w[1];
        } else if (i==m-1) {
            return fc*w[m-1]+hey*w[m]+hex*w[-1];
        } else return fc*(w[m-1]+w[m+1])+fey*w[m]+hex*(w[-1]+w[1]);
    } else if(ij>=mn-m) {
        if(i==0) {
            return fc*w[-m+1]+hey*w[-m]+hex*w[1];
        } else if(i==m-1) {
            return fc*w[-m-1]+hey*w[-m]+hex*w[-1];
        } else return fc*(w[-m-1]+w[-m+1])+fey*w[-m]+hex*(w[-1]+w[1]);
    } else {
        if(i==0) {
            return fc*(w[-m+1]+w[m+1])+hey*(w[-m]+w[m])+fex*w[1];
        } else if(i==m-1) {
            return fc*(w[-m-1]+w[m-1])+hey*(w[-m]+w[m])+fex*w[-1];
        } else return fc*(w[-m-1]+w[-m+1]+w[m-1]+w[m+1])
                      +fey*(w[-m]+w[m])+fex*(w[-1]+w[1]);
    }
}

// Explicit instantiation
#include "tgmg.cc"
#include "tgmg_debug.cc"
template class tgmg<poisson_fem,double,double>;
template class tgmg_base<tgmg_level<double,double>,double,double>;
template void tgmg_base<poisson_fem,double,double>::rat();
template void tgmg_base<poisson_fem,double,double>::jacobi();
template void tgmg_base<poisson_fem,double,double>::gauss_seidel();
template void tgmg_base<poisson_fem,double,double>::gauss_seidel_reverse();
template void tgmg_base<poisson_fem,double,double>::output(char const*,double*,double,double,double,double);
