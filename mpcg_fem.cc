#include "mpcg_fem.hh"
#include <cstdlib>

/** An alternative constructor that independently sets up the multigrid
 * class for testing purposes.
 * \param[in] (m_,n_) the dimensions of the grid.
 * \param[in] (x_prd_,y_prd_) the periodicity of the grid.
 * \param[in] (dx,dy) the grid spacings in the x and y directions, respectively. */
mpcg_fem::mpcg_fem(int m_,int n_,bool x_prd_,bool y_prd_,double ax_,double dx_,double ay_,double dy_) :
    conj_grad(m_*n_), m(m_), n(n_), mn(m_*n_), x_prd(x_prd_), y_prd(y_prd_),
    ax(ax_), ay(ay_), dx(dx_), dy(dy_), pf(m_,n_,x_prd_,y_prd_,false,dx,dy),
    mg(pf,NULL,NULL) {
    pf.set_pre_block(0,m,0,30,1);
    pf.z=mg.z=mg.mg[0]->y=new double[mn];
    pf.b=mg.b=new double[mn];
    mg.setup();
    delete [] pf.z;
    delete [] pf.b;
}

void mpcg_fem::mul_A(double *in,double *out) {
    pf.z=in;
#pragma omp parallel for
    for(int j=0;j<mn;j+=m) {
        for(int ij=j,i=0;i<m;i++,ij++)
            out[ij]=pf.mul_a(i,ij)+pf.a_cc(i,ij)*in[ij];
    }
}

/** Adds a Gaussian to the source term
 * \param[in] (gx,gy) the center of the Gaussian.
 * \param[in] r the radius of the Gaussian.
 * \param[in] amp the amplitude of the Gaussian. */
void mpcg_fem::add_gaussian_source_term(double gx,double gy,double r,double amp) {
    double lam=-0.5/(r*r);
#pragma omp parallel for
    for(int j=0;j<n;j++) {
        double y=ay-gy+dy*j,x=ax-gx;
        for(double *bp=b+j*m,*be=bp+m;bp<be;bp++,x+=dx)
            *bp+=amp*exp(lam*(x*x+y*y));
    }
}

void mpcg_fem::M_inv(double *in,double *out) {
    for(int k=0;k<mn;k++) out[k]=0;
    pf.z=mg.z=mg.mg[0]->y=out;
    pf.b=mg.b=in;
    mg.v_cycle(1,1,10,true);
}
