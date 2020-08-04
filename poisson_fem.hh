#ifndef POISSON_FEM_HH
#define POISSON_FEM_HH

#include "tgmg.hh"

struct poisson_fem {
    /** The number of grid points in the horizontal direction. */
    const int m;
    /** The number of grid points in the vertical direction. */
    const int n;
    /** The total number of grid cells. */
    const int mn;
    /** The periodicity in the x direction. */
    const bool x_prd;
    /** The periodicity in the y direction. */
    const bool y_prd;
    /** Whether the source and solution arrays were allocated internally. */
    const bool alloc;
    /** The mode to use for the Gauss-Seidel smoothing. (0=default) */
    static const char gs_mode=0;
    /** The ratio of the vertical to horizontal grid spacings. */
    const double dydx;
    /** The ratio of the horizontal to vertical grid spacings. */
    const double dxdy;
    /** The central term in the finite element stencil. */
    const double fm;
    /** The reciprocal of the central term. */
    const double fm_inv;
    /** The stencil term connecting to a vertical neighbor. */
    const double fey;
    /** Half the stencil term connecting to a vertical neighbor. */
    const double hey;
    /** The stencil term connecting to a horizontal neighbor. */
    const double fex;
    /** Half the stencial term connection to a horizontal neighbor. */
    const double hex;
    /** The stencil term connecting to a diagonal neighbor. */
    const double fc;
    /** Threshold on L_2 norm of residual to terminate the multigrid solve. */
    const double acc;
    /** The solution vector array. */
    double* z;
    /** The source term array. */
    double* b;
    poisson_fem(int m_,int n_,bool x_prd,bool y_prd,bool alloc_,double dx,double dy);
    ~poisson_fem();
    inline bool not_l(int i) {return x_prd||i>0;}
    inline bool not_r(int i) {return x_prd||i<m-1;}
    inline bool not_lr(int i) {return x_prd||(i>0&&i<m-1);}
    inline bool not_d(int ij) {return y_prd||ij>=m;}
    inline bool not_u(int ij) {return y_prd||ij<mn-m;}
    inline bool not_du(int ij) {return y_prd||(ij>=m&&ij<mn-m);}
    inline double a_dl(int i,int ij) {
        return not_d(ij)&&not_l(i)?fc:0;
    }
    inline double a_dr(int i,int ij) {
        return not_d(ij)&&not_r(i)?fc:0;
    }
    inline double a_ul(int i,int ij) {
        return not_u(ij)&&not_l(i)?fc:0;
    }
    inline double a_ur(int i,int ij) {
        return not_u(ij)&&not_r(i)?fc:0;
    }
    inline double a_dc(int i,int ij) {
        return not_d(ij)?(not_lr(i)?fey:hey):0;
    }
    inline double a_uc(int i,int ij) {
        return not_u(ij)?(not_lr(i)?fey:hey):0;
    }
    inline double a_cl(int i,int ij) {
        return not_l(i)?(not_du(ij)?fex:hex):0;
    }
    inline double a_cr(int i,int ij) {
        return not_r(i)?(not_du(ij)?fex:hex):0;
    }
    inline double a_cc(int i,int ij) {
        return (not_lr(i)?1:0.5)*(not_du(ij)?1:0.5);
    }
    inline double inv_cc(int i,int ij,double v) {
        return (not_lr(i)?1:2)*(not_du(ij)?1:2)*v;
    }
    double mul_a(int i,int ij);
};

#endif
