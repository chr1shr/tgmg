#ifndef VPOISS_FEM_HH
#define VPOISS_FEM_HH

#include "tgmg.hh"

struct vpoiss_fem {
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
    /** An array for enforcing Dirichlet boundary conditions in the solution.
     */
    bool* const pfix;
    /** The solution vector array. */
    double* z;
    /** The source term array. */
    double* b;
    /** The varying field. */
    double* const c;
    vpoiss_fem(int m_,int n_,bool x_prd_,bool y_prd_,bool alloc_,double dx,double dy);
    ~vpoiss_fem();
    void set_dirichlet();
    void set_pre_block(int li,int ui,int lj,int uj,double v);
    inline bool not_l(int i) {return x_prd||i>0;}
    inline bool not_r(int i) {return x_prd||i<m-1;}
    inline bool not_lr(int i) {return x_prd||(i>0&&i<m-1);}
    inline bool not_d(int ij) {return y_prd||ij>=m;}
    inline bool not_u(int ij) {return y_prd||ij<mn-m;}
    inline bool not_du(int ij) {return y_prd||(ij>=m&&ij<mn-m);}
    inline double c_dl(int i,int ij) {return c[ij];}
    inline double c_dr(int i,int ij) {return not_r(i)?c[ij+1]:c[ij+1-m];}
    inline double c_ul(int i,int ij) {return not_u(ij)?c[ij+m]:c[ij+m-mn];}
    inline double c_ur(int i,int ij) {return c[ij+(not_r(i)?1:1-m)+(not_u(ij)?m:m-mn)];}
    inline double a_dl(int i,int ij) {return not_d(ij)&&not_l(i)?fc*c_dl(i,ij):0;}
    inline double a_dr(int i,int ij) {return not_d(ij)&&not_r(i)?fc*c_dr(i,ij):0;}
    inline double a_ul(int i,int ij) {return not_u(ij)&&not_l(i)?fc*c_ul(i,ij):0;}
    inline double a_ur(int i,int ij) {return not_u(ij)&&not_r(i)?fc*c_ur(i,ij):0;}
    inline double a_dc(int i,int ij) {
        return not_d(ij)?((not_l(i)?c_dl(i,ij):0)+(not_r(i)?c_dr(i,ij):0))*hey:0;
    }
    inline double a_uc(int i,int ij) {
        return not_u(ij)?((not_l(i)?c_ul(i,ij):0)+(not_r(i)?c_ur(i,ij):0))*hey:0;
    }
    inline double a_cl(int i,int ij) {
        return not_l(i)?((not_d(ij)?c_dl(i,ij):0)+(not_u(ij)?c_ul(i,ij):0))*hex:0;
    }
    inline double a_cr(int i,int ij) {
        return not_r(i)?((not_d(ij)?c_dr(i,ij):0)+(not_u(ij)?c_ur(i,ij):0))*hex:0;
    }
    inline double a_cc(int i,int ij) {
        return (not_d(ij)?(not_u(ij)?(not_l(i)?(c_ul(i,ij)+c_dl(i,ij)):0)
                                     +(not_r(i)?(c_ur(i,ij)+c_dr(i,ij)):0)
                                    :(not_l(i)?c_dl(i,ij):0)+(not_r(i)?c_dr(i,ij):0))
                         :(not_l(i)?c_ul(i,ij):0)+(not_r(i)?c_ur(i,ij):0))*fm;
    }
    inline double inv_cc(int i,int ij,double v) {return v/a_cc(i,ij);}
    double mul_a(int i,int ij);
};

#endif
