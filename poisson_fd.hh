#ifndef POISSON_FD_HH
#define POISSON_FD_HH

#include "tgmg.hh"

// Multisetup structure for Poisson problem
struct poisson_fd {
    /** Grid dimensions. */
    const int m;
    const int n;
    /** Total number of gridpoints. */
    const int mn;
    /** Periodicity in the x and y directions. */
    static const bool x_prd=false;
    static const bool y_prd=false;
    /** The mode to use for the Gauss-Seidel smoothing. (0=default) */
    static const char gs_mode=0;
    /** The minimum x coordinate of the grid. */
    const double ax;
    /** The minimum y coordinate of the grid. */
    const double ay;
    /** Grid spacings in the x and y directions. */
    const double dx,dy;
    /** Stencil entries. */
    const double hcc,hcc_inv,hxx,hyy;
    /** Threshold on L_2 norm of residual to terminate the multigrid solve. */
    const double acc;
    /** The solution vector array. */
    double* const z;
    /** The source term array. */
    double* const b;
    poisson_fd(const int m_,const int n_,const double ax_,const double bx_,const double ay_,const double by_);
    ~poisson_fd();
    void gaussian_source_term(double gx,double dy,double lam,double amp);
    /** Function to determine whether a grid point is on the edge or not.
     */
    inline bool edge(int i,int ij) {return i==0||i==m-1||ij>=mn-m||ij<m;}
    /** Functions to specify the corner stencil entries. */
    inline double a_dl(int i,int ij) {return 0;}
    inline double a_dr(int i,int ij) {return 0;}
    inline double a_ul(int i,int ij) {return 0;}
    inline double a_ur(int i,int ij) {return 0;}
    /** Functions to specify the vertical stencil entries. */
    inline double a_dc(int i,int ij) {return edge(i,ij)?0:hyy;}
    inline double a_uc(int i,int ij) {return edge(i,ij)?0:hyy;}
    /** Functions to specify the horizontal stencil entries. */
    inline double a_cl(int i,int ij) {return edge(i,ij)?0:hxx;}
    inline double a_cr(int i,int ij) {return edge(i,ij)?0:hxx;}
    /** Function to specify the central stencil entry (on the diagonal of
     * the linear system). */
    inline double a_cc(int i,int ij) {return hcc;}
    /** Function to multiply by the reciprocal of the central stencil
     * entry. This is specified as a separate function for computational
     * efficiency. */
    inline double inv_cc(int i,int ij,double v) {return hcc_inv*v;}
    /** Calculates the ith component of the multiplication (A-D)z, needed
     * in the Gauss--Seidel smoothing iteration. */
    inline double mul_a(int i,int ij) {
        return edge(i,ij)?0:hxx*(z[ij+1]+z[ij-1])+hyy*(z[ij+m]+z[ij-m]);
    }
};

#endif
