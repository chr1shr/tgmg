#ifndef MPCG_FEM
#define MPCG_FEM

#include "conj_grad.hh"
#include "poisson_fem.hh"
#include "vpoiss_fem.hh"
#include "tgmg.hh"

class mpcg_fem : public conj_grad {
    public:
        /** The horizontal grid size. */
        const int m;
        /** The vertical grid size. */
        const int n;
        /** The total number of grid points. */
        const int mn;
        /** The periodicity in the x direction. */
        const bool x_prd;
        /** The periodicity in the y direction. */
        const bool y_prd;
        /** The minimum x coordinate of the grid. */
        const double ax;
        /** The minimum y coordinate of the grid. */
        const double ay;
        /** Grid spacings in the x and y directions. */
        const double dx,dy;
        mpcg_fem(int m_,int n_,bool x_prd_,bool y_prd_,double ax,double dx,double ay,double dy);
        ~mpcg_fem() {}
        virtual void mul_A(double *in,double *out);
        virtual void M_inv(double *in,double *out);
        void add_gaussian_source_term(double gx,double gy,double r,double amp);
        inline void clear_b() {
            for(int ij=0;ij<mn;ij++) b[ij]=0;
        }
        inline void clear_x() {
            for(int ij=0;ij<mn;ij++) x[ij]=0;
        }
        inline void save_fields() {
            mg.b=b;
            mg.z=x;
            mg.output_b("b.0",ax,dx,ay,dy);
            mg.output_z("z.0",ax,dx,ay,dy);
        }
    private:
        vpoiss_fem pf;
        tgmg<vpoiss_fem,double,double> mg;
};

#endif
