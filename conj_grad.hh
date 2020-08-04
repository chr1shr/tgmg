#ifndef CONJ_GRAD_HH
#define CONJ_GRAD_HH

#include <cstring>

class conj_grad {
    public:
        const int dof;
        double *x;
        double *b;
        conj_grad(int dof_);
        ~conj_grad();
        void solve_cg(bool verbose=false);
        void solve_pre_cg(bool verbose=false);
        virtual void mul_A(double *in,double *out) = 0;
        virtual void M_inv(double *in,double *out);
    protected:
        inline void daxpy(double alpha,double *in,double *out) {
#pragma omp parallel for
            for(int i=0;i<dof;i++) out[i]+=alpha*in[i];
        }
        inline void dscal(double alpha,double *x) {
#pragma omp parallel for
            for(double *xp=x;xp<x+dof;xp++) *xp*=alpha;
        }
        inline void copy(double *in,double *out) {
            memcpy(out,in,dof*sizeof(double));
        }
        inline double dot(double *x,double *y) {
            double s=*x*(*y);
#pragma omp parallel for reduction(+:s)
            for(int i=1;i<dof;i++) s+=x[i]*y[i];
            return s;
        }
    private:
        double *rk;
        double *pk;
        double *zk;
        double *yk;
};

#endif
