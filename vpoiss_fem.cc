#include "tgmg.hh"
#include "vpoiss_fem.hh"

/** An alternative constructor that independently sets up the multigrid
 * class for testing purposes.
 * \param[in] (m_,n_) the dimensions of the grid.
 * \param[in] (x_prd_,y_prd_) the periodicity of the grid.
 * \param[in] alloc_ whether to internally allocate the source and solution
 *                   arrays.
 * \param[in] (dx,dy) the grid spacings in the x and y directions,
 *                    respectively. */
vpoiss_fem::vpoiss_fem(int m_,int n_,bool x_prd_,bool y_prd_,bool alloc_,double dx,double dy) :
    m(m_), n(n_), mn(m*n), x_prd(x_prd_), y_prd(y_prd_), alloc(alloc_),
    dydx(dy/dx), dxdy(dx/dy), fm(1./3.*(dxdy+dydx)), fm_inv(1.0/fm),
    fey(1./3.*(-2*dxdy+dydx)), hey(0.5*fey), fex(1./3.*(-2*dydx+dxdy)),
    hex(0.5*fex), fc(-1./6.*(dxdy+dydx)), acc(tgmg_accuracy(1.,1e4)),
    pfix(new bool[mn]), z(alloc?new double[mn]:NULL),
    b(alloc?new double[mn]:NULL), c(new double[mn]) {
    for(int i=0;i<mn;i++) c[i]=1.;
    for(int i=0;i<mn;i++) pfix[i]=false;
}

/** The class destructor frees the dynamically allocated memory. */
vpoiss_fem::~vpoiss_fem() {
    delete [] c;
    if(alloc) {
        delete [] b;
        delete [] z;
    }
    delete [] pfix;
}

/** Sets all the edges to have Dirichlet boundary conditions. */
void vpoiss_fem::set_dirichlet() {
    bool *pp=pfix;
    while(pp<pfix+m) *(pp++)=true;
    while(pp<pfix+mn-m) {*pp=pp[m-1]=true;pp+=m;}
    while(pp<pfix+mn) *(pp++)=true;
}

/** Sets a block of the prefactor array to a given value.
 * \param[in] (li,ui) the range of x indices to set.
 * \param[in] (lj,uj) the range of y indices to set.
 * \param[in] v the value to set. */
void vpoiss_fem::set_pre_block(int li,int ui,int lj,int uj,double v) {
    for(int j=lj;j<uj;j++) for(int i=li;i<ui;i++) c[j*m+i]=v;
}

double vpoiss_fem::mul_a(int i,int ij) {
    if(pfix[ij]) return 0;
    double *w=z+ij;
    const double cdl=c_dl(i,ij),cdr=c_dr(i,ij),
                 cul=c_ul(i,ij),cur=c_ur(i,ij);

    // Interior
    if(i>0&&i<m-1&&ij>=m&&ij<mn-m)
        return cul*(w[-1]*hex+w[+m-1]*fc+w[+m]*hey)
              +cdl*(w[-1]*hex+w[-m-1]*fc+w[-m]*hey)
              +cur*(w[+1]*hex+w[+m+1]*fc+w[+m]*hey)
              +cdr*(w[+1]*hex+w[-m+1]*fc+w[-m]*hey);

    // Compute constants and memory shifts for x-periodicity
    int sl=-1,sr=1;
    double lmu=1,rmu=1,ans,
           hex_l=hex*(cdl+cul),hex_r=hex*(cdr+cur),
           hey_d=hey*(cdl+cdr),hey_u=hey*(cul+cur);
    if(i==0) {
        sl+=m;
        if(!x_prd) {lmu=0,hey_d=hey*cdr,hey_u=hey*cur;}
    } else if(i==m-1) {
        sr-=m;
        if(!x_prd) {rmu=0;hey_d=hey*cdl,hey_u=hey*cul;}
    }

    // Assemble terms, taking into account y-periodicity
    if(ij>=m) {
        ans=fc*(lmu*w[-m+sl]*cdl+rmu*w[-m+sr]*cdr)+hey_d*w[-m];
    } else if(y_prd) {
        ans=fc*(lmu*w[mn-m+sl]*cdl+rmu*w[mn-m+sr]*cdr)+hey_d*w[mn-m];
    } else {
        ans=0,hex_l=hex*cul,hex_r=hex*cur;
    }

    if(ij<mn-m) {
        ans+=fc*(lmu*w[m+sl]*cul+rmu*w[m+sr]*cur)+hey_u*w[m];
    } else if(y_prd) {
        ans+=fc*(lmu*w[m-mn+sl]*cul+rmu*w[m-mn+sr]*cur)+hey_u*w[m-mn];
    } else {
        hex_l=hex*cdl,hex_r=hex*cdr;
    }
    return ans+hex_l*lmu*w[sl]+hex_r*rmu*w[sr];
}

// Explicit instantiation
#include "tgmg.cc"
#include "tgmg_debug.cc"
template class tgmg<vpoiss_fem,double,double>;
template class tgmg_base<tgmg_level<double,double>,double,double>;
template void tgmg_base<vpoiss_fem,double,double>::jacobi();
template void tgmg_base<vpoiss_fem,double,double>::gauss_seidel();
template void tgmg_base<vpoiss_fem,double,double>::gauss_seidel_reverse();
template void tgmg_base<vpoiss_fem,double,double>::output(char const*,double*,double,double,double,double);
