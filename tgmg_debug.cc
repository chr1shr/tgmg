#include "tgmg.hh"

/** Prints the calculated matrix entries. This function is mainly used for
 * diagnostic purposes and debugging. */
template<class V,class M>
void tgmg_level<V,M>::print_rat() {
    int i,j,k;
    M* sp=s;
    for(j=0;j<n;j++) {
        for(i=0;i<m;i++,sp+=10) {
            printf("%d %d [",i,j);
            for(k=0;k<8;k++) printf("%g,",sp[k]);
            printf("%g] {%g}\n",sp[8],sp[9]);
        }
    }
}

/** Checks that the problem class is symmetric. */
template<class S,class V,class M>
void tgmg<S,V,M>::check_symmetric() {
    bool sym=true,ll,rr;
    int ri,li;
    for(int j=0,ij=0;j<n;j++) for(int i=0;i<m;i++,ij++) {

        // Check for which neighbors are available in the left and right
        // directions
        if(i>0) {li=-1;ll=true;}
        else if(x_prd) {li=m-1;ll=true;}
        else ll=false;
        if(i<m-1) {ri=1;rr=true;}
        else if(x_prd) {ri=1-m;rr=true;}
        else rr=false;

        // Check the four matrix entries to the up and to the right, so that
        // every pair of entries is considered exactly once
        if(rr) sym&=sym_check(i,ij+ri,q.a_cr(i,ij),q.a_cl(i+ri,ij+ri),"cl","cr");
        if(j<n-1||y_prd) {
            int uij=ij+(j<n-1?m:m-mn);
            if(ll) sym&=sym_check(ij,uij+li,q.a_ul(i,ij),q.a_dr(i+li,uij+li),"ul","dr");
            sym&=sym_check(ij,uij,q.a_uc(i,ij),q.a_dc(i,uij),"uc","dc");
            if(rr) sym&=sym_check(ij,uij+ri,q.a_ur(i,ij),q.a_dl(i+ri,uij+ri),"ur","dl");
        }
    }
    printf("The matrix is %ssymmetric\n",sym?"":"not ");
}

/** Checks to see if two matrix entries are consistent with being symmetric.
 * \param[in] ij the first grid point index.
 * \param[in] ij2 the second grid point index.
 * \param[in] v the value of the first matrix entry.
 * \param[in] v the value of the second matrix entry.
 * \param[in] p the two character code for the first matrix entry.
 * \param[in] p2 the two character code for the second matrix entry.
 * \return Whether the matrix entries are consistent with symmetry or not. */
template<class S,class V,class M>
bool tgmg<S,V,M>::sym_check(int ij,int ij2,M a,M a2,const char* p,const char* p2) {
    M da=a-a2;
    if(mod_sq(da)>q.acc) {
        printf("a_%s(%d,%d)=%.10g a_%s(%d,%d)=%.10g\n",p,ij%m,ij/m,a,p2,ij2%m,ij2/m,a2);
        return false;
    }
    return true;
}

/** Checks that the mul_a routine in the problem class is consistent with the
 * specification of individual matrix terms. This routine has to clear the
 * solution array. */
template<class S,class V,class M>
void tgmg<S,V,M>::check_consistent() {
    bool cons=true,ll,rr;
    int ri=0,li=0;
    clear_z();
    for(int j=0,ij=0;j<n;j++) for(int i=0;i<m;i++,ij++) {

        // Check for which neighbors are available in the left and right
        if(i>0) {li=-1;ll=true;}
        else if(x_prd) {li=m-1;ll=true;}
        else ll=false;
        if(i<m-1) {ri=1;rr=true;}
        else if(x_prd) {ri=1-m;rr=true;}
        else rr=false;

        // Test lower matrix entries
        if(j>0||y_prd) {
            int dij=ij-(j>0?m:m-mn);
            if(ll) cons&=cons_check(i,ij,dij+li,q.a_dl(i,ij),"dl");
            cons&=cons_check(i,ij,dij,q.a_dc(i,ij),"dc");
            if(rr) cons&=cons_check(i,ij,dij+ri,q.a_dr(i,ij),"dr");
        }

        // Test middle matrix entries
        if(ll) cons&=cons_check(i,ij,ij+li,q.a_cl(i,ij),"cl");
        if(rr) cons&=cons_check(i,ij,ij+ri,q.a_cr(i,ij),"cr");

        // Test upper matrix entries
        if(j<n-1||y_prd) {
            int uij=ij+(j<n-1?m:m-mn);
            if(ll) cons&=cons_check(i,ij,uij+li,q.a_ul(i,ij),"ul");
            cons&=cons_check(i,ij,uij,q.a_uc(i,ij),"uc");
            if(rr) cons&=cons_check(i,ij,uij+ri,q.a_ur(i,ij),"ur");
        }
    }
    printf("The mul_a routine is %sconsistent with the matrix terms\n",cons?"":"not ");
}

/** Checks to see if a matrix entry is consistent with the mul_a routine. It
 * assumes that the solution array is zero, and it leaves the array in this
 * condition upon completion.
 * \param[in] i the horizontal co-ordinate of the grid point to consider.
 * \param[in] ij the grid point index.
 * \param[in] ij2 the grid point referenced by the matrix entry.
 * \param[in] v the value of the matrix entry.
 * \param[in] p the two character code for the matrix entry.
 * \return Whether the matrix entry is consistent or not. */
template<class S,class V,class M>
bool tgmg<S,V,M>::cons_check(int i,int ij,int ij2,M a,const char* p) {
    z[ij2]=V(1.);V v=a*V(1.),v2=q.mul_a(i,ij),dv=v2-v;z[ij2]=0.;
    if(mod_sq(dv)>q.acc) {
        printf("a_%s at (%d,%d)=%.10g but mul_a gives %.10g\n",p,ij%m,ij/m,v,v2);
        return false;
    }
    return true;
}
