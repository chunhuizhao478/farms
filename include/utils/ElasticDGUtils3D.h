#pragma once
#include "Moose.h"
#include "libmesh/dense_matrix.h"

/** 
 * Utilities for 3D elastodynamic discontinuous Galerkin interface computations
 * 
 * Implements coordinate transformations and material property calculations required for 
 * exact Godunov flux in discontinuous Galerkin elastic wave propagation.
 *
 * Key components:
 * - Orthonormal basis construction: Build local (n,t1,t2) from interface normal
 * - Rotation matrices: Transform between global (x,y,z) and local coordinates  
 * - Vector transformation: v_local = T^T * v_global, v_global = T * v_local
 * - Tensor transformation: σ_local = T^T * σ_global * T, σ_global = T * σ_local * T^T
 * - Acoustic impedances: Z_p = ρc_p, Z_s = ρc_s for P and S wave modes
 *
 * Based on standard elastic wave DG formulations with references:
 * [1] Stress tensor rotations: https://web.mit.edu/course/3/3.11/www/modules/trans.pdf
 * [2] Godunov elastic DG: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011JB008857
 * [3] 3D implementation: Pelties et al. (2012) GJI, doi:10.1111/j.1365-246X.2012.05626.x
 */
namespace ElasticDG3D
{
  using DM = DenseMatrix<Real>;

  // Build orthonormal basis (n,t1,t2) from interface normal vector
  // Uses Gram-Schmidt-like process to ensure numerical stability
  inline void orthonormal_basis(const RealVectorValue & n_in, RealVectorValue & n,
                                RealVectorValue & t1, RealVectorValue & t2)
  {
    n = n_in; n /= n.norm();
    // pick a reference not parallel to n
    RealVectorValue a = std::abs(n(0)) < 0.9 ? RealVectorValue(1,0,0) : RealVectorValue(0,1,0);
    t1 = n.cross(a); t1 /= t1.norm();
    t2 = n.cross(t1); t2 /= t2.norm();
  }

  // Build rotation matrix T with columns [n, t1, t2] and its transpose TT  
  inline void build_T(const RealVectorValue & n, const RealVectorValue & t1, const RealVectorValue & t2,
                      DM & T, DM & TT)
  {
    T.resize(3,3); TT.resize(3,3);
    // columns are (n, t1, t2)
    T(0,0)=n(0);  T(1,0)=n(1);  T(2,0)=n(2);
    T(0,1)=t1(0); T(1,1)=t1(1); T(2,1)=t1(2);
    T(0,2)=t2(0); T(1,2)=t2(1); T(2,2)=t2(2);
    // TT = T^T
    for (unsigned i=0;i<3;i++) for (unsigned j=0;j<3;j++) TT(i,j)=T(j,i);
  }

  // Transform vector from global to local coordinates: v_local = T^T * v_global
  inline void vec_to_local(const DM & TT, const RealVectorValue & v, Real & vn, Real & vt1, Real & vt2)
  {
    vn  = TT(0,0)*v(0) + TT(0,1)*v(1) + TT(0,2)*v(2);
    vt1 = TT(1,0)*v(0) + TT(1,1)*v(1) + TT(1,2)*v(2);
    vt2 = TT(2,0)*v(0) + TT(2,1)*v(1) + TT(2,2)*v(2);
  }

  // Transform vector from local to global coordinates: v_global = T * v_local
  inline void vec_to_global(const DM & T, Real vn, Real vt1, Real vt2, RealVectorValue & v)
  {
    v(0) = T(0,0)*vn + T(0,1)*vt1 + T(0,2)*vt2;
    v(1) = T(1,0)*vn + T(1,1)*vt1 + T(1,2)*vt2;
    v(2) = T(2,0)*vn + T(2,1)*vt1 + T(2,2)*vt2;
  }

  // Transform stress tensor from global to local coordinates: σ_local = T^T * σ_global * T
  inline void ten_to_local(const DM & T, const DM & TT,
                           Real sxx, Real sxy, Real sxz,
                           Real syy, Real syz, Real szz,
                           // outputs (s_nn, s_nt1, s_nt2, s_t1t1, s_t1t2, s_t2t2)
                           Real & snn, Real & snt1, Real & snt2,
                           Real & st1t1, Real & st1t2, Real & st2t2)
  {
    DM Sg(3,3), tmp(3,3), Sl(3,3);
    Sg(0,0)=sxx; Sg(0,1)=sxy; Sg(0,2)=sxz;
    Sg(1,0)=sxy; Sg(1,1)=syy; Sg(1,2)=syz;
    Sg(2,0)=sxz; Sg(2,1)=syz; Sg(2,2)=szz;
    // Sl = T^T Sg T = TT * Sg * T
    // tmp = TT * Sg
    for (unsigned i=0;i<3;i++) for (unsigned j=0;j<3;j++)
      tmp(i,j) = TT(i,0)*Sg(0,j) + TT(i,1)*Sg(1,j) + TT(i,2)*Sg(2,j);
    // Sl = tmp * T
    for (unsigned i=0;i<3;i++) for (unsigned j=0;j<3;j++)
      Sl(i,j) = tmp(i,0)*T(0,j) + tmp(i,1)*T(1,j) + tmp(i,2)*T(2,j);
    snn   = Sl(0,0); snt1  = Sl(0,1); snt2  = Sl(0,2);
    st1t1 = Sl(1,1); st1t2 = Sl(1,2); st2t2 = Sl(2,2);
  }

  // Transform stress tensor from local to global coordinates: σ_global = T * σ_local * T^T
  inline void ten_to_global(const DM & T, const DM & TT,
                            Real snn, Real snt1, Real snt2,
                            Real st1t1, Real st1t2, Real st2t2,
                            Real & sxx, Real & sxy, Real & sxz,
                            Real & syy, Real & syz, Real & szz)
  {
    DM Sl(3,3), tmp(3,3), Sg(3,3);
    Sl(0,0)=snn; Sl(0,1)=snt1; Sl(0,2)=snt2;
    Sl(1,0)=snt1; Sl(1,1)=st1t1; Sl(1,2)=st1t2;
    Sl(2,0)=snt2; Sl(2,1)=st1t2; Sl(2,2)=st2t2;
    // tmp = T * Sl
    for (unsigned i=0;i<3;i++) for (unsigned j=0;j<3;j++)
      tmp(i,j) = T(i,0)*Sl(0,j) + T(i,1)*Sl(1,j) + T(i,2)*Sl(2,j);
    // Sg = tmp * T^T
    for (unsigned i=0;i<3;i++) for (unsigned j=0;j<3;j++)
      Sg(i,j) = tmp(i,0)*TT(j,0) + tmp(i,1)*TT(j,1) + tmp(i,2)*TT(j,2);
    sxx = Sg(0,0); sxy = Sg(0,1); sxz = Sg(0,2);
    syy = Sg(1,1); syz = Sg(1,2); szz = Sg(2,2);
  }

  inline void impedances(Real lambda, Real mu, Real rho, Real & cp, Real & cs, Real & Zp, Real & Zs)
  {
    cp = std::sqrt((lambda + 2.0*mu)/rho);  // P speed
    cs = std::sqrt(mu/rho);                 // S speed
    Zp = rho * cp; Zs = rho * cs;           // impedances
  }
}