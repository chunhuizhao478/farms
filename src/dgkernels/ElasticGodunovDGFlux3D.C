#include "ElasticGodunovDGFlux3D.h"
using namespace ElasticDG3D;

registerMooseObject("farmsApp", ElasticGodunovDGFlux3D);

InputParameters ElasticGodunovDGFlux3D::validParams()
{
  InputParameters params = DGKernel::validParams();
  // velocities
  params.addRequiredCoupledVar("ux","x-velocity"); params.addRequiredCoupledVar("uy","y-velocity");
  params.addRequiredCoupledVar("uz","z-velocity");
  // stresses
  params.addRequiredCoupledVar("sxx",""); params.addRequiredCoupledVar("syy","");
  params.addRequiredCoupledVar("szz",""); params.addRequiredCoupledVar("sxy","");
  params.addRequiredCoupledVar("sxz",""); params.addRequiredCoupledVar("syz","");
  // materials
  params.addRequiredParam<MaterialPropertyName>("lambda_name","λ");
  params.addRequiredParam<MaterialPropertyName>("mu_name","μ");
  params.addRequiredParam<MaterialPropertyName>("rho_name","ρ");
  params.addParam<std::vector<SubdomainName>>("skip_block_pair",{},"Exactly two block names: if face is between them skip this flux (fault handled elsewhere)");
  return params;
}

ElasticGodunovDGFlux3D::ElasticGodunovDGFlux3D(const InputParameters & params)
  : DGKernel(params),
    _ux(coupledValue("ux")), _ux_neighbor(coupledNeighborValue("ux")),
    _uy(coupledValue("uy")), _uy_neighbor(coupledNeighborValue("uy")),
    _uz(coupledValue("uz")), _uz_neighbor(coupledNeighborValue("uz")),
    _sxx(coupledValue("sxx")), _sxx_neighbor(coupledNeighborValue("sxx")),
    _syy(coupledValue("syy")), _syy_neighbor(coupledNeighborValue("syy")),
    _szz(coupledValue("szz")), _szz_neighbor(coupledNeighborValue("szz")),
    _sxy(coupledValue("sxy")), _sxy_neighbor(coupledNeighborValue("sxy")),
    _sxz(coupledValue("sxz")), _sxz_neighbor(coupledNeighborValue("sxz")),
    _syz(coupledValue("syz")), _syz_neighbor(coupledNeighborValue("syz")),
    _lambda(getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("lambda_name"))),
    _lambda_n(getNeighborMaterialProperty<Real>(getParam<MaterialPropertyName>("lambda_name"))),
    _mu(getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("mu_name"))),
    _mu_n(getNeighborMaterialProperty<Real>(getParam<MaterialPropertyName>("mu_name"))),
    _rho(getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("rho_name"))),
    _rho_n(getNeighborMaterialProperty<Real>(getParam<MaterialPropertyName>("rho_name")))
{
  if (isParamValid("skip_block_pair"))
  {
    const auto & names = getParam<std::vector<SubdomainName>>("skip_block_pair");
    if (names.size()==2)
    { _use_skip_pair = true; _skip_a = _mesh.getSubdomainID(names[0]); _skip_b = _mesh.getSubdomainID(names[1]); }
  }
}

Real ElasticGodunovDGFlux3D::computeQpResidual(Moose::DGResidualType type)
{
  // Skip faces if between specified block pair (fault handled by separate kernel)
  if (_use_skip_pair)
  {
    const Elem * e = _current_elem;      // minus side
    const Elem * n = _neighbor_elem;     // plus side (provided by DGKernel)
    if (e && n)
    {
      const SubdomainID a = e->subdomain_id();
      const SubdomainID b = n->subdomain_id();
      if ((a == _skip_a && b == _skip_b) || (a == _skip_b && b == _skip_a))
        return 0.0; // skip flux, handled by specialized fault kernel
    }
  }

  // Implementation follows Pelties et al. (2012) "Three-dimensional dynamic rupture simulation 
  // with a high-order discontinuous Galerkin method on unstructured tetrahedral meshes"
  // Geophysical Journal International, doi:10.1111/j.1365-246X.2012.05626.x
  // Specifically implements the Godunov flux from equation (13) and surrounding formulation

  // 1) Build local coordinate system (n,t1,t2) and rotation matrices
  // Following Section 2.2: Local coordinate transformation for interface flux computation
  RealVectorValue n_glob = _normals[_qp], n,t1,t2;
  orthonormal_basis(n_glob, n, t1, t2);
  DenseMatrix<Real> T, TT; build_T(n,t1,t2,T,TT);

  // 2) Compute acoustic impedances for both sides of the interface
  // Each element uses its own material properties as per standard DG formulation
  // Left side (element): Z_p^-, Z_s^- 
  Real cp_m, cs_m, Zp_m, Zs_m; 
  impedances(_lambda[_qp], _mu[_qp], _rho[_qp], cp_m, cs_m, Zp_m, Zs_m);
  
  // Right side (neighbor): Z_p^+, Z_s^+
  Real cp_p, cs_p, Zp_p, Zs_p; 
  impedances(_lambda_n[_qp], _mu_n[_qp], _rho_n[_qp], cp_p, cs_p, Zp_p, Zs_p);

  // 3) Transform stress tensors and velocity vectors to local coordinate system
  // Transform from global (x,y,z) to local (n,t1,t2) coordinates
  // This enables application of 1D Riemann solver in normal direction
  Real snn_m,snt1_m,snt2_m, st11_m,st12_m,st22_m;
  Real snn_p,snt1_p,snt2_p, st11_p,st12_p,st22_p;
  ten_to_local(T,TT,_sxx[_qp],_sxy[_qp],_sxz[_qp],_syy[_qp],_syz[_qp],_szz[_qp],
               snn_m,snt1_m,snt2_m,st11_m,st12_m,st22_m);
  ten_to_local(T,TT,_sxx_neighbor[_qp],_sxy_neighbor[_qp],_sxz_neighbor[_qp],
               _syy_neighbor[_qp],_syz_neighbor[_qp],_szz_neighbor[_qp],
               snn_p,snt1_p,snt2_p,st11_p,st12_p,st22_p);

  Real un_m,ut1_m,ut2_m, un_p,ut1_p,ut2_p;
  vec_to_local(TT, RealVectorValue(_ux[_qp],_uy[_qp],_uz[_qp]), un_m, ut1_m, ut2_m);
  vec_to_local(TT, RealVectorValue(_ux_neighbor[_qp],_uy_neighbor[_qp],_uz_neighbor[_qp]),
               un_p, ut1_p, ut2_p);

  // 4) Godunov (Riemann) states in local coordinates  (Eq. 13)
  // Paper notation: (+) = plus side, (-) = minus side. Here we map: (-)= element (m), (+)= neighbor (p).
  // Eq. (13) lines reproduced (local frame where first index is normal direction):
  //   2 σ_nn^G = (σ_nn^+ + σ_nn^-) + ρ c_p (u_n^- - u_n^+)        (13a)
  //   2 σ_nt1^G = (σ_nt1^+ + σ_nt1^-) + (μ / c_s)(u_t1^- - u_t1^+) (13b)
  //   2 σ_nt2^G = (σ_nt2^+ + σ_nt2^-) + (μ / c_s)(u_t2^- - u_t2^+) (13c)
  //   2 u_n^G   = (u_n^+ + u_n^-) + (1/(ρ c_p))(σ_nn^- - σ_nn^+)   (13d)
  //   2 u_t1^G  = (u_t1^+ + u_t1^-) + (c_s/μ)(σ_nt1^- - σ_nt1^+)   (13e)
  //   2 u_t2^G  = (u_t2^+ + u_t2^-) + (c_s/μ)(σ_nt2^- - σ_nt2^+)   (13f)
  // Assumption: material is locally homogeneous across the face (as in the paper). We therefore
  // use minus-side (element) properties (ρ, c_p, μ, c_s). For heterogeneous interfaces a symmetric
  // impedance form should replace these expressions (not implemented here).

  // (13d) & (13a): P-wave (normal) Godunov states
  const Real unG   = 0.5 * ( (un_p + un_m) + (1.0/(cp_m * _rho[_qp])) * (snn_m - snn_p) );
  const Real snnG  = 0.5 * ( (snn_p + snn_m) + (cp_m * _rho[_qp])     * (un_m - un_p) );

  // (13e,13b) S-wave (t1) Godunov states
  const Real ut1G  = 0.5 * ( (ut1_p + ut1_m) + (cs_m / _mu[_qp]) * (snt1_m - snt1_p) );
  const Real snt1G = 0.5 * ( (snt1_p + snt1_m) + (_mu[_qp] / cs_m) * (ut1_m - ut1_p) );

  // (13f,13c) S-wave (t2) Godunov states
  const Real ut2G  = 0.5 * ( (ut2_p + ut2_m) + (cs_m / _mu[_qp]) * (snt2_m - snt2_p) );
  const Real snt2G = 0.5 * ( (snt2_p + snt2_m) + (_mu[_qp] / cs_m) * (ut2_m - ut2_p) );

  // 5) Physical fluxes evaluated at Godunov state (after Eq. 13)
  // Governing first-order system in local frame (normal derivative part only):
  //   ∂σ_nn/∂t + (λ+2μ) ∂u_n/∂n = 0,
  //   ∂σ_t1t1/∂t + λ ∂u_n/∂n = 0,
  //   ∂σ_t2t2/∂t + λ ∂u_n/∂n = 0,
  //   ∂σ_nt1/∂t + μ ∂u_t1/∂n = 0,
  //   ∂σ_nt2/∂t + μ ∂u_t2/∂n = 0,
  //   ∂u_n/∂t  + (1/ρ) ∂σ_nn/∂n = 0,
  //   ∂u_t1/∂t + (1/ρ) ∂σ_nt1/∂n = 0,
  //   ∂u_t2/∂t + (1/ρ) ∂σ_nt2/∂n = 0.
  // Hence numerical flux components (evaluated at Godunov state) are:
  //   F_n(σ_nn)   = (λ+2μ) u_n^G, etc.,  F_n(u_n) = σ_nn^G / ρ. (Consistent with linear elasticity)
  // NOTE: Using minus-side coefficients (assumed homogeneous patch). For jumps in material
  // coefficients, a two-sided consistent flux (e.g., average or impedance-weighted) is needed.
  const Real g_nn   = (_lambda[_qp] + 2.0*_mu[_qp]) * unG;   // (λ+2μ) u_n^G
  const Real g_t1t1 = (_lambda[_qp])                * unG;   // λ u_n^G  
  const Real g_t2t2 = (_lambda[_qp])                * unG;   // λ u_n^G
  const Real g_nt1  = (_mu[_qp])                    * ut1G;  // μ u_t1^G
  const Real g_nt2  = (_mu[_qp])                    * ut2G;  // μ u_t2^G

  // Velocity fluxes
  const Real f_un = snnG/_rho[_qp];
  const Real f_ut1 = snt1G/_rho[_qp];
  const Real f_ut2 = snt2G/_rho[_qp];

  // 6) Transform numerical fluxes back to global coordinate system
  // Equation (16): Apply inverse rotation to obtain flux in global (x,y,z) frame
  Real f_sxx, f_sxy, f_sxz, f_syy, f_syz, f_szz;
  ten_to_global(T,TT, g_nn, g_nt1, g_nt2, g_t1t1, /*st1t2*/0.0, g_t2t2,
                f_sxx,f_sxy,f_sxz,f_syy,f_syz,f_szz);

  // Transform velocity fluxes back to global coordinates
  RealVectorValue flux_vel_glob; 
  vec_to_global(T, f_un, f_ut1, f_ut2, flux_vel_glob);
  const Real f_ux = flux_vel_glob(0), f_uy = flux_vel_glob(1), f_uz = flux_vel_glob(2);

  // 7) Assemble DG contribution (standard interior face term)
  // Weak form (interface part):  ∫_Γ F_n^*( U^- , U^+ ) ( ψ^- - ψ^+ ) dΓ
  // MOOSE DGKernel convention: element side gets -F* ψ^-, neighbor +F* ψ^+.
  // (Consistent with sign usage in other DG kernels.)
  const std::string & var = _var.name();
  Real F = 0.0;
  if (var=="sxx") F=f_sxx; else if (var=="sxy") F=f_sxy; else if (var=="sxz") F=f_sxz;
  else if (var=="syy") F=f_syy; else if (var=="syz") F=f_syz; else if (var=="szz") F=f_szz;
  else if (var=="ux") F=f_ux; else if (var=="uy") F=f_uy; else if (var=="uz") F=f_uz;

  return (type==Moose::Element ? -_test[_i][_qp]*F : _test_neighbor[_i][_qp]*F);
}

Real ElasticGodunovDGFlux3D::computeQpJacobian(Moose::DGJacobianType type)
{
  // Analytical Jacobian of the Godunov flux w.r.t. the current variable
  // Mirrors computeQpResidual with linearized contributions per side and variable.

  // Build local frame and rotation matrices
  RealVectorValue n_glob = _normals[_qp], n, t1, t2;
  orthonormal_basis(n_glob, n, t1, t2);
  DenseMatrix<Real> T, TT; build_T(n, t1, t2, T, TT);

  // Minus-side material properties are used in flux evaluation (per implementation in residual)
  Real cp_m, cs_m, Zp_m, Zs_m;
  impedances(_lambda[_qp], _mu[_qp], _rho[_qp], cp_m, cs_m, Zp_m, Zs_m);

  const Real rho = _rho[_qp];
  const Real mu  = _mu[_qp];
  const Real lam = _lambda[_qp];

  // Helper to compute dF/dU_side for the current variable.
  // side_sign = +1 for minus (element) side, -1 for plus (neighbor) side.
  auto dF_dU_for_side = [&](int side_sign) -> Real
  {
    const std::string & var = _var.name();

    // Stress variables: map unit perturbation in the chosen component to local (snn, snt1, snt2)
    auto stress_component_to_local = [&](Real & dsnn, Real & dsnt1, Real & dsnt2)
    {
      Real sxx=0, sxy=0, sxz=0, syy=0, syz=0, szz=0;
      if (var=="sxx") sxx=1.0; else if (var=="sxy") sxy=1.0; else if (var=="sxz") sxz=1.0;
      else if (var=="syy") syy=1.0; else if (var=="syz") syz=1.0; else if (var=="szz") szz=1.0;
      Real st11, st12, st22; // unused here but required by API
      ten_to_local(T, TT, sxx, sxy, sxz, syy, syz, szz, dsnn, dsnt1, dsnt2, st11, st12, st22);
    };

    // Velocity variables: map unit perturbation in global velocity to local (un, ut1, ut2)
    auto velocity_component_to_local = [&](Real & dun, Real & dut1, Real & dut2)
    {
      RealVectorValue du(0,0,0);
      if (var=="ux") du(0)=1.0; else if (var=="uy") du(1)=1.0; else if (var=="uz") du(2)=1.0;
      vec_to_local(TT, du, dun, dut1, dut2);
    };

    // Select between stress and velocity flux derivatives
    if (var=="sxx" || var=="sxy" || var=="sxz" || var=="syy" || var=="syz" || var=="szz")
    {
      // Local stress increments from unit global component
      Real dsnn=0.0, dsnt1=0.0, dsnt2=0.0; stress_component_to_local(dsnn, dsnt1, dsnt2);

      // Godunov coupling factors (linear, side dependent)
      const Real d_unG = 0.5 * side_sign * (1.0/(rho * cp_m)) * dsnn;
      const Real dg_nn   = (lam + 2.0*mu) * d_unG;
      const Real dg_t1t1 = lam * d_unG;
      const Real dg_t2t2 = lam * d_unG;
      const Real dg_nt1  = 0.5 * side_sign * cs_m * dsnt1;
      const Real dg_nt2  = 0.5 * side_sign * cs_m * dsnt2;

      // Map local stress-flux increments back to global components
      Real df_sxx, df_sxy, df_sxz, df_syy, df_syz, df_szz;
      ten_to_global(T, TT, dg_nn, dg_nt1, dg_nt2, dg_t1t1, /*st1t2*/0.0, dg_t2t2,
                    df_sxx, df_sxy, df_sxz, df_syy, df_syz, df_szz);

      if (var=="sxx") return df_sxx;
      if (var=="sxy") return df_sxy;
      if (var=="sxz") return df_sxz;
      if (var=="syy") return df_syy;
      if (var=="syz") return df_syz;
      /* var=="szz" */ return df_szz;
    }
    else
    {
      // Velocity variables: local velocity increments
      Real dun=0.0, dut1=0.0, dut2=0.0; velocity_component_to_local(dun, dut1, dut2);

      // Godunov velocity flux derivatives (linear, side dependent)
      const Real df_un  = 0.5 * side_sign * cp_m * dun;
      const Real df_ut1 = 0.5 * side_sign * cs_m * dut1;
      const Real df_ut2 = 0.5 * side_sign * cs_m * dut2;

      // Back to global velocity flux components
      RealVectorValue df_vel_glob; vec_to_global(T, df_un, df_ut1, df_ut2, df_vel_glob);

      if (var=="ux") return df_vel_glob(0);
      if (var=="uy") return df_vel_glob(1);
      /* var=="uz" */ return df_vel_glob(2);
    }
  };

  // Assemble contribution per DGJacobianType
  switch (type)
  {
    case Moose::ElementElement:
    {
      const Real dF = dF_dU_for_side(+1);
      return -_test[_i][_qp] * dF * _phi[_j][_qp];
    }
    case Moose::NeighborNeighbor:
    {
      const Real dF = dF_dU_for_side(-1);
      return _test_neighbor[_i][_qp] * dF * _phi_neighbor[_j][_qp];
    }
    case Moose::ElementNeighbor:
    {
      const Real dF = dF_dU_for_side(-1);
      return -_test[_i][_qp] * dF * _phi_neighbor[_j][_qp];
    }
    case Moose::NeighborElement:
    {
      const Real dF = dF_dU_for_side(+1);
      return _test_neighbor[_i][_qp] * dF * _phi[_j][_qp];
    }
  }

  mooseError("Internal error.");
}

Real ElasticGodunovDGFlux3D::computeQpOffDiagJacobian(Moose::DGJacobianType /*type*/, unsigned int /*jvar*/)
{
  // TODO: Provide off-diagonal coupling terms. Currently using finite difference / zero contribution.
  return 0.0;
}
