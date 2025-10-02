#include "FaultSlipWeakeningDG3DMaterial.h"
using namespace ElasticDG3D;

registerMooseObject("farmsApp", FaultSlipWeakeningDG3DMaterial);

InputParameters FaultSlipWeakeningDG3DMaterial::validParams()
{
  InputParameters params = InterfaceMaterial::validParams();
  params.addRequiredCoupledVar("ux","x velocity"); params.addRequiredCoupledVar("uy","y velocity"); params.addRequiredCoupledVar("uz","z velocity");
  params.addRequiredCoupledVar("sxx",""); params.addRequiredCoupledVar("syy",""); params.addRequiredCoupledVar("szz","");
  params.addRequiredCoupledVar("sxy",""); params.addRequiredCoupledVar("sxz",""); params.addRequiredCoupledVar("syz","");
  // No external slip variable coupling; slip maintained internally as stateful material property
  params.addRequiredParam<MaterialPropertyName>("lambda_name","λ prop name");
  params.addRequiredParam<MaterialPropertyName>("mu_name","μ prop name");
  params.addRequiredParam<MaterialPropertyName>("rho_name","ρ prop name");
  params.addRequiredParam<Real>("mu_d","dynamic friction");
  params.addRequiredParam<Real>("Dc","critical slip");
  params.addParam<Real>("mu_s",0.0,"static friction (if no mu_s_aux)");
  params.addParam<Real>("tau0_t1",0.0,"initial shear traction offset t1 (if no tau0_t1_aux)");
  params.addParam<Real>("tau0_t2",0.0,"initial shear traction offset t2");
  params.addParam<Real>("sigma0n",0.0,"initial normal stress offset (if no sigma0n_aux)");
  params.addParam<VariableName>("mu_s_aux","Aux variable for static friction (optional)");
  params.addParam<VariableName>("tau0_t1_aux","Aux variable initial shear t1 (optional)");
  params.addParam<VariableName>("sigma0n_aux","Aux variable normal stress offset (optional)");
  return params;
}

// NOTE: Initializer order MUST match member declaration order in header to avoid -Wreorder (treated as error)
FaultSlipWeakeningDG3DMaterial::FaultSlipWeakeningDG3DMaterial(const InputParameters & p)
  : InterfaceMaterial(p),
    // Coupled primary variables (minus / plus)
    _ux(coupledValue("ux")), _ux_n(coupledNeighborValue("ux")),
    _uy(coupledValue("uy")), _uy_n(coupledNeighborValue("uy")),
    _uz(coupledValue("uz")), _uz_n(coupledNeighborValue("uz")),
    _sxx(coupledValue("sxx")), _sxx_n(coupledNeighborValue("sxx")),
    _syy(coupledValue("syy")), _syy_n(coupledNeighborValue("syy")),
    _szz(coupledValue("szz")), _szz_n(coupledNeighborValue("szz")),
    _sxy(coupledValue("sxy")), _sxy_n(coupledNeighborValue("sxy")),
    _sxz(coupledValue("sxz")), _sxz_n(coupledNeighborValue("sxz")),
    _syz(coupledValue("syz")), _syz_n(coupledNeighborValue("syz")),
    // Elastic material properties (minus / plus)
    _lambda(getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("lambda_name"))),
    _lambda_n(getNeighborMaterialProperty<Real>(getParam<MaterialPropertyName>("lambda_name"))),
    _mu(getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("mu_name"))),
    _mu_n(getNeighborMaterialProperty<Real>(getParam<MaterialPropertyName>("mu_name"))),
    _rho(getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("rho_name"))),
    _rho_n(getNeighborMaterialProperty<Real>(getParam<MaterialPropertyName>("rho_name"))),
    // Friction constants
    _mu_s_const(getParam<Real>("mu_s")), _mu_d(getParam<Real>("mu_d")), _Dc(getParam<Real>("Dc")),
    _tau0_t1_const(getParam<Real>("tau0_t1")), _tau0_t2_const(getParam<Real>("tau0_t2")), _sigma0n_const(getParam<Real>("sigma0n")),
    // Stateful slip (declare AFTER friction constants per header order)
    _slip(declareProperty<Real>("fault_slip")), _slip_old(getMaterialPropertyOld<Real>("fault_slip")),
    // Flux properties
    _flux_sxx(declareProperty<Real>("fault_flux_sxx")), _flux_sxy(declareProperty<Real>("fault_flux_sxy")), _flux_sxz(declareProperty<Real>("fault_flux_sxz")),
    _flux_syy(declareProperty<Real>("fault_flux_syy")), _flux_syz(declareProperty<Real>("fault_flux_syz")), _flux_szz(declareProperty<Real>("fault_flux_szz")),
    _flux_ux(declareProperty<Real>("fault_flux_ux")), _flux_uy(declareProperty<Real>("fault_flux_uy")), _flux_uz(declareProperty<Real>("fault_flux_uz")),
    // Save properties
    _traction_t1(declareProperty<Real>("fault_tau_t1")), _traction_t2(declareProperty<Real>("fault_tau_t2")), _traction_mag(declareProperty<Real>("fault_tau_mag")),
    _normal_stress(declareProperty<Real>("fault_sigma_n")), _mu_f_prop(declareProperty<Real>("fault_mu_f")),
    _sr_t1(declareProperty<Real>("fault_slip_rate_t1")), _sr_t2(declareProperty<Real>("fault_slip_rate_t2")), _sr_mag(declareProperty<Real>("fault_slip_rate_mag")),
    _ut1_star(declareProperty<Real>("fault_ut1_star")), _ut2_star(declareProperty<Real>("fault_ut2_star"))
{
  if (isParamValid("mu_s_aux")) _mu_s_aux = &coupledValue("mu_s_aux");
  if (isParamValid("tau0_t1_aux")) _tau0_t1_aux = &coupledValue("tau0_t1_aux");
  if (isParamValid("sigma0n_aux")) _sigma0n_aux = &coupledValue("sigma0n_aux");
}

void FaultSlipWeakeningDG3DMaterial::initQpStatefulProperties()
{ _slip[_qp] = 0.0; }

void FaultSlipWeakeningDG3DMaterial::computeQpProperties()
{
  // Build local frame
  RealVectorValue n_glob = _normals[_qp], n,t1,t2; orthonormal_basis(n_glob,n,t1,t2);
  DenseMatrix<Real> T,TT; build_T(n,t1,t2,T,TT);

  // Minus-side material and impedances
  const Real lam_m = _lambda[_qp]; const Real mu_m = _mu[_qp]; const Real rho_m = _rho[_qp];
  Real cp_m, cs_m, Zp_m, Zs_m; impedances(lam_m, mu_m, rho_m, cp_m, cs_m, Zp_m, Zs_m);
  // Plus-side (not currently used for heterogeneous formulation except diagnostics)
  Real cp_p, cs_p, Zp_p, Zs_p; impedances(_lambda_n[_qp], _mu_n[_qp], _rho_n[_qp], cp_p, cs_p, Zp_p, Zs_p);

  // Rotate stresses & velocities to local
  Real snn_m,snt1_m,snt2_m, st11_m,st12_m,st22_m;
  Real snn_p,snt1_p,snt2_p, st11_p,st12_p,st22_p;
  ten_to_local(T,TT,_sxx[_qp],_sxy[_qp],_sxz[_qp],_syy[_qp],_syz[_qp],_szz[_qp],
               snn_m,snt1_m,snt2_m,st11_m,st12_m,st22_m);
  ten_to_local(T,TT,_sxx_n[_qp],_sxy_n[_qp],_sxz_n[_qp],_syy_n[_qp],_syz_n[_qp],_szz_n[_qp],
               snn_p,snt1_p,snt2_p,st11_p,st12_p,st22_p);

  Real un_m,ut1_m,ut2_m, un_p,ut1_p,ut2_p;
  vec_to_local(TT, RealVectorValue(_ux[_qp],_uy[_qp],_uz[_qp]), un_m, ut1_m, ut2_m);
  vec_to_local(TT, RealVectorValue(_ux_n[_qp],_uy_n[_qp],_uz_n[_qp]), un_p, ut1_p, ut2_p);

  // Godunov states
  const Real unG   = 0.5 * ( (un_p + un_m) + (1.0/(cp_m * rho_m)) * (snn_m - snn_p) );
  const Real snnG  = 0.5 * ( (snn_p + snn_m) + (cp_m * rho_m)     * (un_m - un_p) );
  // Tangential Godunov stresses (velocities not used directly later, omit ut1G/ut2G to silence -Wunused warnings)
  const Real snt1G = 0.5 * ( (snt1_p + snt1_m) + (mu_m / cs_m) * (ut1_m - ut1_p) );
  const Real snt2G = 0.5 * ( (snt2_p + snt2_m) + (mu_m / cs_m) * (ut2_m - ut2_p) );

  // Friction overrides
  Real mu_s_val = _mu_s_const; if (_mu_s_aux) mu_s_val = (*_mu_s_aux)[_qp];
  Real tau0_t1_val = _tau0_t1_const; if (_tau0_t1_aux) tau0_t1_val = (*_tau0_t1_aux)[_qp];
  Real sigma0n_val = _sigma0n_const; if (_sigma0n_aux) sigma0n_val = (*_sigma0n_aux)[_qp];
  const Real mu_f = weakeningMu(_slip_old[_qp], mu_s_val); // weakening uses OLD slip
  const Real sigmaN = snnG + sigma0n_val;
  Real tauG_t1 = snt1G + tau0_t1_val;
  Real tauG_t2 = snt2G + _tau0_t2_const;
  const Real tauG_mag = std::sqrt(tauG_t1*tauG_t1 + tauG_t2*tauG_t2);
  const Real strength = mu_f * std::abs(sigmaN);
  Real tauimp_t1 = tauG_t1, tauimp_t2 = tauG_t2;
  if (tauG_mag > strength && tauG_mag > 0.0)
  { const Real scale = strength / tauG_mag; tauimp_t1*=scale; tauimp_t2*=scale; }

  // Radiation damping updates
  const Real factor = cs_m / mu_m;
  const Real dut1_p =  factor * (tauimp_t1 - snt1_p);
  const Real dut1_m = -factor * (tauimp_t1 - snt1_m);
  const Real dut2_p =  factor * (tauimp_t2 - snt2_p);
  const Real dut2_m = -factor * (tauimp_t2 - snt2_m);
  const Real ut1_tilde_p = ut1_p + dut1_p;
  const Real ut1_tilde_m = ut1_m + dut1_m;
  const Real ut2_tilde_p = ut2_p + dut2_p;
  const Real ut2_tilde_m = ut2_m + dut2_m;
  const Real ut1_star = 0.5*(ut1_tilde_p + ut1_tilde_m);
  const Real ut2_star = 0.5*(ut2_tilde_p + ut2_tilde_m);
  const Real sr_t1 = ut1_tilde_m - ut1_tilde_p; // slip-rate components
  const Real sr_t2 = ut2_tilde_m - ut2_tilde_p;
  const Real sr_mag = std::sqrt(sr_t1*sr_t1 + sr_t2*sr_t2);

  // Compute Slip
  _slip[_qp] = _slip_old[_qp] + sr_mag * _dt;

  // Local fluxes (same sign convention as kernel implementation earlier)
  const Real g_nn   = -(lam_m + 2.0*mu_m) * unG;
  const Real g_t1t1 = -(lam_m)            * unG;
  const Real g_t2t2 = -(lam_m)            * unG;
  const Real g_nt1  = -(mu_m)             * ut1_star;
  const Real g_nt2  = -(mu_m)             * ut2_star;
  const Real t_n = snnG; const Real t_t1 = tauimp_t1; const Real t_t2 = tauimp_t2;

  // Rotate stress fluxes back
  Real f_sxx,f_sxy,f_sxz,f_syy,f_syz,f_szz;
  ten_to_global(T,TT, g_nn, g_nt1, g_nt2, g_t1t1, /*st1t2*/0.0, g_t2t2,
                f_sxx,f_sxy,f_sxz,f_syy,f_syz,f_szz);
  RealVectorValue traction_glob; vec_to_global(T, t_n, t_t1, t_t2, traction_glob);
  const Real f_ux = -traction_glob(0), f_uy = -traction_glob(1), f_uz = -traction_glob(2);

  _flux_sxx[_qp]=f_sxx; _flux_sxy[_qp]=f_sxy; _flux_sxz[_qp]=f_sxz;
  _flux_syy[_qp]=f_syy; _flux_syz[_qp]=f_syz; _flux_szz[_qp]=f_szz;
  _flux_ux[_qp]=f_ux;   _flux_uy[_qp]=f_uy;   _flux_uz[_qp]=f_uz;
  _traction_t1[_qp]=tauimp_t1; _traction_t2[_qp]=tauimp_t2; _traction_mag[_qp]=std::sqrt(tauimp_t1*tauimp_t1 + tauimp_t2*tauimp_t2);
  _normal_stress[_qp]=sigmaN; _mu_f_prop[_qp]=mu_f;
  _sr_t1[_qp]=sr_t1; _sr_t2[_qp]=sr_t2; _sr_mag[_qp]=sr_mag;
  _ut1_star[_qp]=ut1_star; _ut2_star[_qp]=ut2_star;
}
