#include "FaultGodunovStateUO.h"
#include "MooseMesh.h"

// Use elastic DG utility free functions in namespace ElasticDG3D
using namespace ElasticDG3D;

registerMooseObject("farmsApp", FaultGodunovStateUO);

InputParameters FaultGodunovStateUO::validParams()
{
  InputParameters params = InterfaceUserObject::validParams();
  params.addRequiredCoupledVar("ux", ""); params.addRequiredCoupledVar("uy", ""); params.addRequiredCoupledVar("uz", "");
  params.addRequiredCoupledVar("sxx", ""); params.addRequiredCoupledVar("syy", ""); params.addRequiredCoupledVar("szz", "");
  params.addRequiredCoupledVar("sxy", ""); params.addRequiredCoupledVar("sxz", ""); params.addRequiredCoupledVar("syz", "");
  params.addRequiredCoupledVar("slip", "accumulated slip magnitude (Aux)");
  params.addRequiredParam<MaterialPropertyName>("lambda_name", "λ");
  params.addRequiredParam<MaterialPropertyName>("mu_name", "μ");
  params.addRequiredParam<MaterialPropertyName>("rho_name", "ρ");
  params.addRequiredParam<Real>("mu_s","static friction");
  params.addRequiredParam<Real>("mu_d","dynamic friction");
  params.addRequiredParam<Real>("Dc","critical slip");
  params.addParam<Real>("sigma0n",0.0,"initial normal stress offset");
  params.addParam<Real>("tau0_t1",0.0,"initial shear traction offset along t1");
  params.addParam<Real>("tau0_t2",0.0,"initial shear traction offset along t2");
  return params;
}

FaultGodunovStateUO::FaultGodunovStateUO(const InputParameters & params)
 : InterfaceUserObject(params),
   _lambda(getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("lambda_name"))),
   _lambda_n(getNeighborMaterialProperty<Real>(getParam<MaterialPropertyName>("lambda_name"))),
   _mu(getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("mu_name"))),
   _mu_n(getNeighborMaterialProperty<Real>(getParam<MaterialPropertyName>("mu_name"))),
   _rho(getMaterialPropertyByName<Real>(getParam<MaterialPropertyName>("rho_name"))),
   _rho_n(getNeighborMaterialProperty<Real>(getParam<MaterialPropertyName>("rho_name"))),
   _ux(coupledValue("ux")), _ux_n(coupledNeighborValue("ux")),
   _uy(coupledValue("uy")), _uy_n(coupledNeighborValue("uy")),
   _uz(coupledValue("uz")), _uz_n(coupledNeighborValue("uz")),
   _sxx(coupledValue("sxx")), _sxx_n(coupledNeighborValue("sxx")),
   _syy(coupledValue("syy")), _syy_n(coupledNeighborValue("syy")),
   _szz(coupledValue("szz")), _szz_n(coupledNeighborValue("szz")),
   _sxy(coupledValue("sxy")), _sxy_n(coupledNeighborValue("sxy")),
   _sxz(coupledValue("sxz")), _sxz_n(coupledNeighborValue("sxz")),
   _syz(coupledValue("syz")), _syz_n(coupledNeighborValue("syz")),
   _slip_old(coupledValueOld("slip")),
   _mu_s(getParam<Real>("mu_s")), _mu_d(getParam<Real>("mu_d")), _Dc(getParam<Real>("Dc")),
   _tau0_t1(getParam<Real>("tau0_t1")), _tau0_t2(getParam<Real>("tau0_t2")), _sigma0n(getParam<Real>("sigma0n"))
{}

void FaultGodunovStateUO::initialize()
{
  _states.clear();
  _lookup.clear();
}

void FaultGodunovStateUO::execute()
{
  // Loop over current face quadrature points; InterfaceUserObject provides _q_point and _normals
  const unsigned short n_qp = _q_point.size();
  const Elem * elem = _current_elem; // minus side element
  const unsigned short side = _current_side;

  for (unsigned short qp=0; qp<n_qp; ++qp)
  {
    // 1) Local frame
    RealVectorValue n_glob = _normals[qp], n,t1,t2; orthonormal_basis(n_glob,n,t1,t2);
    DenseMatrix<Real> T,TT; build_T(n,t1,t2,T,TT);

    // 2) Impedances (minus/plus)
    Real cp_m, cs_m, Zp_m, Zs_m; impedances(_lambda[qp], _mu[qp], _rho[qp], cp_m, cs_m, Zp_m, Zs_m);
    Real cp_p, cs_p, Zp_p, Zs_p; impedances(_lambda_n[qp], _mu_n[qp], _rho_n[qp], cp_p, cs_p, Zp_p, Zs_p);

    // 3) Traces -> local
    Real snn_m,snt1_m,snt2_m, st11_m,st12_m,st22_m;
    Real snn_p,snt1_p,snt2_p, st11_p,st12_p,st22_p;
    ten_to_local(T,TT,_sxx[qp],_sxy[qp],_sxz[qp],_syy[qp],_syz[qp],_szz[qp],
                 snn_m,snt1_m,snt2_m,st11_m,st12_m,st22_m);
    ten_to_local(T,TT,_sxx_n[qp],_sxy_n[qp],_sxz_n[qp],_syy_n[qp],_syz_n[qp],_szz_n[qp],
                 snn_p,snt1_p,snt2_p,st11_p,st12_p,st22_p);

    Real un_m,ut1_m,ut2_m, un_p,ut1_p,ut2_p;
    vec_to_local(TT, RealVectorValue(_ux[qp],_uy[qp],_uz[qp]), un_m, ut1_m, ut2_m);
    vec_to_local(TT, RealVectorValue(_ux_n[qp],_uy_n[qp],_uz_n[qp]), un_p, ut1_p, ut2_p);

    // 4) Godunov states (Eq.13)
    const Real unG   = 0.5 * ( (un_p + un_m) + (1.0/(cp_m * _rho[qp])) * (snn_m - snn_p) );
    const Real snnG  = 0.5 * ( (snn_p + snn_m) + (cp_m * _rho[qp])     * (un_m - un_p) );
    const Real ut1G  = 0.5 * ( (ut1_p + ut1_m) + (cs_m / _mu[qp]) * (snt1_m - snt1_p) );
    const Real snt1G = 0.5 * ( (snt1_p + snt1_m) + (_mu[qp] / cs_m) * (ut1_m - ut1_p) );
    const Real ut2G  = 0.5 * ( (ut2_p + ut2_m) + (cs_m / _mu[qp]) * (snt2_m - snt2_p) );
    const Real snt2G = 0.5 * ( (snt2_p + snt2_m) + (_mu[qp] / cs_m) * (ut2_m - ut2_p) );

    // 5) Friction projection
    const Real mu_f = mu_slipweak(_slip_old[qp]);
    const Real sigmaN_eff = snnG + _sigma0n;
    Real tauG_t1 = snt1G + _tau0_t1;
    Real tauG_t2 = snt2G + _tau0_t2;
    const Real tauG_mag = std::sqrt(tauG_t1*tauG_t1 + tauG_t2*tauG_t2);
    const Real strength = mu_f * std::abs(sigmaN_eff);
    Real tauimp_t1 = tauG_t1, tauimp_t2 = tauG_t2;
    if (tauG_mag > strength && tauG_mag>0.0)
    { const Real scale = strength / tauG_mag; tauimp_t1*=scale; tauimp_t2*=scale; }

    // 6) Radiation damping (Eq.15–16)
    const Real factor = cs_m / _mu[qp];
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

  // 7) Slip-rate components in local tangential directions (minus - plus) using star velocities
  const Real sr_t1 = ut1_tilde_m - ut1_tilde_p; // relative tangential velocity t1
  const Real sr_t2 = ut2_tilde_m - ut2_tilde_p; // relative tangential velocity t2

  // Store
  State st{ snnG, snt1G, snt2G, unG, ut1_star, ut2_star, tauimp_t1, tauimp_t2, sigmaN_eff, sr_t1, sr_t2 };
    unsigned long long key = pack(elem->id(), side, qp);
    _lookup[key] = _states.size();
    _states.push_back(st);
  }
}

void FaultGodunovStateUO::threadJoin(const UserObject & y)
{
  const auto & o = static_cast<const FaultGodunovStateUO &>(y);
  const std::size_t offset = _states.size();
  _states.insert(_states.end(), o._states.begin(), o._states.end());
  for (const auto & kv : o._lookup)
    _lookup.emplace(kv.first, kv.second + offset);
}

const FaultGodunovStateUO::State & FaultGodunovStateUO::get(dof_id_type e, unsigned short side, unsigned short qp) const
{
  auto it = _lookup.find(pack(e,side,qp));
  if (it == _lookup.end())
    mooseError("FaultGodunovStateUO: state not found for elem=", e, " side=", side, " qp=", qp);
  return _states[it->second];
}

bool FaultGodunovStateUO::has(dof_id_type e, unsigned short side, unsigned short qp) const
{
  return _lookup.count(pack(e,side,qp));
}
