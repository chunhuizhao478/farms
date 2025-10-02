#include "ElasticVelocityFlux3D.h"

registerMooseObject("farmsApp", ElasticVelocityFlux3D);

InputParameters ElasticVelocityFlux3D::validParams()
{
  InputParameters p = Kernel::validParams();
  p.addRequiredCoupledVar("sxx", "");
  p.addRequiredCoupledVar("syy", "");
  p.addRequiredCoupledVar("szz", "");
  p.addRequiredCoupledVar("sxy", "");
  p.addRequiredCoupledVar("sxz", "");
  p.addRequiredCoupledVar("syz", "");
  p.addRequiredParam<std::string>("component", "Velocity component: ux, uy, or uz");
  return p;
}

ElasticVelocityFlux3D::ElasticVelocityFlux3D(const InputParameters & p)
  : Kernel(p),
    _sxx(coupledValue("sxx")), _syy(coupledValue("syy")), _szz(coupledValue("szz")),
    _sxy(coupledValue("sxy")), _sxz(coupledValue("sxz")), _syz(coupledValue("syz")),
    _comp([&](){ const std::string c=getParam<std::string>("component");
      if (c=="ux") return 0u; if (c=="uy") return 1u; if (c=="uz") return 2u;
      mooseError("ElasticVelocityFlux3D invalid component ", c); }()),
    _sxx_var(coupled("sxx")), _syy_var(coupled("syy")), _szz_var(coupled("szz")),
    _sxy_var(coupled("sxy")), _sxz_var(coupled("sxz")), _syz_var(coupled("syz"))
{}

Real ElasticVelocityFlux3D::computeQpResidual()
{
  // Build the stress column corresponding to the selected velocity component.
  // We treat the governing equation as rho * dt(u) - sum_j d_j sigma_{j x} = 0
  // so dt(u) + div(F) = 0 with F_j = sigma_{j x}. The residual contribution is
  // grad(test) dot F; adjust signs if your existing kernels use a different convention.
  Real Fx=0.0, Fy=0.0, Fz=0.0;
  switch (_comp)
  {
    case 0: // ux equation uses s_xx, s_xy, s_xz
      Fx = _sxx[_qp]; Fy = _sxy[_qp]; Fz = _sxz[_qp]; break;
    case 1: // uy equation uses s_xy, s_yy, s_yz
      Fx = _sxy[_qp]; Fy = _syy[_qp]; Fz = _syz[_qp]; break;
    case 2: // uz equation uses s_xz, s_yz, s_zz
      Fx = _sxz[_qp]; Fy = _syz[_qp]; Fz = _szz[_qp]; break;
    default:
      mooseError("ElasticVelocityFlux3D missing flux definition for component index ", _comp);
  }
  const Real term = _grad_test[_i][_qp](0)*Fx + _grad_test[_i][_qp](1)*Fy + _grad_test[_i][_qp](2)*Fz;
  return term;
}

Real ElasticVelocityFlux3D::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Residual R = grad(test) dot F, where F components are selected stresses.
  // For ux equation: F = (s_xx, s_xy, s_xz)
  // For uy equation: F = (s_xy, s_yy, s_yz)
  // For uz equation: F = (s_xz, s_yz, s_zz)
  // Derivatives dR/dsigma = grad_component(test) * phi_j when the stress participates.

  const Real gx = _grad_test[_i][_qp](0);
  const Real gy = _grad_test[_i][_qp](1);
  const Real gz = _grad_test[_i][_qp](2);

  Real coeff = 0.0;
  if (_comp == 0) // ux
  {
    if (jvar == _sxx_var) coeff = gx;
    else if (jvar == _sxy_var) coeff = gy;
    else if (jvar == _sxz_var) coeff = gz;
  }
  else if (_comp == 1) // uy
  {
    if (jvar == _sxy_var) coeff = gx;
    else if (jvar == _syy_var) coeff = gy;
    else if (jvar == _syz_var) coeff = gz;
  }
  else if (_comp == 2) // uz
  {
    if (jvar == _sxz_var) coeff = gx;
    else if (jvar == _syz_var) coeff = gy;
    else if (jvar == _szz_var) coeff = gz;
  }
  else
    mooseError("ElasticVelocityFlux3D invalid component index ", _comp);

  if (coeff == 0.0)
    return 0.0;
  return coeff * _phi[_j][_qp];
}
