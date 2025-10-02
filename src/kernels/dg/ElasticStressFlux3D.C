#include "ElasticStressFlux3D.h"

registerMooseObject("farmsApp", ElasticStressFlux3D);

InputParameters ElasticStressFlux3D::validParams()
{
  InputParameters p = Kernel::validParams();
  p.addRequiredCoupledVar("ux", "");
  p.addRequiredCoupledVar("uy", "");
  p.addRequiredCoupledVar("uz", "");
  p.addRequiredParam<MaterialPropertyName>("lambda_name", "lambda material property name");
  p.addRequiredParam<MaterialPropertyName>("mu_name", "shear modulus material property name");
  p.addRequiredParam<std::string>("component", "Stress component: sxx syy szz sxy sxz syz");
  return p;
}

ElasticStressFlux3D::ElasticStressFlux3D(const InputParameters & p)
  : Kernel(p),
    _ux(coupledValue("ux")), _uy(coupledValue("uy")), _uz(coupledValue("uz")),
    _lambda(getMaterialProperty<Real>(getParam<MaterialPropertyName>("lambda_name"))),
    _mu(getMaterialProperty<Real>(getParam<MaterialPropertyName>("mu_name"))),
    _comp([&](){
      const std::string c = getParam<std::string>("component");
      if (c=="sxx") return 0u; if (c=="syy") return 1u; if (c=="szz") return 2u; //0u means unsigned int u
      if (c=="sxy") return 3u; if (c=="sxz") return 4u; if (c=="syz") return 5u;
      mooseError("ElasticStressFlux3D invalid component ", c);
    }()),
    _ux_var(coupled("ux")), _uy_var(coupled("uy")), _uz_var(coupled("uz"))
{}

Real ElasticStressFlux3D::computeQpResidual()
{
  // Material coefficients at this quadrature point
  const Real lam = _lambda[_qp];
  const Real mu  = _mu[_qp];

  // Flux vector F whose divergence appears with a plus sign in the strong form
  //  dt(s) + div(F) = 0. We evaluate -grad(test) dot F (volume contribution after integration by parts).
  Real Fx=0.0, Fy=0.0, Fz=0.0;
  switch (_comp)
  {
    case 0: // s_xx equation
      Fx = (lam+2*mu)*_ux[_qp]; Fy = lam*_uy[_qp];      Fz = lam*_uz[_qp];      break;
    case 1: // s_yy equation
      Fx = lam*_ux[_qp];        Fy = (lam+2*mu)*_uy[_qp]; Fz = lam*_uz[_qp];    break;
    case 2: // s_zz equation
      Fx = lam*_ux[_qp];        Fy = lam*_uy[_qp];      Fz = (lam+2*mu)*_uz[_qp]; break;
    case 3: // s_xy equation: dt(s_xy) - mu*(d_x v + d_y u) = 0, so Fx = mu*uy, Fy = mu*ux
      Fx = mu*_uy[_qp];         Fy = mu*_ux[_qp];       break;
    case 4: // s_xz equation: dt(s_xz) - mu*(d_x w + d_z u) = 0
      Fx = mu*_uz[_qp];         Fz = mu*_ux[_qp];       break;
    case 5: // s_yz equation: dt(s_yz) - mu*(d_y w + d_z v) = 0
      Fy = mu*_uz[_qp];         Fz = mu*_uy[_qp];       break;
    default:
      mooseError("ElasticStressFlux3D missing flux definition for component index ", _comp);
  }

  // grad(test) dot F
  const Real divF = _grad_test[_i][_qp](0)*Fx + _grad_test[_i][_qp](1)*Fy + _grad_test[_i][_qp](2)*Fz;
  return divF;
}

Real ElasticStressFlux3D::computeQpJacobian()
{
  // Flux depends only on coupled velocities, not on the stress variable itself in this first-order form.
  return 0.0;
}

Real ElasticStressFlux3D::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Diagonal Jacobian is zero (flux depends on velocities only), so only provide
  // off-diagonal couplings with ux, uy, uz.
  const Real lam = _lambda[_qp];
  const Real mu  = _mu[_qp];

  Real dFx_du = 0.0, dFy_du = 0.0, dFz_du = 0.0; // partials wrt selected velocity

  // Select which velocity this jvar corresponds to and set the appropriate flux derivatives
  if (jvar == _ux_var)
  {
    switch (_comp)
    {
      case 0: dFx_du = (lam + 2.0*mu);            break; // sxx: Fx = (lam+2mu) ux
      case 1: dFx_du = lam;                       break; // syy: Fx = lam ux
      case 2: dFx_du = lam;                       break; // szz: Fx = lam ux
      case 3: dFy_du = mu;                        break; // sxy: Fy = mu ux
      case 4: dFz_du = mu;                        break; // sxz: Fz = mu ux
      case 5: /* syz: no ux */                    break;
    }
  }
  else if (jvar == _uy_var)
  {
    switch (_comp)
    {
      case 0: dFy_du = lam;                       break; // sxx: Fy = lam uy
      case 1: dFy_du = (lam + 2.0*mu);            break; // syy: Fy = (lam+2mu) uy
      case 2: dFy_du = lam;                       break; // szz: Fy = lam uy
      case 3: dFx_du = mu;                        break; // sxy: Fx = mu uy
      case 4: /* sxz: no uy */                    break;
      case 5: dFz_du = mu;                        break; // syz: Fz = mu uy
    }
  }
  else if (jvar == _uz_var)
  {
    switch (_comp)
    {
      case 0: dFz_du = lam;                       break; // sxx: Fz = lam uz
      case 1: dFz_du = lam;                       break; // syy: Fz = lam uz
      case 2: dFz_du = (lam + 2.0*mu);            break; // szz: Fz = (lam+2mu) uz
      case 3: /* sxy: no uz */                    break;
      case 4: dFx_du = mu;                        break; // sxz: Fx = mu uz
      case 5: dFy_du = mu;                        break; // syz: Fy = mu uz
    }
  }
  else
    return 0.0; // not a coupled velocity variable

  // Chain rule: dR = grad(test) dot dF, with dF = [dFx_du, dFy_du, dFz_du] * phi_j
  const Real dR = _grad_test[_i][_qp](0) * dFx_du +
                  _grad_test[_i][_qp](1) * dFy_du +
                  _grad_test[_i][_qp](2) * dFz_du;
  return dR * _phi[_j][_qp];
}
