
#include "INSmassFluid.h"
#include "Function.h"

registerMooseObject("farmsApp", INSmassFluid);

InputParameters
INSmassFluid::validParams()
{
  InputParameters params = INSBase::validParams();

  params.addClassDescription("This class computes the mass equation residual and Jacobian "
                             "contributions for the incompressible Navier-Stokes momentum "
                             "equation.");
  return params;
}

INSmassFluid::INSmassFluid(const InputParameters & parameters)
  : INSBase(parameters)
{
}

Real
INSmassFluid::computeQpResidual()
{
  // (div u) * q

  Real r = (_grad_u_vel[_qp](0) + _grad_v_vel[_qp](1) + _grad_w_vel[_qp](2)) * _test[_i][_qp];

  return r;
}

Real
INSmassFluid::computeQpJacobian()
{
  // Derivative wrt to p is zero
  Real r = 0;

  return r;
}

Real
INSmassFluid::computeQpOffDiagJacobian(const unsigned int jvar)
{
  if (jvar == _u_vel_var_number)
  {
    Real jac = _grad_phi[_j][_qp](0) * _test[_i][_qp];
    return jac;
  }

  else if (jvar == _v_vel_var_number)
  {
    Real jac = _grad_phi[_j][_qp](1) * _test[_i][_qp];
    return jac;
  }

  else if (jvar == _w_vel_var_number)
  {
    Real jac = _grad_phi[_j][_qp](2) * _test[_i][_qp];
    return jac;
  }

  else
    return 0.0;
}
