#include "RhoTimeDerivative.h"

registerMooseObject("farmsApp", RhoTimeDerivative);

InputParameters RhoTimeDerivative::validParams()
{
  InputParameters p = TimeDerivative::validParams();
  p.addRequiredParam<MaterialPropertyName>("rho_name", "density material property");
  return p;
}

RhoTimeDerivative::RhoTimeDerivative(const InputParameters & p)
  : TimeDerivative(p),
    _rho(getMaterialProperty<Real>(getParam<MaterialPropertyName>("rho_name")))
{}

Real RhoTimeDerivative::computeQpResidual()
{ return _rho[_qp] * TimeDerivative::computeQpResidual(); }

Real RhoTimeDerivative::computeQpJacobian()
{ return _rho[_qp] * TimeDerivative::computeQpJacobian(); }
