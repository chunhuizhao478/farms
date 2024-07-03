#include "FarmsCZM.h"
registerMooseObject("farmsApp", FarmsCZM);

InputParameters
FarmsCZM::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addClassDescription(
      "Interface Kernel for Linear Slip Weakening Friction Law");
  params.addRequiredCoupledVar("displacements", "the string containing displacement variables");
  return params;
}

FarmsCZM::FarmsCZM(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _traction_on_interface(getMaterialPropertyByName<RealVectorValue>("traction_on_interface")),
    _material_tangent_modulus_on_interface(
        getMaterialPropertyByName<RealTensorValue>("material_tangent_modulus_on_interface")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp)
{
  for (const auto i : make_range(_ndisp))
  {
    _disp_var[i] = coupled("displacements", i);
    if (_disp_var[i] == _var.number())
      _component = i;
  }
}

Real
FarmsCZM::computeQpResidual(Moose::DGResidualType type)
{
  Real r(0.0);
  
  r = _traction_on_interface[_qp](_component);
  switch (type)
  {
    case Moose::Element:
      r *= _test[_i][_qp];
      break;
    case Moose::Neighbor:
      r *= -_test_neighbor[_i][_qp];
      break;
  }

  return r;
}

Real
FarmsCZM::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;
  return jac;
}