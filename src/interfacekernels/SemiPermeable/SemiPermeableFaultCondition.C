#include "SemiPermeableFaultCondition.h"
#include "InterfaceKernel.h"

registerMooseObject("farmsApp", SemiPermeableFaultCondition);

InputParameters
SemiPermeableFaultCondition::validParams()
{
    InputParameters params = InterfaceKernel::validParams();
    params.addCoupledVar("pressure_main","Pressure at main side of the fault");
    params.addCoupledVar("pressure_sec","Pressure at secondary side of the fault");
    params.addParam<Real>("Transmissibility", 1.0, "Transmissibility of the fault boundary");
    return params;
}
SemiPermeableFaultCondition::SemiPermeableFaultCondition(const InputParameters &parameters)
  : InterfaceKernel(parameters), 
  _trans(getParam<Real>("Transmissibility")),
  _pore_pressure_main(coupledValue("pressure_main")),
  _pore_pressure_secondary(coupledNeighborValue("pressure_sec"))
{
}
Real
SemiPermeableFaultCondition::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;
  switch (type)
  {

    // Weak form for primary domain is: (test, T * (P_main - p_secondary) )
    case Moose::Element:
      r = _test[_i][_qp] * ( _u[_qp]  - _trans[_qp] * (_pore_pressure_main[_qp] - _pore_pressure_secondary[_qp]));
      break;

    // Similarly, weak form for secondary domain is: (test, - T * (P_main - p_secondary) )
    // flip the sign because the direction is opposite.
    case Moose::Neighbor:
      r = -_test_neighbor[_i][_qp] * ( _u[_qp]  - _trans[_qp] * (_pore_pressure_main[_qp] - _pore_pressure_secondary[_qp]));
      break;
  }
  return r;
}

Real
SemiPermeableFaultCondition::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;
  switch (type)
  {
    case Moose::ElementElement:
      jac = _test[_i][_qp] * _phi[_j][_qp];
      break;
    case Moose::NeighborNeighbor:
      jac = -_test_neighbor[_i][_qp] * _phi_neighbor[_j][_qp];
      break;
    case Moose::NeighborElement:
      jac = -_test_neighbor[_i][_qp] * _phi[_j][_qp];
      break;
    case Moose::ElementNeighbor:
      jac = _test[_i][_qp] *  _phi_neighbor[_j][_qp];
      break;
  }
  return jac;
}