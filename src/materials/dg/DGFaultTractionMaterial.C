//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DGFaultTractionMaterial.h"

registerMooseObject("farmsApp", DGFaultTractionMaterial);

InputParameters
DGFaultTractionMaterial::validParams()
{
  InputParameters params = InterfaceMaterial::validParams();
  params.addClassDescription(
      "Computes fault traction for DG SEAS simulations. Provides prescribed or "
      "linear slip-dependent traction for testing. For full SEAS, use rate-state friction.");

  params.addRequiredCoupledVar("displacement", "The displacement variable");

  params.addParam<Real>("tau0", 0.0, "Background/initial shear traction");
  params.addParam<Real>("dtau_ds", 0.0,
                        "Derivative of traction w.r.t. slip (positive = slip-strengthening, "
                        "negative = slip-weakening)");

  params.addParam<MaterialPropertyName>("fault_traction_name", "fault_traction",
                                        "Name of the fault traction material property");
  params.addParam<MaterialPropertyName>("dtraction_dslip_name", "dtraction_dslip",
                                        "Name of the traction derivative material property");

  return params;
}

DGFaultTractionMaterial::DGFaultTractionMaterial(const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _u(coupledValue("displacement")),
    _u_neighbor(coupledNeighborValue("displacement")),
    _tau0(getParam<Real>("tau0")),
    _dtau_ds(getParam<Real>("dtau_ds")),
    _fault_traction(
        declarePropertyByName<Real>(getParam<MaterialPropertyName>("fault_traction_name"))),
    _dtraction_dslip(
        declarePropertyByName<Real>(getParam<MaterialPropertyName>("dtraction_dslip_name")))
{
}

void
DGFaultTractionMaterial::computeQpProperties()
{
  // Compute slip as jump in displacement
  // Convention: slip = u_elem - u_neighbor (positive when element side moves in +z)
  const Real slip = _u[_qp] - _u_neighbor[_qp];

  // Compute traction using linear model: τ = τ0 + (dτ/ds) * s
  // For testing:
  //   - dtau_ds = 0: constant traction
  //   - dtau_ds > 0: slip-strengthening (stable)
  //   - dtau_ds < 0: slip-weakening (unstable, needs care)
  _fault_traction[_qp] = _tau0 + _dtau_ds * slip;

  // Derivative for Jacobian
  _dtraction_dslip[_qp] = _dtau_ds;
}
