//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "FarmsExternalWork.h"

registerMooseObject("farmsApp", FarmsExternalWork);

InputParameters
FarmsExternalWork::validParams()
{
  InputParameters params = NodalPostprocessor::validParams();
  params.addClassDescription("This class computes the total external work. The power expenditure "
                             "(rate of external work) is defined as $\\mathcal{P}^\\text{ext} = "
                             "\\int_\\bodyboundary \\bft \\cdot \\dot{\\bs{\\phi}} \\diff{A}$. The "
                             "power expenditure is integrated in time to get the total work.");
  params.addRequiredCoupledVar("forces",
                               "The reaction forces associated with each of the displacement");
  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system");
  return params;
}

FarmsExternalWork::FarmsExternalWork(const InputParameters & parameters)
  : NodalPostprocessor(parameters),
    _sum(0),
    _ndisp(coupledComponents("displacements")),
    _u_dots(coupledDots("displacements")),
    _nforce(coupledComponents("forces")),
    _forces(coupledValues("forces")),
    _sum_old(getPostprocessorValueOldByName(name()))
{
  for (unsigned int i = _ndisp; i < 3; ++i)
    _u_dots.push_back(&_zero);

  for (unsigned int i = _nforce; i < 3; ++i)
    _forces.push_back(&_zero);
}

void
FarmsExternalWork::initialize()
{
  _sum = 0;
}

void
FarmsExternalWork::execute()
{
  _sum += computeQpValue();
}

Real
FarmsExternalWork::computeQpValue()
{
  RealVectorValue u_dot((*_u_dots[0])[_qp], (*_u_dots[1])[_qp], (*_u_dots[2])[_qp]);
  RealVectorValue force((*_forces[0])[_qp], (*_forces[1])[_qp], (*_forces[2])[_qp]);

  return force * u_dot;
}

Real
FarmsExternalWork::getValue() const
{
  return _sum * _dt + _sum_old;
}

void
FarmsExternalWork::finalize()
{
  gatherSum(_sum);
}

void
FarmsExternalWork::threadJoin(const UserObject & y)
{
  const FarmsExternalWork & pps = static_cast<const FarmsExternalWork &>(y);
  _sum += pps._sum;
}
