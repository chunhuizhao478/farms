

#include "ComputeFluidVelocityGradient.h"
#include "Assembly.h"
#include "libmesh/quadrature.h"

registerMooseObject("farmsApp", ComputeFluidVelocityGradient);

InputParameters
ComputeFluidVelocityGradient::validParams()
{
  InputParameters params = DarcyFluidVelocity::validParams();
  params.addClassDescription("Compute a darcy velocity divergence.");
  return params;
}

ComputeFluidVelocityGradient::ComputeFluidVelocityGradient(const InputParameters & parameters)
  : DarcyFluidVelocity(parameters)
{
}

void
ComputeFluidVelocityGradient::computeProperties()
{
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    
    const auto grad_tensor = RankTwoTensor ::initializeFromRows(
        (*_grad_vel[0])[_qp], (*_grad_vel[1])[_qp], (*_grad_vel[2])[_qp]);

        _darcy_vel_grad [_qp] = grad_tensor;
  }
}

