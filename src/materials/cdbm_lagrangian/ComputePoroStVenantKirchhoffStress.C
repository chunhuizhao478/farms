//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputePoroStVenantKirchhoffStress.h"

registerMooseObject("farmsApp", ComputePoroStVenantKirchhoffStress);

InputParameters
ComputePoroStVenantKirchhoffStress::validParams()
{
  InputParameters params = ComputeLagrangianStressPK1::validParams();
  params.addParam<MaterialPropertyName>(
    "elasticity_tensor", "elasticity_tensor", "The name of the elasticity tensor.");
     params.addRequiredCoupledVar(
    "porepressure", 
    "The pore pressure appropriate for the simulation geometry and coordinate system");
  
  return params;
}

ComputePoroStVenantKirchhoffStress::ComputePoroStVenantKirchhoffStress(const InputParameters & parameters)
  : ComputeLagrangianStressPK1(parameters),
    GuaranteeConsumer(this),
    _elasticity_tensor_name(getParam<MaterialPropertyName>("elasticity_tensor")),
    _elasticity_tensor(getMaterialProperty<RankFourTensor>(_elasticity_tensor_name)),
    _pore_pressure(coupledValue("porepressure")),
    _E(declareProperty<RankTwoTensor>(_base_name + "green_lagrange_strain")),
    _S(declareProperty<RankTwoTensor>(_base_name + "pk2_stress")),  
    _C(declareProperty<RankFourTensor>(_base_name + "pk2_jacobian")),
    _pk1_off_diag_jacobian(declareProperty<RankTwoTensor>(_base_name + "pk1_off_diag_jacobian")),
    _I1(declareProperty<Real>(_base_name + "first_elastic_strain_invariant")),
    _I1_small(declareProperty<Real>(_base_name + "I1")), 
    _dI1dF(declareProperty<RankTwoTensor>(_base_name + "first_elastic_strain_invariant_derivative")), 
    _biot_coeff_eff(getMaterialProperty<Real>(_base_name + "biot_coefficient_effective"))
{
}

void
ComputePoroStVenantKirchhoffStress::initialSetup()
{
  // Enforce isotropic elastic tensor
  if (!hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
    mooseError("ComputePoroStVenantKirchhoffStress requires an isotropic elasticity tensor");
}

void
ComputePoroStVenantKirchhoffStress::initQpStatefulProperties()
{
  _S[_qp].zero();
  _C[_qp].zero();
  _dI1dF[_qp].zero();
  _pk1_off_diag_jacobian[_qp].zero();
}

void
ComputePoroStVenantKirchhoffStress::computeQpPK2Stress()
{
  // Hyperelasticity is weird, we need to branch on the type of update if we
  // want a truly linear model
  //
  // This is because we need to drop quadratic terms for the linear update to
  // use a linear strain measure
  // Jacobian is the same for both the small and Green-Lagrange strains
  _C[_qp] = _elasticity_tensor[_qp];

  // Get the right strain
  RankTwoTensor strain;
  if (_large_kinematics) // Large deformations = Green-Lagrange strain
    strain = _E[_qp];
  else // Small deformations = linear strain
    strain = 0.5 * (_F[_qp] + _F[_qp].transpose()) - RankTwoTensor::Identity();

  _I1[_qp] = strain.trace();
  _I1_small[_qp] = _I1[_qp];

  // The stress update is linear with the correct strains/frame
  _S[_qp] = _C[_qp] * strain - _biot_coeff_eff[_qp] * _pore_pressure[_qp] * RankTwoTensor::Identity();
}

void
ComputePoroStVenantKirchhoffStress::computeQpPK1Stress()
{
  // Calculate the green-lagrange strain for the benefit of the subclasses
  _E[_qp] = 0.5 * (_F[_qp].transpose() * _F[_qp] - RankTwoTensor::Identity());

  // PK2 update
  computeQpPK2Stress();

  // Complicated wrapping, see documentation
  if (_large_kinematics)
  {
    _pk1_stress[_qp] = _F[_qp] * _S[_qp];
    
    usingTensorIndices(i_, j_, k_, l_);
    RankFourTensor dE =
      0.5 * (RankTwoTensor::Identity().times<i_, l_, j_, k_>(_F[_qp].transpose()) +
             _F[_qp].transpose().times<i_, k_, j_, l_>(RankTwoTensor::Identity()));
    
    _pk1_jacobian[_qp] = RankTwoTensor::Identity().times<i_, k_, j_, l_>(_S[_qp].transpose()) +
                         (_C[_qp] * dE).singleProductI(_F[_qp]);

    // dI1/dF = d(tr(E))/dF = d(tr(0.5*(F^T*F - I)))/dF = F
    _dI1dF[_qp] = _F[_qp];

    // Compute off-diagonal Jacobian: ∂P/∂p = ∂(F*S)/∂p = F * (∂S/∂p)
    // From PK2 stress: S = C*E - α*p*I, so ∂S/∂p = -α*I
    RankTwoTensor dS_dp = -_biot_coeff_eff[_qp] * RankTwoTensor::Identity();

    _pk1_off_diag_jacobian[_qp] = _F[_qp] * dS_dp;
  }
  else
  {
    _pk1_stress[_qp] = _S[_qp];
    _pk1_jacobian[_qp] = _C[_qp];
    
    // For small deformations: P = S, so ∂P/∂p = ∂S/∂p = -α*I
    _pk1_off_diag_jacobian[_qp] = -_biot_coeff_eff[_qp] * RankTwoTensor::Identity();
    
    // For small deformations: E ≈ 0.5*(F + F^T) - I, so dI1/dF = I (identity tensor)
    _dI1dF[_qp] = RankTwoTensor::Identity();
  }
}
