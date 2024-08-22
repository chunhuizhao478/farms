//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeLagrangianDamageBreakageStressPK2.h"


registerMooseObject("farmsApp", ComputeLagrangianDamageBreakageStressPK2);

InputParameters
ComputeLagrangianDamageBreakageStressPK2::validParams()
{
  InputParameters params = ComputeLagrangianStressPK1::validParams();
  return params;
}

ComputeLagrangianDamageBreakageStressPK2::ComputeLagrangianDamageBreakageStressPK2(const InputParameters & parameters)
  : ComputeLagrangianStressPK1(parameters),
  _Fp(declareProperty<RankTwoTensor>(_base_name + "plastic_deformation_gradient")),
  _Tau(declareProperty<RankTwoTensor>(_base_name + "deviatroic_stress")),
  _Ee(declareProperty<RankTwoTensor>(_base_name + "green_lagrange_elastic_strain")),
  _I1(declareProperty<Real>(_base_name + "first_elastic_strain_invariant")),
  _I2(declareProperty<Real>(_base_name + "second_elastic_strain_invariant")),
  _xi(declareProperty<Real>(_base_name + "strain_invariant_ratio")),
  _S(declareProperty<RankTwoTensor>(_base_name + "pk2_stress")),
  _C(declareProperty<RankFourTensor>(_base_name + "pk2_jacobian")),
  _shear_modulus(getMaterialProperty<Real>("shear_modulus")),
  _damaged_modulus(getMaterialProperty<Real>("damaged_modulus")),
  _B_breakagevar(getMaterialProperty<Real>("B_damagedvar")),
  _B_breakagevar_old(getMaterialPropertyOldByName<Real>("B_damagedvar")),
  _Tau_old(getMaterialPropertyOldByName<RankTwoTensor>("deviatroic_stress")),
  _Fp_old(getMaterialPropertyOldByName<RankTwoTensor>("plastic_deformation_gradient")),
  _lambda(getParam<Real>("lambda")),
  _C_g(getParam<Real>("C_g")),
  _m1(getParam<Real>("m1")),
  _m2(getParam<Real>("m2")),
  _a0(getParam<Real>("a0")),
  _a1(getParam<Real>("a1")),
  _a2(getParam<Real>("a2")),
  _a3(getParam<Real>("a3")),
  _dim(_mesh.dimension())
{
}

void
ComputeLagrangianDamageBreakageStressPK2::computeQpPK1Stress()
{
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
  }
  // Small deformations all are equivalent
  else
  {
    _pk1_stress[_qp] = _S[_qp];
    _pk1_jacobian[_qp] = _C[_qp];
  }
}

void
ComputeLagrangianDamageBreakageStressPK2::computeQpPK2Stress()
{
  /* Evaluate Fp */
  RankTwoTensor Fp_updated = computeQpFp();

  /* Compute Fe */
  RankTwoTensor Fe = _F[_qp] * Fp_updated.inverse();

  /* Compute Ee */
  RankTwoTensor Ee = 0.5 * (Fe.transpose() * Fe - RankTwoTensor::Identity());

  /* Compute I1 */
  Real I1 = Ee.trace();

  /* Compute I2 */
  Real I2 = 0.0;
  for (unsigned int i = 0; i < 3; ++i){
    for (unsigned int j = 0; j < 3; ++j){
      I2 += Ee(i,j) * Ee(i,j);
    }
  }

  /* Compute xi */
  Real xi = I1 / std::sqrt(I2);

  /* Compute stress */
  RankTwoTensor sigma_s = (_lambda - _damaged_modulus[_qp] / xi) * I1 * RankTwoTensor::Identity() + (2 * _shear_modulus[_qp] - _damaged_modulus[_qp] * xi) * Ee;
  RankTwoTensor sigma_b = (2 * _a2 + _a1 / xi + 3 * _a3 * xi) * I1 * RankTwoTensor::Identity() + (2 * _a0 + _a1 * xi - _a3 * std::pow(xi, 3)) * Ee;
  RankTwoTensor sigma_total = (1 - _B_breakagevar[_qp]) * sigma_s + _B_breakagevar[_qp] * sigma_b;

  /* Compute tangent */
  std::vector<Real> components(21);
  RankFourTensor tangent(components, RankFourTensor::symmetric21);
}

RankTwoTensor
ComputeLagrangianDamageBreakageStressPK2::computeQpFp()
{
  //Apply power operation on every element of Tau
  RankTwoTensor Tau_old_power2;
  for (unsigned int i = 0; i < _dim; i++){
    for (unsigned int j = 0; j < _dim; j++){
      Tau_old_power2(i,j) = std::pow(_Tau_old[_qp](i,j), _m2);
    }
  }

  //Compute Plastic Deformation Rate Tensor Dp at t_{n+1} using quantities from t_{n}
  RankTwoTensor Dp = _C_g * std::pow(_B_breakagevar_old[_qp], _m1) * Tau_old_power2; 

  //Compute Cp = I - Dp dt
  RankTwoTensor Cp = RankTwoTensor::Identity() - Dp * _dt;

  //Use Implicit Euler Integration, Update Fp
  RankTwoTensor Fp_updated = Cp.inverse() * _Fp_old[_qp];

  //Update material property
  _Fp[_qp] = Fp_updated;

  return Fp_updated;
}