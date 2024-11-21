#include "ViscoelasticStressKernel.h"
#include "ElasticityTensorTools.h"

registerMooseObject("farmsApp", ViscoelasticStressKernel);

InputParameters
ViscoelasticStressKernel::validParams()
{
  InputParameters params = StressDivergenceTensors::validParams();
  params.addClassDescription("Compute Viscoelastic Stress Residual using eta_kv/G * stress_rate");
  params.addRequiredParam<Real>("kelvin_voigt_viscosity", "Kelvin-Voigt viscosity parameter (eta_kv)");
  params.addRequiredParam<Real>("shear_modulus", "Shear modulus (G)");
  return params;
}

ViscoelasticStressKernel::ViscoelasticStressKernel(const InputParameters & parameters)
  : StressDivergenceTensors(parameters),
    _stress_older(getMaterialPropertyOlderByName<RankTwoTensor>(_base_name + "stress")),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _eta_kv(getParam<Real>("kelvin_voigt_viscosity")),
    _G(getParam<Real>("shear_modulus")),
    _dim(_mesh.dimension())
{
}

Real
ViscoelasticStressKernel::computeQpResidual()
{
  // Get current and old stress
  RankTwoTensor current_stress = _stress[_qp];
  RankTwoTensor old_stress = _stress_older[_qp];

  // Calculate isochoric parts
  Real current_vol = current_stress.trace() / (_dim == 2 ? 2.0 : 3.0);
  Real old_vol = old_stress.trace() / (_dim == 2 ? 2.0 : 3.0);

  RankTwoTensor current_iso = current_stress;
  RankTwoTensor old_iso = old_stress;

  // Remove volumetric components
  for (unsigned int i = 0; i < _dim; ++i)
  {
    current_iso(i,i) -= current_vol;
    old_iso(i,i) -= old_vol;
  }

  // Calculate viscous stress: η_kv/G * stress_rate
  RankTwoTensor viscous_stress = (current_iso - old_iso);
  viscous_stress *= _eta_kv / (_G * _dt);

  return viscous_stress.row(_component) * _grad_test[_i][_qp];
}

Real
ViscoelasticStressKernel::computeQpJacobian()
{
  // Get the full elasticity jacobian from StressDivergenceTensors
  Real jac = StressDivergenceTensors::computeQpJacobian();

  // Get the Jacobian multiplier tensor which contains elasticity tensor
  const RankFourTensor & C = _Jacobian_mult[_qp];

  // Calculate the volumetric part to be removed
  Real vol_jac = 0.0;
  for (unsigned int i = 0; i < _dim; ++i)
    vol_jac += C(i,i,_component,_component) * _grad_phi[_j][_qp](i);
  vol_jac /= (_dim == 2 ? 2.0 : 3.0);

  // Remove volumetric component
  jac -= vol_jac * _grad_test[_i][_qp](_component);

  // Scale by viscosity/shear ratio and time
  jac *= _eta_kv / (_G * _dt);

  return jac;
}

Real
ViscoelasticStressKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    if (jvar == _disp_var[i])
    {
      // Get the full elasticity off-diagonal jacobian
      Real jac = StressDivergenceTensors::computeQpOffDiagJacobian(jvar);

      // Get the Jacobian multiplier tensor
      const RankFourTensor & C = _Jacobian_mult[_qp];

      // Calculate volumetric part
      Real vol_jac = 0.0;
      for (unsigned int k = 0; k < _dim; ++k)
        vol_jac += C(k,k,_component,i) * _grad_phi[_j][_qp](i);
      vol_jac /= (_dim == 2 ? 2.0 : 3.0);

      // Remove volumetric component
      jac -= vol_jac * _grad_test[_i][_qp](_component);

      // Scale by viscosity/shear ratio and time
      jac *= _eta_kv / (_G * _dt);

      return jac;
    }
  }
  return 0.0;
}