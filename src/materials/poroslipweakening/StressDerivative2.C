// from ElasticEnergyMaterial kernel

#include "StressDerivative2.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"


registerMooseObject("farmsApp", StressDerivative2);

InputParameters
StressDerivative2::validParams()
{
  InputParameters params = DerivativeFunctionMaterialBase::validParams();
  params.addClassDescription("Stress derivative.");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addCoupledVar("args", "Vector of variable arguments of the displacements");
  return params;
}

StressDerivative2::StressDerivative2(const InputParameters & parameters)
  : DerivativeFunctionMaterialBase(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _dstress(declarePropertyByName<RankTwoTensor>(_base_name + "dstress")),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")),
    _Jacobian_mult(getMaterialProperty<RankFourTensor>(_base_name + "Jacobian_mult")),
    _strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "elastic_strain"))
{
  _dstrain.resize(_nargs);
  _delasticity_tensor.resize(_nargs);

  // fetch stress and elasticity tensor derivatives (in simple eigenstrain models this is is only
  // w.r.t. 'displacements')
  for (unsigned int i = 0; i < _nargs; ++i)
  {
    _dstrain[i] = &getMaterialPropertyDerivativeByName<RankTwoTensor>(_base_name + "elastic_strain",
                                                                      _arg_names[i]);
    _delasticity_tensor[i] = &getMaterialPropertyDerivativeByName<RankFourTensor>(
        _base_name + "elasticity_tensor", _arg_names[i]);

  }
}

void
StressDerivative2::initialSetup()
{
  validateCoupling<RankTwoTensor>(_base_name + "elastic_strain");
  validateCoupling<RankFourTensor>(_base_name + "elasticity_tensor");
}


Real
StressDerivative2::computeDF(unsigned int i_var)
{
  unsigned int i = argIndex(i_var);

   _dstress[_qp] =  ((*_delasticity_tensor[i])[_qp] * _strain[_qp]) + (_elasticity_tensor[_qp] * (*_dstrain[i])[_qp]) ;
  // _dstress[_qp] =  _Jacobian_mult[_qp] * (*_dstrain[i])[_qp] ;

  // RankTwoTensor dstrain;

  // // Apply the strain-displacement derivatives for the 2D case
  // dstrain(0, 0) = (i_var == 0) ? 1 : 0.0;  // ∂ε_xx/∂u_x
  // dstrain(1, 1) = (i_var == 1) ? 1: 0.0;  // ∂ε_yy/∂u_y
  // dstrain(0, 1) = (i_var == 0) ? 0.5 : 0.0;  // ∂ε_xy/∂u_x
  // dstrain(1, 0) = (i_var == 1) ? 0.5 : 0.0;  // ∂ε_xy/∂u_y

  // // Compute the stress derivative
  // _dstress[_qp] = _Jacobian_mult[_qp] * dstrain;

  return 0.0;

}
