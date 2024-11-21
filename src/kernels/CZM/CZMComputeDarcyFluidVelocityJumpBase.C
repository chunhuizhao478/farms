
#include "CZMComputeDarcyFluidVelocityJumpBase.h"
#include "CohesiveZoneModelTools.h"

InputParameters
CZMComputeDarcyFluidVelocityJumpBase::validParams()
{
  InputParameters params = InterfaceMaterial::validParams();
  params.addClassDescription("Base class used to compute the fluid velocity jump across a czm "
                             "interface in local coordinates");
  params.addRequiredCoupledVar("fluid_vel",
                               "The string of fluid velocity suitable for the problem statement");
  params.suppressParameter<bool>("use_displaced_mesh");
  // MooseEnum interface_type("Permeable Impermeable Darcy Hydrofault", "Impermeable" );
  // params.addParam<MooseEnum>("interface_type", interface_type, "specify the fluid interface type");
  params.addParam<std::string>("base_name", "Material property base name");

  return params;
}

CZMComputeDarcyFluidVelocityJumpBase::CZMComputeDarcyFluidVelocityJumpBase(
    const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
                   ? getParam<std::string>("base_name") + "_"
                   : ""),
    _nvel(coupledComponents("fluid_vel")),
    _fluid_vel(3),
    _fluid_vel_neighbor(3),
    _fluid_vel_jump_global(declareGenericPropertyByName<RealVectorValue, false>(
        _base_name + "fluid_vel_jump_global")),
    _interface_fluid_vel_jump(declareGenericPropertyByName<RealVectorValue, false>(
        _base_name + "interface_fluid_vel_jump")),
    // _interface_type(getParam<MooseEnum>("interface_type").getEnum<InterfaceType>()),
    _czm_total_rotation(
        declareGenericPropertyByName<RankTwoTensor, false>(_base_name + "czm_total_rotation"))
{
  // Enforce consistency
  if (_nvel != _mesh.dimension())
    paramError("fluid_vel", "Number of fluid velocity must match problem dimension.");

  if (_nvel > 3 || _nvel < 1)
    mooseError("the CZM material requires 1, 2 or 3 fluid velocity variables");

    for (unsigned int i = 0; i < _nvel; ++i)
      {
        _fluid_vel[i] = &coupledGenericValue<false>("fluid_vel", i);
        _fluid_vel_neighbor[i] = &coupledGenericNeighborValue<false>("fluid_vel", i);
      }
      for (unsigned int i = _nvel; i < 3; i++)
      {
         _fluid_vel[i] = &_zero;
         _fluid_vel_neighbor[i] = &_zero;
      }

  // switch (_interface_type)
  // {
  //   case InterfaceType::Permeable:
  //   {
  //     for (unsigned int i = 0; i < _nvel; ++i)
  //     {
  //       _fluid_vel[i] = &coupledGenericValue<false>("fluid_vel", i);
  //       _fluid_vel_neighbor[i] = &coupledGenericNeighborValue<false>("fluid_vel", i);
  //     }

  //     for (unsigned int i = _nvel; i < 3; i++)
  //     {
  //         _fluid_vel[i] = &_zero;
  //         _fluid_vel_neighbor[i] = &_zero;
  //     }
  //     break;
  //   }
  //   case InterfaceType::Impermeable:
  //   {
  //     for (unsigned int i = 0; i < _nvel; ++i)
  //     {
  //       _fluid_vel[i] = &_zero;
  //       _fluid_vel_neighbor[i] = &_zero;
  //     }

  //     for (unsigned int i = _nvel; i < 3; i++)
  //     {
  //         _fluid_vel[i] = &_zero;
  //         _fluid_vel_neighbor[i] = &_zero;
  //     }
  //     break;
  //   }
  //   case InterfaceType::Hydrofault:
  //   {
  //     for (unsigned int i = 0; i < _nvel; ++i)
  //     {
  //       _fluid_vel[i] = &coupledGenericValue<false>("fluid_vel", i);
  //       _fluid_vel_neighbor[i] = &coupledGenericNeighborValue<false>("fluid_vel", i);
  //     }

  //     for (unsigned int i = _nvel; i < 3; i++)
  //     {
  //         _fluid_vel[i] = &_zero;
  //         _fluid_vel_neighbor[i] = &_zero;
  //     }
  //     break;
  //   }
  // }

}

void
CZMComputeDarcyFluidVelocityJumpBase::initQpStatefulProperties()
{
  // requried to promote _interface_fluid_vel_jump to stateful in case someone needs it
  _interface_fluid_vel_jump[_qp] = 0;
}

void
CZMComputeDarcyFluidVelocityJumpBase::computeQpProperties()
{

  // computing the fluid velocity  jump
  for (unsigned int i = 0; i < _nvel; i++)
    _fluid_vel_jump_global[_qp](i) = (*_fluid_vel_neighbor[i])[_qp] - (*_fluid_vel[i])[_qp];
  for (unsigned int i = _nvel; i < 3; i++)
    _fluid_vel_jump_global[_qp](i) = 0;

  computeRotationMatrices();
  computeLocalFluidVelocityJump();
}

void
CZMComputeDarcyFluidVelocityJumpBase::computeRotationMatrices()
{
  _czm_total_rotation[_qp] =
      CohesiveZoneModelTools::computeReferenceRotation(_normals[_qp], _mesh.dimension());
}

