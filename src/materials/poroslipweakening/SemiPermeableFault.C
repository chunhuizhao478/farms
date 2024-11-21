/* 
Material Description of Slip Weakening Friction v3
Generalize the computation of sticking traction using consistent displacement jump and nodal reaction forces
*/

#include "SemiPermeableFault.h"
#include "InterfaceKernel.h"
#include "FEProblem.h"

registerMooseObject("farmsApp", SemiPermeableFault);

InputParameters
SemiPermeableFault::validParams()
{
  InputParameters params = InterfaceMaterial::validParams();
  params.addClassDescription("Leaky fault model song and rudnicki 2017.");
  params.addCoupledVar("pressure_secondary","Pressure at postive side of the fault");
  params.addCoupledVar("pressure_main","Pressure at minus side of the fault");
  params.addParam<std::string>("base_name", "Material property base name");
  return params;
}

SemiPermeableFault::SemiPermeableFault(const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
                   ? getParam<std::string>("base_name") + "_"
                   : ""),
    _Transmissibility(getMaterialPropertyByName<Real>(_base_name + "Transmissibility")),
    _interface_pressure_secondary(coupledNeighborValue("pressure_secondary")),
    _interface_pressure_main(coupledValue("pressure_main")),
    _across_fault_flux(declarePropertyByName<Real>(_base_name + "across_fault_flux"))
    
{
}

void
SemiPermeableFault::computeQpProperties()
{   

  // if(_fe_problem.getCurrentExecuteOnFlag()=="LINEAR")
  //  {

  // Semi Permeable Fault condition
   _across_fault_flux[_qp] = - _Transmissibility[_qp] * (_interface_pressure_main[_qp] - _interface_pressure_secondary[_qp]) ;

  // std::cout << _interface_pressure_main[_qp] << std::endl;


  // }
}