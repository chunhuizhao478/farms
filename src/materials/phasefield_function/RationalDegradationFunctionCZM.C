//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "RationalDegradationFunctionCZM.h"

registerMooseObject("farmsApp", RationalDegradationFunctionCZM);

InputParameters
RationalDegradationFunctionCZM::validParams()
{
  InputParameters params = DegradationFunctionBase::validParams();
  params.addClassDescription(
      "Defines the rational degradation function $g(d) = "
      "(1-d)^p/((1-d)^p+a1*d*(1+a2*d+a2*a3*d^2))*(1-eta)+eta");

  params.set<std::string>("expression") =
      "(1-d)^p/((1-d)^p+a1*d*(1+a2*d+a2*a3*d^2))*(1-eta)+eta";
  params.set<std::vector<std::string>>("material_property_names") = {"a1","p", "a2", "a3", "eta"};
  return params;
}

RationalDegradationFunctionCZM::RationalDegradationFunctionCZM(const InputParameters & parameters)
  : DegradationFunctionBase(parameters)
{
}
