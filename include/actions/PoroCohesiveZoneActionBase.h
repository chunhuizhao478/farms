
#pragma once

#include "Action.h"

class PoroCohesiveZoneActionBase : public Action
{
public:
  static InputParameters validParams();

  PoroCohesiveZoneActionBase(const InputParameters & params);

  ///@{ output methods
  static MultiMooseEnum outputPropertiesType();
  static MultiMooseEnum materialOutputOrders();
  static MultiMooseEnum materialOutputFamilies();
  ///@}

  ///@{ table data for output generation
  static const std::map<std::string, std::string> _real_vector_cartesian_component_table;
  static const std::map<std::string, std::string> _real_pressure;
  static const std::map<std::string, std::pair<std::string, std::vector<std::string>>>
      _vector_direction_table;
  static const std::vector<char> _component_table;
  ///@}
};