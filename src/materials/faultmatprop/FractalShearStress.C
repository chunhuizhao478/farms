//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FractalShearStress.h"
#include "DelimitedFileReader.h"

/**
 *  Created by Chunhui Zhao, Apr 12th, 2025
 *  Material used in Generate Fractal Shear Stress Distribution
 */
registerMooseObject("farmsApp", FractalShearStress);

InputParameters
FractalShearStress::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material to create fractal shear stress distribution");
  params.addRequiredParam<std::string>("csv_file", "Relative/Absolute Path to the CSV file containing the fractal shear stress data");
  return params;
}

FractalShearStress::FractalShearStress(const InputParameters & parameters)
  : Material(parameters),
  _fractal_shear_stress(declareProperty<Real>("fractal_shear_stress")),
  _fractal_shear_stress_old(getMaterialPropertyOldByName<Real>("fractal_shear_stress")),
  _csv_file(getParam<std::string>("csv_file"))
{
}

void
FractalShearStress::initQpStatefulProperties()
{
    // Get current quadrature point coordinates
    Real x_coord = _q_point[_qp](0);
    Real y_coord = _q_point[_qp](1);

    // Load CSV file
    MooseUtils::DelimitedFileReader reader(_csv_file);
    reader.read();

    // Get data vectors from CSV
    auto qp_coords_x             = reader.getData("x");            // vector<double>
    auto qp_coords_y             = reader.getData("y");            // vector<double>
    auto qp_coords_shear_stress  = reader.getData("shear_stress");   // vector<double>

    // Ensure we have data to process
    if (qp_coords_x.empty() || qp_coords_y.empty() || qp_coords_shear_stress.empty())
    {
        mooseError("CSV file did not load all required data (x, y, shear_stress)!");
    }

    // Number of data ptrs in the CSV file
    const unsigned num_ptrs = qp_coords_x.size();

    // Find the closest point in the CSV file using Euclidean distance
    Real min_dist = std::numeric_limits<Real>::max();
    unsigned closest_index = 0;
    for (unsigned i = 0; i < num_ptrs; i++)
    {
        Real dx = qp_coords_x[i] - x_coord;
        Real dy = qp_coords_y[i] - y_coord;
        Real dist = dx * dx + dy * dy;  // No need to compute sqrt for comparison
        if (dist < min_dist)
        {
            min_dist = dist;
            closest_index = i;
        }
    }

    // Assign the shear stress value from the closest point to the material property
    _fractal_shear_stress[_qp] = qp_coords_shear_stress[closest_index];
}

void
FractalShearStress::computeQpProperties()
{
  // If not at the first time step, use the previous value
  _fractal_shear_stress[_qp] = _fractal_shear_stress_old[_qp];
}