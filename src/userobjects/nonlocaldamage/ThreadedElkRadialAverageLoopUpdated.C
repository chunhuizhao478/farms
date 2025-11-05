//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ThreadedElkRadialAverageLoopUpdated.h"
#include "Function.h"

// Make newer nanoflann API spelling compatible with older nanoflann
// versions
#if NANOFLANN_VERSION < 0x150
namespace nanoflann
{
typedef SearchParams SearchParameters;

template <typename T, typename U>
using ResultItem = std::pair<T, U>;
}
#endif

ThreadedElkRadialAverageLoopUpdated::ThreadedElkRadialAverageLoopUpdated(ElkRadialAverageUpdated & green) : _radavg(green) {}

// Splitting Constructor
ThreadedElkRadialAverageLoopUpdated::ThreadedElkRadialAverageLoopUpdated(const ThreadedElkRadialAverageLoopUpdated & x,
                                                     Threads::split /*split*/)
  : _radavg(x._radavg)
{
}

void
ThreadedElkRadialAverageLoopUpdated::operator()(const QPDataRangeUpdated & qpdata_range)
{
  // fetch data from parent
  const auto radius = _radavg._radius;
  const auto & qp_data = _radavg._qp_data;
  const auto & kd_tree = _radavg._kd_tree;
  const auto & weights_type = _radavg._weights_type;

  /*----------------------------------------------*/
  // get characteristic length scale
  const auto length_scale = _radavg._l;
  /*----------------------------------------------*/

  // tree search data structures
  std::vector<nanoflann::ResultItem<std::size_t, Real>> ret_matches;
  nanoflann::SearchParameters search_params;

  // result map entry
  const auto end_it = _radavg._average.end();
  auto it = end_it;

  // iterate over qp range
  for (auto && local_qp : qpdata_range)
  {
    // Look up result map iterator only if we enter a new element. this saves a bunch
    // of map lookups because same element entries are consecutive in the qp_data vector.
    if (it == end_it || it->first != local_qp._elem_id)
      it = _radavg._average.find(local_qp._elem_id);

    // initialize result entry
    mooseAssert(it != end_it, "Current element id not found in result set.");
    auto & sum = it->second[local_qp._qp];
    sum = 0.0;

    ret_matches.clear();
    std::size_t n_result =
        kd_tree->radiusSearch(&(local_qp._q_point(0)), radius * radius, ret_matches, search_params);
    Real total_vol = 0.0;
    Real weight = 1.0;
    /*----------------------------------------------*/
    //Real k = std::pow(6*sqrt(libMesh::pi), 1.0/3.0);
    /*----------------------------------------------*/
    for (std::size_t j = 0; j < n_result; ++j)
    {
      const auto & other_qp = qp_data[ret_matches[j].first];
      const Real distance_sq = ret_matches[j].second;
      const Real distance = std::sqrt(distance_sq);

      switch (weights_type)
      {
        case ElkRadialAverageUpdated::WeightsType::CONSTANT:
          break;

        case ElkRadialAverageUpdated::WeightsType::LINEAR:
          weight = radius - distance;
          break;

        case ElkRadialAverageUpdated::WeightsType::COSINE:
          weight = std::cos(distance / radius * libMesh::pi) + 1.0;
          break;
        /*-------------------------------------------------------------------------------------------------*/
        case ElkRadialAverageUpdated::WeightsType::BAZANT:
        {
          // Bažant nonlocal kernel (Eq. (4) in Bažant & Pijaudier-Cabot, 1988): exp[-(2 r / l)^2]
          const Real inv_l2 = 1.0 / (length_scale * length_scale);
          weight = std::exp(-4.0 * distance_sq * inv_l2);
          break;
        }
       /*--------------------------------------------------------------------------------------------------*/
        case ElkRadialAverageUpdated::WeightsType::BAZANT3D:
        {
          // 3D Bažant kernel (Eq. (5) in Bažant & Pijaudier-Cabot, 1988): exp[-(k r / l)^2] with k^2 = (6 sqrt(pi))^{2/3}
          static const Real bazant3d_coeff = std::pow(6.0 * std::sqrt(libMesh::pi), 2.0 / 3.0);
          const Real inv_l2 = 1.0 / (length_scale * length_scale);
          weight = std::exp(-bazant3d_coeff * distance_sq * inv_l2);
          break;
        }
       /*--------------------------------------------------------------------------------------------------*/
      }

      sum += other_qp._value * other_qp._volume * weight;
      total_vol += other_qp._volume * weight;
    }
    sum /= total_vol;
  }
}
