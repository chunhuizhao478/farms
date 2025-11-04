//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElkRadialAverageUpdated.h"

#include "libmesh/nanoflann.hpp"

using QPDataRangeUpdated =
    StoredRange<std::vector<ElkRadialAverageUpdated::QPData>::const_iterator, ElkRadialAverageUpdated::QPData>;

/**
 * ElkRadialAverageUpdated threaded loop
 */
class ThreadedElkRadialAverageLoopUpdated
{
public:
  ThreadedElkRadialAverageLoopUpdated(ElkRadialAverageUpdated &);

  /// Splitting constructor
  ThreadedElkRadialAverageLoopUpdated(const ThreadedElkRadialAverageLoopUpdated & x, Threads::split split);

  /// dummy virtual destructor
  virtual ~ThreadedElkRadialAverageLoopUpdated() {}

  /// parens operator with the code that is executed in threads
  void operator()(const QPDataRangeUpdated & range);

  /// thread join method
  virtual void join(const ThreadedElkRadialAverageLoopUpdated & /*x*/) {}

protected:
  /// rasterizer to manage the sample data
  ElkRadialAverageUpdated & _radavg;

  /// ID number of the current thread
  THREAD_ID _tid;

private:
};
