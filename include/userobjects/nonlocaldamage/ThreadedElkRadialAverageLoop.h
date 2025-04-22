//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElkRadialAverage.h"

#include "libmesh/nanoflann.hpp"

using QPDataRange =
    StoredRange<std::vector<ElkRadialAverage::QPData>::const_iterator, ElkRadialAverage::QPData>;

/**
 * ElkRadialAverage threaded loop
 */
class ThreadedElkRadialAverageLoop
{
public:
  ThreadedElkRadialAverageLoop(ElkRadialAverage &);

  /// Splitting constructor
  ThreadedElkRadialAverageLoop(const ThreadedElkRadialAverageLoop & x, Threads::split split);

  /// dummy virtual destructor
  virtual ~ThreadedElkRadialAverageLoop() {}

  /// parens operator with the code that is executed in threads
  void operator()(const QPDataRange & range);

  /// thread join method
  virtual void join(const ThreadedElkRadialAverageLoop & /*x*/) {}

protected:
  /// rasterizer to manage the sample data
  ElkRadialAverage & _radavg;

  /// ID number of the current thread
  THREAD_ID _tid;

private:
};