//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// Modern MOOSE main using moose::main<TApp>.
// This avoids printing anything to stdout before the Language Server
// establishes the LSP connection, and it handles initialization/teardown.

#include "MooseMain.h"
#include "farmsApp.h"   // Your app class header (defines farmsTestApp)

int main(int argc, char * argv[])
{
  // moose::main<T>() performs MPI/init, registers apps, builds the MooseApp, and runs it.
  return Moose::main<farmsApp>(argc, argv);
}
