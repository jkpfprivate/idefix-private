// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************


#include <unistd.h>

#include <sys/time.h>
#include <stdlib.h>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

//#include <petsc.h>
//#include <slepc.h>

#include "idefix.hpp"
#include "menhir.hpp"
#include "newton.hpp"
#include "profiler.hpp"
#include "input.hpp"
#include "units.hpp"
#include "grid.hpp"
#include "gridHost.hpp"
#include "fluid.hpp"
#include "dataBlock.hpp"
#include "timeIntegrator.hpp"
#include "setup.hpp"
#include "output.hpp"
#ifdef WITH_MPI
#include "mpi.hpp"
#endif

Menhir::Menhir(int argc, char* argv[]) {
  Initialise(argc, argv);
}

void Menhir::Initialise(int argc, char* argv[]) {
    idfx::initialize();
    ///////////////////////////////
    // Idefix Initialization
    ///////////////////////////////

    input = std::make_unique<Input>(argc, argv);
    input->PrintLogo();
    idfx::cout << "Main: initialization stage." << std::endl;

    haveDNS = input->haveDNS;
    haveNewton = input->haveNewton;
    haveStability = input->haveStability;
    haveContinuation = input->haveContinuation;

    // Init the units when needed
    idfx::units.Init(*input);

    // Allocate the grid on device
    grid = std::make_unique<Grid>(*input);
    // Allocate the grid image on host
    gridHost = std::make_unique<GridHost>(*grid);

    // Actually make the grid on host and sync it on the device
    gridHost->MakeGrid(*input);
    gridHost->SyncToDevice();

    // instantiate required objects.

    data = std::make_unique<DataBlock>(*grid, *input);
    Tint = std::make_unique<TimeIntegrator>(*input, *data);
    #ifdef WITH_PYTHON
      pydefix = std::make_unique<Pydefix>(*input);
    #endif
    output = std::make_unique<Output>(*input, *data);
    mysetup = std::make_unique<Setup>(*input, *grid, *data, *output);

    idfx::cout << "Main: initialisation finished." << std::endl;

    char host[1024];
    gethostname(host,1024);

    idfx::cout << "Main: running on " << std::string(host) << std::endl;

    ///////////////////////////////
    // Show configuration
    ///////////////////////////////
    if(initKokkosBeforeMPI) {
      idfx::cout << "Main: detected your configuration needed Kokkos to be initialised before MPI. "
                 << std::endl;
    }
    input->ShowConfig();
    idfx::units.ShowConfig();
    grid->ShowConfig();
    data->ShowConfig();
    Tint->ShowConfig();
    #ifdef WITH_PYTHON
    pydefix->ShowConfig();
    #endif

    ///////////////////////////////
    // Initial conditions (or restart)
    ///////////////////////////////
    // Are we restarting?
    if(input->restartRequested) {
      if(input->forceInitRequested) {
        #ifdef WITH_PYTHON
          if(pydefix->haveInitflow) {
            idfx::pushRegion("Pydefix::Initflow");
            pydefix->InitFlow(*data);
          } else {
            idfx::pushRegion("Setup::Initflow");
            mysetup->InitFlow(*data);
          }
          data->DeriveVectorPotential();
          idfx::popRegion();
        #else
          idfx::pushRegion("Setup::Initflow");
          mysetup->InitFlow(*data);
          data->DeriveVectorPotential();
          idfx::popRegion();
        #endif
      }
      idfx::cout << "Main: Restarting from dump file."  << std::endl;
      bool restartSuccess = output->RestartFromDump(*data,input->restartFileNumber);
      if(!restartSuccess) {
        idfx::cout << "Main: restart aborted." << std::endl;
        input->restartRequested = false;
      } else {
        data->SetBoundaries();
      }
    }
    if(!input->restartRequested) {
      idfx::cout << "Main: Creating initial conditions." << std::endl;
      #ifdef WITH_PYTHON
        if(pydefix->haveInitflow) {
          idfx::pushRegion("Pydefix::Initflow");
          pydefix->InitFlow(*data);
        } else {
          idfx::pushRegion("Setup::Initflow");
          mysetup->InitFlow(*data);
        }
      #else
        idfx::pushRegion("Setup::Initflow");
        mysetup->InitFlow(*data);
      #endif
      idfx::popRegion();
      data->DeriveVectorPotential();   // This does something only when evolveVectorPotential is on
      data->SetBoundaries();
      data->Validate();
      output->CheckForWrites(*data);
    }
}

void Menhir::PerformDNS(real stopping_time) {
    ///////////////////////////////
    // Main Loop
    ///////////////////////////////
    idfx::cout << "Main: Cycling Time Integrator..." << std::endl;

    output->ResetTimer();

    int tstop;
    if (stopping_time < 0.) tstop = input->Get<real>("TimeIntegrator","tstop",0);
    else tstop = stopping_time;

    while(data->t < tstop) {
      if(tstop-data->t < data->dt) data->dt = tstop-data->t;
      try {
        Tint->Cycle(*data);
      } catch(std::exception &e) {
        idfx::cout << "Main: WARNING! Caught an exception in TimeIntegrator." << std::endl;
        #ifdef WITH_MPI
          if(!Mpi::CheckSync(5)) {
            std::stringstream message;
            message << "A non-synchronous exception was raised in TimeIntegrator:" << std::endl;
            message << e.what();
            message << std::endl << "No emergency output can be produced." << std::endl;
            IDEFIX_ERROR(message);
          }
        #endif
        idfx::cout << e.what() << std::endl;
        idfx::cout << "Main: attempting to save the current state for inspection." << std::endl;
        output->ForceWriteVtk(*data);
        idfx::cout << "Main: Aborting current calculation." << std::endl;
        returnCode = 1;
        break;
      }
      output->CheckForWrites(*data);
      if(input->CheckForAbort() || Tint->CheckForMaxRuntime() ) {
        idfx::cout << "Main: Saving current state and aborting calculation." << std::endl;
        output->ForceWriteDump(*data);
        returnCode = -1;
        break;
      }
      if(input->maxCycles>=0) {
        if(Tint->GetNCycles() >= input->maxCycles) {
          idfx::cout << "Main: Reached maximum number of integration cycles." << std::endl;
          break;
        }
      }
    }

    int n_days{0}, n_hours{0}, n_minutes{0}, n_seconds{0};
    div_t divres;
    divres = div(timer.seconds(), 86400);
    n_days = divres.quot;
    divres = div(divres.rem, 3600);
    n_hours = divres.quot;
    divres = div(divres.rem, 60);
    n_minutes = divres.quot;
    n_seconds = divres.rem;

    double perfs = timer.seconds() / grid->np_int[IDIR] / grid->np_int[JDIR]
                            / grid->np_int[KDIR] / Tint->GetNCycles() * idfx::psize;

    idfx::cout << "Main: Reached t=" << data->t << std::endl;
    idfx::cout << "Main: Completed in ";
    if (n_days > 0) {
      idfx::cout << n_days << " day";
      if (n_days != 1) {
        idfx::cout << "s";
      }
      idfx::cout << " ";
    }
    if (n_hours > 0) {
      idfx::cout << n_hours << " hour";
      if (n_hours != 1) {
        idfx::cout << "s";
      }
      idfx::cout << " ";
    }
    if (n_minutes > 0) {
      idfx::cout << n_minutes << " minute";
      if (n_minutes != 1) {
        idfx::cout << "s";
      }
      idfx::cout << " ";
    }
    idfx::cout << n_seconds << " second";
    if (n_seconds != 1) {
      idfx::cout << "s";
    }
    idfx::cout << " ";
    idfx::cout << "and " << Tint->GetNCycles() << " cycle";
    if (Tint->GetNCycles() != 1) {
      idfx::cout << "s";
    }
    idfx::cout << std::endl;
    idfx::cout << "Main: ";
    idfx::cout << "Perfs are " << std::scientific << 1/perfs << " cell updates/second" << std::endl;
    #ifdef WITH_MPI
      idfx::cout << "MPI overhead represents "
                 << static_cast<int>(100.0*idfx::mpiCallsTimer/timer.seconds())
                 << "% of total run time." << std::endl;
    #endif

    idfx::cout << "Outputs represent "
               << static_cast<int>(100.0*output->GetTimer()/timer.seconds())
              << "% of total run time." << std::endl;
    // Show profiler output
    idfx::prof.Show();
}

void Menhir::PerformNewton() {
  idfx::cout << "I've been asked to perform Newton but it's not yet coded." << std::endl;
}

void Menhir::PerformStability() {
  idfx::cout << "I've been asked to perform a stability analysis but it's not yet coded." << std::endl;
}

void Menhir::PerformContinuation() {
  idfx::cout << "I've been asked to perform a continuation but it's not yet coded." << std::endl;
}

void Menhir::Finalise() {
}

