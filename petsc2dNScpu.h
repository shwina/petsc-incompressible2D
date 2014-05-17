#pragma once

#include <petscvec.h>
#include <petscmat.h>
#include <petscdm.h>
#include <petscdmda.h>

#include "boundary.h"
#include "params.h"
/*
    PETSC-based 2D Navier Stokes solver (CPU version)

*/

class Petsc2dNScpu{

public:
    Petsc2dNScpu(Boundaries&, Params&);
    void initialise();
    void takeStep();
    void finalise();
    void writeResults();
    ~Petsc2dNScpu();

private:
    void initialiseArrays();
    void makePoissonMatrix();
    void calculateRHS1();
    void calculateIntermediateVelocities();
    void calculateRHS2();
    void solvePoissonSystem();

Boundaries& boundaries;
Params& params;

/* PETSC field vectors: */
Vec u_local;
Vec u_global;
Vec v_local;
Vec v_global;
Vec p_local;
Vec p_global;

// Poisson matrix and associated RHS:
Mat poissonmatrix;
Vec rhs2;

// DM is a PETSc provided object that handles
// communication of data between processes
// for fields defined on rectangular, regular grids

// Let all field vectors share the same DM:
DM da;

PetscErrorCode ierr;
};