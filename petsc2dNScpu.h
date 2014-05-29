#pragma once

#include <petscvec.h>
#include <petscmat.h>
#include <petscdm.h>
#include <petscksp.h>
#include <petscdmda.h>

#include "boundary.h"
#include "params.h"

#include <gtest/gtest_prod.h>

/*
    PETSC-based 2D Navier Stokes solver (CPU version)

*/

class Petsc2dNScpu{
public:
    Petsc2dNScpu(const Boundaries&, const Params&);
    void initialise();
    void takeStep();
    void finalise();
    void writeResults();
    ~Petsc2dNScpu();


private:
    FRIEND_TEST(Petsc2dNScpuTest, IntermediateVelocityTest);
    FRIEND_TEST(Petsc2dNScpuTest, PoissonConstructionTest);
    void initialiseArrays();
    void makePoissonMatrix();
    void calculateIntermediateVelocities();
    void calculateRHS2();
    void solvePoissonSystem();
    void updateVelocities();
    void updateGhosts();

    const Boundaries& boundaries;
    const Params& params;
    
    // DM is a PETSc provided object that handles
    // communication of data between processes
    // for fields defined on rectangular, regular grids

    // Let all field vectors share the same DM:
    DM da;

    // Linear solver context:
    KSP ksp;

    // Null space context:
    MatNullSpace nullspace;

    // PETSC field vectors:
    Vec u_local;
    Vec u_global;
    Vec v_local;
    Vec v_global;
    Vec p_local;
    Vec p_global;

    // RHS for computing intermediate velocities:
    Vec rhs1x;
    Vec rhs1y;

    // Poisson matrix and associated RHS:
    Mat poissonmatrix;
    Vec rhs2;
};