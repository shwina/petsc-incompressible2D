#include "params.h"
#include "boundary.h"
#include "petsc2dNScpu.h"

#include <iostream>

int main(int argc, char** args){

    PetscInitialize(&argc, &args, NULL, NULL);

    // define the problem:

    Params params;
    Boundaries boundaries;

    params.N_x = 128;
    params.N_y = 128;
    params.dx  = 0.0025;
    params.dy  = 0.0025;
    params.Re  = 2000;
    params.dt  = 0.0001;
    params.nsteps = 100;

/*
    boundaries.left  = new Dirichlet(0, 0);
    boundaries.right = new Dirichlet(0, 0);
    boundaries.top   = new Dirichlet(0.001, 0);
    boundaries.bottom = new Dirichlet(0, 0);
*/

    Petsc2dNScpu solver(boundaries, params);
    solver.initialise();

    for (int i = 0; i < params.nsteps; i++){
        solver.takeStep();
    }

    solver.writeResults();
    solver.finalise();
    PetscFinalize();
}
