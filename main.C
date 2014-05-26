#include "params.h"
#include "boundary.h"
#include "petsc2dNScpu.h"

#include <iostream>

int main(){

    // define the problem:

    Params params;
    Boundaries boundaries;

    params.N_x = 32;
    params.N_y = 32;
    params.dx  = 0.01;
    params.dy  = 0.01;
    params.Re  = 2000;
    params.dt  = 0.00001;
    params.nsteps = 20;

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
}
