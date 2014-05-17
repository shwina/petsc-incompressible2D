#include "params.h"
#include "boundary.h"
#include "petsc2dNScpu.h"

int main(){

    // define the problem:

    Params params;
    Boundaries boundaries;

    params.N_x = 5;
    params.N_y = 5;
    params.dx  = 1;
    params.dy  = 1;
    params.Re  = 2;
    params.dt  = 0.1;
    params.nsteps = 1;

    double x = 5;
    boundaries.left  = new Dirichlet(x);
    boundaries.right = new Neumann(x);

    Petsc2dNScpu solver(boundaries, params);
    solver.initialise();
}