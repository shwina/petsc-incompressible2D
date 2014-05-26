#include <petscsys.h>

#pragma once

class Params{

public:
    // grid:
    int N_x, N_y;
    PetscReal dx, dy;
    
    // simulation parameters:
    PetscReal Re, dt;
    int nsteps;
};