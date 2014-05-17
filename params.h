

#pragma once

class Params{

public:

    // grid:
    int N_x, N_y;
    double dx, dy;
    
    // simulation parameters:
    double Re, dt;
    int nsteps;
};