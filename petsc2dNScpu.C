#include <petsc2dNScpu.h>
#include <iostream>

/*
    2-D Incompressible Flow solver based
    on the fractional step method.

    A "padding layer" is introduced to keep
    the u, v, and p arrays identically sized.
*/



Petsc2dNScpu::Petsc2dNScpu(Boundaries& _boundaries, Params& _params) :
    boundaries(_boundaries), params(_params){} 


void Petsc2dNScpu::initialise(){
    PetscInitialize(NULL, NULL, NULL, NULL); 
    initialiseArrays();
    makePoissonMatrix();

    // remove these:
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "mat.txt", &viewer);
    MatView(poissonmatrix, viewer);


}

void Petsc2dNScpu::takeStep(){}
void Petsc2dNScpu::writeResults(){}

void Petsc2dNScpu::initialiseArrays(){

    // create the DA (distributed array) object, with DOF 1 and stencil width 1

    DMDACreate2d(PETSC_COMM_WORLD, 
                        DMDA_BOUNDARY_GHOSTED, DMDA_BOUNDARY_GHOSTED,
                        DMDA_STENCIL_STAR, 
                        params.N_x, params.N_y, 
                        PETSC_DECIDE, PETSC_DECIDE, 
                        1,
                        1, 
                        PETSC_NULL, PETSC_NULL,
                        &da); 

    // create field vectors (global and local copies)

    DMCreateLocalVector(da, &u_local); 
    DMCreateLocalVector(da, &v_local); 
    DMCreateLocalVector(da, &p_local); 

    DMCreateGlobalVector(da, &u_global); 
    DMCreateGlobalVector(da, &v_global); 
    DMCreateGlobalVector(da, &p_global); 

    // poisson matrix and RHS:

    DMCreateMatrix(da, MATMPIAIJ, &poissonmatrix); 
    DMCreateGlobalVector(da, &rhs2); 

    // Set the initial values:

    VecSet(u_local, 0); 
    VecSet(v_local, 0); 
    VecSet(rhs2, 0); 
}

void Petsc2dNScpu::makePoissonMatrix(){
    
    PetscInt x, y, m, n;
    MatStencil row, col[5];
    PetscReal v[5];
    PetscInt numi, numj, num;

    PetscReal dx, dy, N_x, N_y;
    dx = params.dx;
    dy = params.dy;
    N_x = params.N_x;
    N_y = params.N_y; 

    // get the corner indices of the process (global) 
    // and the widths in each direction:
    
    DMDAGetCorners(da, &x, &y, NULL, &m, &n, NULL);

    // loop over the grid points belonging to this process:
    for (int j = y; j < y+n; j++){
        for (int i = x; i < x+m; i++){

            // the row of the matrix corresponds
            // to the (i, j) grid point
            
            row.i = i; row.j = j;
            numi = 0; numj = 0; num = 0;

            if (i == 0 || j == 0 || i == N_x-1 || j == N_y-1){
                // this is padding layer:
               
                col[0].i = i; col[0].j = j; v[0] = 0; // arbitrarily set to 0.
                num++;
                MatSetValuesStencil(poissonmatrix, 1, &row, num, col, v, INSERT_VALUES); 
            }

            else if(i == 1 || j == 1 || i == N_x-2 || j == N_y-2){ 
                // boundaries:

                if (i != 1){
                    col[num].i = i-1; // this point is 
                    col[num].j = j;   // inside domain
                    v[num] = -1/(dx*dx);
                    numi++; num++;
                }

                if (i != N_x-2){
                    col[num].i = i+1;
                    col[num].j = j;
                    v[num] = -1/(dx*dx);
                    numi++; num++;
                }

                if (j != 1){
                    col[num].j = j-1;
                    col[num].i = i;
                    v[num] = -1/(dy*dy);
                    numj++; num++;
                }

                if (j != N_y-2){
                    col[num].j = j+1;
                    col[num].i = i;
                    v[num] = -1/(dy*dy);
                    numj++; num++;
                }

                col[num].i = i; col[num].j = j;
                v[num] = numi*(1/(dx*dx)) + numj*(1/(dy*dy));
                num++; // we want absolute value

                MatSetValuesStencil(poissonmatrix, 1, &row, num, col, v, INSERT_VALUES);
                

            }

            else{
                col[0].i = i-1; col[0].j =   j; v[0] = -1/(dx*dx);
                col[1].i = i+1; col[1].j =   j; v[1] = -1/(dx*dx);
                col[2].i =   i; col[2].j =   j; v[2] =  2/(dx*dx) + 2/(dy*dy);
                col[3].i =   i; col[3].j = j-1; v[3] = -1/(dy*dy);
                col[4].i =   i; col[4].j = j+1; v[4] = -1/(dy*dy);
                MatSetValuesStencil(poissonmatrix, 1, &row, 5, col, v, INSERT_VALUES); 
                
            }
        }
    }

    MatAssemblyBegin(poissonmatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(poissonmatrix, MAT_FINAL_ASSEMBLY);
}

void Petsc2dNScpu::calculateRHS1(){}
void Petsc2dNScpu::calculateIntermediateVelocities(){}
void calculateRHS2(){}

Petsc2dNScpu::~Petsc2dNScpu(){

    VecDestroy(&u_local); 
    VecDestroy(&v_local); 
    VecDestroy(&p_local); 

    VecDestroy(&u_global); 
    VecDestroy(&v_global); 
    VecDestroy(&p_global); 

    MatDestroy(&poissonmatrix);

    PetscFinalize();
}





