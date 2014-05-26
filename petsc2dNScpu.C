#include <petsc2dNScpu.h>
#include <iostream>
#include <petscdm.h>
#include <petscsys.h>



/*
    2-D Incompressible Flow solver based
    on the fractional step method on a staggered
    grid (with padding)

*/

Petsc2dNScpu::Petsc2dNScpu(const Boundaries& _boundaries, const Params& _params) :
    boundaries(_boundaries), params(_params){} 


void Petsc2dNScpu::initialise(){
    PetscInitialize(NULL, NULL, NULL, NULL); 
    initialiseArrays();
    makePoissonMatrix();
}

void Petsc2dNScpu::takeStep(){
    updateGhosts();
    calculateIntermediateVelocities();
    updateGhosts();
    calculateRHS2();
    solvePoissonSystem();
    updateVelocities();    
}


void Petsc2dNScpu::writeResults(){

    PetscViewer viewer;

    DMLocalToGlobalBegin(da, u_local, INSERT_VALUES, u_global);
    DMLocalToGlobalEnd(da, u_local, INSERT_VALUES,u_global);

    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "u.txt", &viewer);
    VecView(u_global, viewer);


    DMLocalToGlobalBegin(da, v_local, INSERT_VALUES, v_global);
    DMLocalToGlobalEnd(da, v_local, INSERT_VALUES, v_global);

    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "v.txt", &viewer);
    VecView(v_global, viewer);
}


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

    // right hand side for intermediate velocities:

    DMCreateGlobalVector(da, &rhs1x);
    DMCreateGlobalVector(da, &rhs1y);

    // poisson matrix and RHS:

    DMCreateMatrix(da, MATMPIAIJ, &poissonmatrix); 
    DMCreateGlobalVector(da, &rhs2); 

    // Linear solver:
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPCG);

    // Set the initial values:

    /* TODO incorporate initial conditions*/
    VecSet(u_local, 0); 
    VecSet(v_local, 0); 
    VecSet(rhs2, 0);

    // Null space:
    MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace);
}

void Petsc2dNScpu::makePoissonMatrix(){

    
    // Based on the following PETSc example code:
    // http://www.mcs.anl.gov/petsc/petsc-current/src/ksp/ksp/examples/tutorials/ex34.c.html

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
                col[0].i = i; col[0].j = j; v[0] = 1/(dx*dy); // we want decent scaling
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

    KSPSetOperators(ksp, poissonmatrix, poissonmatrix, SAME_NONZERO_PATTERN);
}


void Petsc2dNScpu::calculateIntermediateVelocities(){

    PetscReal dx, dy, N_x, N_y;
    dx = params.dx;
    dy = params.dy;
    N_x = params.N_x;
    N_y = params.N_y; 

    // compute the convection and diffusion terms:

    PetscReal **u_l, **v_l, **rx, **ry;
    PetscInt x, y, m, n;
    
    DMDAVecGetArray(da, u_local, &u_l);
    DMDAVecGetArray(da, v_local, &v_l);
    DMDAVecGetArray(da, rhs1x,   &rx);
    DMDAVecGetArray(da, rhs1y,   &ry);

    DMDAGetCorners(da, &x, &y, NULL, &m, &n, NULL);

    // loop over the grid points belonging to this process:
    for (int i = y; i < y+n; i++){
        for (int j = x; j < x+m; j++){

            // u is padded on the right, and we don't update
            // boundaries
            if (i == 0 || j == 0){
                // do nothing
            }

            else{

                // account for u-padding:
                if (j < N_x - 2){

                    // u-velocity update:
                    rx[i][j] =  ((u_l[i][j] + u_l[i][j-1])*(u_l[i][j] + u_l[i][j-1]) -
                                 (u_l[i][j+1] + u_l[i][j])*(u_l[i][j+1] + u_l[i][j]))/(4*dx) +

                                ((u_l[i][j] + u_l[i-1][j])*(v_l[i-1][j] + v_l[i-1][j+1]) - 
                                 (u_l[i+1][j] + u_l[i][j])*(v_l[i][j] + v_l[i][j+1]))/(4*dy) +

                                (1./(params.Re))*
                                ((u_l[i][j+1] - 2*u_l[i][j] + u_l[i][j-1])/(dx*dx) +
                                 (u_l[i+1][j] - 2*u_l[i][j] + u_l[i-1][j])/(dy*dy));
                }

                // account for v-padding:
                if (i < N_y - 2){

                    // v-velocity update:
                    ry[i][j]  = ((v_l[i][j] + v_l[i-1][j])*(v_l[i][j] + v_l[i-1][j]) - 
                                 (v_l[i+1][j] + v_l[i][j])*(v_l[i+1][j] + v_l[i][j]))/(4*dy) + 

                                ((v_l[i][j-1] + v_l[i][j])*(u_l[i][j-1] + u_l[i+1][j-1]) -
                                 (v_l[i][j] + v_l[i][j+1])*(u_l[i][j] + u_l[i+1][j]))/(4*dx) + 

                                (1./(params.Re))*
                                ((v_l[i][j+1] - 2*v_l[i][j] + v_l[i][j-1])/(dx*dx) + 
                                 (v_l[i+1][j] - 2*v_l[i][j] + v_l[i-1][j])/(dy*dy));
                } 
            } 
        } 
    }

    // compute intermediate velocities from the RHS:
    for (int i = y; i < y+n; i++){
        for (int j = x; j < x+m; j++){

            if (i == 0 || j == 0){
                // do nothing
            }


            else{

                // account for u-padding:
                if (j < N_x - 2){
                    u_l[i][j] += params.dt*rx[i][j];
                }

                // account for v-padding:
                if (i < N_y - 2){
                    v_l[i][j] += params.dt*ry[i][j];
                } 
            } 
        } 
    }        

    DMDAVecRestoreArray(da, u_local, &u_l);
    DMDAVecRestoreArray(da, v_local, &v_l);
    DMDAVecRestoreArray(da, rhs1x,   &rx);
    DMDAVecRestoreArray(da, rhs1y,   &ry);
}

void Petsc2dNScpu::calculateRHS2(){

    PetscReal dx, dy, N_x, N_y;
    dx = params.dx;
    dy = params.dy;
    N_x = params.N_x;
    N_y = params.N_y; 

    PetscReal **u_l, **v_l, **r2;
    PetscInt x, y, m, n;

    DMDAVecGetArray(da, u_local,  &u_l);
    DMDAVecGetArray(da, v_local,  &v_l);
    DMDAVecGetArray(da, rhs2,     &r2);
    
    DMDAGetCorners(da, &x, &y, NULL, &m, &n, NULL);

    // loop over the grid points belonging to this process:
    for (int i = y; i < y+n; i++){
        for (int j = x; j < x+m; j++){
        
            if (i == 0 || j == 0 || i == N_x-1 || j == N_y-1){
                // padding
            }

            else{  
                r2[i][j] = -((u_l[i][j] - u_l[i][j-1])/dx + (v_l[i][j] - v_l[i-1][j])/dy)/params.dt;
            }
        }
    }

    DMDAVecRestoreArray(da, u_local, &u_l);
    DMDAVecRestoreArray(da, v_local, &v_l);
    DMDAVecRestoreArray(da, rhs2,    &r2);
}

void Petsc2dNScpu::solvePoissonSystem(){
    MatNullSpaceRemove(nullspace, rhs2, NULL);
    KSPSolve(ksp, rhs2, p_global);
}

void Petsc2dNScpu::updateVelocities(){

    PetscReal dx, dy, N_x, N_y;
    dx = params.dx;
    dy = params.dy;
    N_x = params.N_x;
    N_y = params.N_y; 

    PetscReal **p_l, **u_l, **v_l;
    PetscInt x, y, m, n;

    DMGlobalToLocalBegin(da, p_global, INSERT_VALUES, p_local);
    DMGlobalToLocalEnd(da, p_global, INSERT_VALUES, p_local);

    DMDAVecGetArray(da, u_local, &u_l);
    DMDAVecGetArray(da, v_local, &v_l);
    DMDAVecGetArray(da, p_local, &p_l);

    DMDAGetCorners(da, &x, &y, NULL, &m, &n, NULL);

    for (int i = y; i < y+n; i++){
        for (int j = x; j < x+m; j++){

            // u is padded on the right, and we don't update
            // boundaries
            if (i == 0 || j == 0){
                // do nothing
            }

            else{

                if (j < N_x - 2){
                    u_l[i][j] -= params.dt*(p_l[i][j] - p_l[i][j-1])/dx;
                }

                if (i < N_y - 2){
                    v_l[i][j] -= params.dt*(p_l[i][j] - p_l[i-1][j])/dy;
                } 
            } 
        } 
    }

    DMDAVecRestoreArray(da, u_local, &u_l);
    DMDAVecRestoreArray(da, v_local, &v_l);
    DMDAVecRestoreArray(da, p_local, &p_l);

    DMLocalToGlobalBegin(da, p_local, INSERT_VALUES, p_global);
    DMLocalToGlobalEnd(da, p_local, INSERT_VALUES, p_global);
}


void Petsc2dNScpu::updateGhosts(){
    // local-to-local update of
    // ghosts:
    DMDALocalToLocalBegin(da, u_local, INSERT_VALUES, u_local);
    DMDALocalToLocalEnd(da, u_local, INSERT_VALUES, u_local);

    DMDALocalToLocalBegin(da, v_local, INSERT_VALUES, v_local);
    DMDALocalToLocalEnd(da, v_local, INSERT_VALUES, v_local);

    PetscReal dx, dy, N_x, N_y;
    dx = params.dx;
    dy = params.dy;
    N_x = params.N_x;
    N_y = params.N_y; 

    PetscReal **u_l, **v_l;
    PetscInt x, y, m, n;   

    DMDAVecGetArray(da, u_local, &u_l);
    DMDAVecGetArray(da, v_local, &v_l);

    DMDAGetCorners(da, &x, &y, NULL, &m, &n, NULL);


    // left and right (only v needs update)
    int i, j;
    for (i = y; i < y+n; i++){
        j = x;
        if(j == 0){
            v_l[i][j] = 0 - v_l[i][j+1];
        }

        j = (x+m)-2;
        if(j == N_x-2){
            v_l[i][j] = 0 - v_l[i][j-1];
        }
    }

    // bottom and top (only u needs update)
    for (j = x; j < x+m; j++){
        
        i = y;
        if (i == 0){
            u_l[i][j] = 0 - u_l[i+1][j];
        }

        i = (y+n)-2;
        if (i == N_y-2){
            u_l[i][j] = 0.00002 - u_l[i-1][j];
        }
    }

    DMDAVecRestoreArray(da, u_local, &u_l);
    DMDAVecRestoreArray(da, v_local, &v_l);
}

Petsc2dNScpu::~Petsc2dNScpu(){

    VecDestroy(&u_local); 
    VecDestroy(&v_local); 
    VecDestroy(&p_local); 

    VecDestroy(&u_global); 
    VecDestroy(&v_global); 
    VecDestroy(&p_global); 

    VecDestroy(&rhs1x);
    VecDestroy(&rhs1y);

    VecDestroy(&rhs2);
    MatDestroy(&poissonmatrix);

    KSPDestroy(&ksp);
    MatNullSpaceDestroy(&nullspace);

    PetscFinalize();
}





