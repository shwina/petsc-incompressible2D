

void convectionKernel_u(PetscReal **u_l, PetscReal **v_l, PetscReal **rx,  PetscInt i, PetscInt j, PetscReal dx, PetscReal dy, PetscReal Re) {

    // Apply the 5-4 pt. convection stencil at some
    // point i,j

    rx[i][j] =  ((u_l[i][j] + u_l[i][j-1])*(u_l[i][j] + u_l[i][j-1]) -
                 (u_l[i][j+1] + u_l[i][j])*(u_l[i][j+1] + u_l[i][j]))/(4*dx) +

                ((u_l[i][j] + u_l[i-1][j])*(v_l[i-1][j] + v_l[i-1][j+1]) - 
                 (u_l[i+1][j] + u_l[i][j])*(v_l[i][j] + v_l[i][j+1]))/(4*dy) +

                (1./Re)*
                ((u_l[i][j+1] - 2*u_l[i][j] + u_l[i][j-1])/(dx*dx) +
                 (u_l[i+1][j] - 2*u_l[i][j] + u_l[i-1][j])/(dy*dy));
}

void convectionKernel_v(PetscReal **u_l, PetscReal **v_l, PetscReal **ry,  PetscInt i, PetscInt j, PetscReal dx, PetscReal dy, PetscReal Re) {

    // Apply the 5-4 pt. convection stencil at some 
    // point i,j

    ry[i][j]  = ((v_l[i][j] + v_l[i-1][j])*(v_l[i][j] + v_l[i-1][j]) - 
                 (v_l[i+1][j] + v_l[i][j])*(v_l[i+1][j] + v_l[i][j]))/(4*dy) + 

                ((v_l[i][j-1] + v_l[i][j])*(u_l[i][j-1] + u_l[i+1][j-1]) -
                 (v_l[i][j] + v_l[i][j+1])*(u_l[i][j] + u_l[i+1][j]))/(4*dx) + 

                (1./Re)*
                ((v_l[i][j+1] - 2*v_l[i][j] + v_l[i][j-1])/(dx*dx) + 
                 (v_l[i+1][j] - 2*v_l[i][j] + v_l[i-1][j])/(dy*dy));
}

void divergenceKernel(PetscReal **u_l, PetscReal **v_l, PetscReal **r2, PetscInt i, PetscInt j, PetscReal dx, PetscReal dy) {
    r2[i][j] =  (u_l[i][j] - u_l[i][j-1])/dx + (v_l[i][j] - v_l[i-1][j])/dy;
}

void gradientKernel_u(PetscReal **p_l, PetscReal **u_l, PetscInt i, PetscInt j, PetscReal dx, PetscReal dt) {
    u_l[i][j] -= dt*(p_l[i][j+1] - p_l[i][j])/dx;
}

void gradientKernel_v(PetscReal **p_l, PetscReal **v_l, PetscInt i, PetscInt j, PetscReal dy, PetscReal dt) {
    v_l[i][j] -= dt*(p_l[i+1][j] - p_l[i][j])/dy;
}

