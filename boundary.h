#include <petscsys.h>
#pragma once

/* List of all boundary condition types */

enum boundary_type{
    DIRICHLET,
    NEUMANN
};


/* An abstract base class for boundary conditions */
class BoundaryCondition{
public:
    virtual ~BoundaryCondition() = 0;
};

/* Boundary condition "container" class for rectangular domain */
class Boundaries{
public:
    BoundaryCondition* left;
    BoundaryCondition* right;
    BoundaryCondition* bottom;
    BoundaryCondition* top;
};


/* Concrete boundary condition classes */
class Dirichlet : public BoundaryCondition{
public:
    Dirichlet(PetscReal x, PetscReal y);
    PetscReal x;
    PetscReal y;
private:
    boundary_type type;
};


class Neumann : public BoundaryCondition{
public:
    Neumann(PetscReal x, PetscReal y);
    PetscReal x;
    PetscReal y;
private:
    boundary_type type;
};
