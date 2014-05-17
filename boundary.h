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

/* Boundary condition "container" class for square domain */
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
    Dirichlet(double& x);

private:
    double& x;
    boundary_type type;
};

class Neumann : public BoundaryCondition{
public:
    Neumann(double& x);

private:
    double& x;
    boundary_type type;
};