#include "boundary.h"

BoundaryCondition::~BoundaryCondition(){}

Dirichlet::Dirichlet(double& _x):x(_x){
    type = DIRICHLET;
}

Neumann::Neumann(double& _x):x(_x){
    type = NEUMANN;
}