#include "boundary.h"

BoundaryCondition::~BoundaryCondition(){}

Dirichlet::Dirichlet(double _x, double _y):x(_x), y(_y){
    type = DIRICHLET;
}

Neumann::Neumann(double _x, double _y):x(_x), y(_y){
    type = NEUMANN;
}