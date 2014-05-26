#include "petsc2dNScpu.h"
#include <gtest/gtest.h>


class Petsc2dNScpuTest : public ::testing::Test{

protected:
  
    Petsc2dNScpu *solver;
    Params params;
    Boundaries boundaries;

    virtual void SetUp(){
        params.dx  = (PetscReal)1;
        params.dy  = (PetscReal)1;
        params.dt  = 1.;
        params.Re  = (PetscReal)1;
        params.N_x = 5;
        params.N_y = 5;

        // arbitrary boundary conditions are
        // fine
        
        PetscReal zero = 0;
        boundaries.left   = new Dirichlet(zero, zero);
        boundaries.right  = new Dirichlet(zero, zero);
        boundaries.top    = new Dirichlet(zero, zero);
        boundaries.bottom = new Dirichlet(zero, zero);
        solver = new Petsc2dNScpu(boundaries, params);
        solver->initialise();
    }

    void TearDown(){
    }   

};


TEST_F(Petsc2dNScpuTest, PoissonConstruction){
    
    int rows[1];
    int cols[1];
    PetscReal vals[1];

    // check that the diagonal contribution
    // by the centre cell is 4

    rows[0] = 12;
    cols[0] = 12;

    MatGetValues(solver->poissonmatrix, 1, rows, 1, cols, vals);
    ASSERT_EQ(vals[0], 4);

    // edge
    rows[0] = 7;
    cols[0] = 7;

    MatGetValues(solver->poissonmatrix, 1, rows, 1, cols, vals);
    ASSERT_EQ(vals[0], 3);

    // corner
    rows[0] = 6;
    cols[0] = 6;

    MatGetValues(solver->poissonmatrix, 1, rows, 1, cols, vals);
    ASSERT_EQ(vals[0], 2);

    // padding
    rows[0] = 0;
    cols[0] = 0;

    MatGetValues(solver->poissonmatrix, 1, rows, 1, cols, vals);
    ASSERT_EQ(vals[0], 1);

}


TEST_F(Petsc2dNScpuTest, IntermediateVelocityTest){

    
    PetscReal **u_l, **v_l;
    PetscInt x, y, m, n;

    DMDAVecGetArray(solver->da, solver->u_local, &u_l);
    DMDAVecGetArray(solver->da, solver->v_local, &v_l);

    DMDAGetCorners(solver->da, &x, &y, NULL, &m, &n, NULL);

    // set up some initial conditions for the
    // velocities:

    for (int i = y; i < y+n; i++){
        for (int j = x; j < x+m; j++){
            u_l[i][j] = (PetscReal)j;
            v_l[i][j] = (PetscReal)i;
        }
    }

    DMDAVecRestoreArray(solver->da, solver->u_local, &u_l);
    DMDAVecRestoreArray(solver->da, solver->v_local, &v_l);

    solver->calculateIntermediateVelocities();

    DMDAVecGetArray(solver->da, solver->u_local, &u_l);
    DMDAVecGetArray(solver->da, solver->u_local, &v_l);
    
    // u-velocity check
    ASSERT_EQ(u_l[2][2], -4); // centre (computed)
    ASSERT_EQ(u_l[2][3],  3); // boundary
    ASSERT_EQ(u_l[2][4],  4); // padding

    // v-velocity check
    ASSERT_EQ(v_l[2][2], -4); // centre (computed)

    DMDAVecRestoreArray(solver->da, solver->u_local, &u_l);
    DMDAVecRestoreArray(solver->da, solver->v_local, &v_l);

}


GTEST_API_ int main(int argc, char **argv) {
  printf("Running main() from gtest_main.cc\n");
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}