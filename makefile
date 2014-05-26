#
CFLAGS		 =
FFLAGS		 =
CPPFLAGS     = -I .
FPPFLAGS     =

GTEST_DIR = /home/ashwin/workspace/gtest-1.7.0
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

CPPFLAGS += -I $(GTEST_DIR)/include -I ..
CXXFLAGS += -g -Wall -Wextra -pthread



include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

ROOT_DIR=/home/ashwin/workspace/projects/Incompressible3D/petsc-cpu

OBJ = $(ROOT_DIR)/petsc2dNScpu.o $(ROOT_DIR)/boundary.o

petsc2dNScpu: $(OBJ) $(ROOT_DIR)/main.o chkopts
	-${CLINKER} -o petsc2dNScpu $(OBJ) $(ROOT_DIR)/main.o $(PETSC_LIB) -I .
