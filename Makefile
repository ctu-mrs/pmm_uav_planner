include Mk/libs_flags.mk


OBJS=pmm_mg_trajectory3d.o pmm_trajectory.o pmm_trajectory3d.o common.o
OBJS_TESTS=
TARGET=main
TEST_TARGET=
OBJ_DIR=obj

include Mk/recipe.mk
