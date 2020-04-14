
GCC_FLAGS="-O3 -std=c++11"
CC="g++"
SRC_DIR = scIB/knn_graph
SRC = $(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/*.cpp))
OBJS = $(patsubst src/%.cpp,$(SRC_DIR)/%.o,$(SRC))

all:
   $CC $GCC_FLAGS ${SRC} -o ${OBJS}
