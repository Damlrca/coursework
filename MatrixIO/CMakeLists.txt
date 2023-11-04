add_library(MatrixIO STATIC MatrixIO.h MatrixIO.c mmio.h mmio.c)

target_include_directories(MatrixIO PUBLIC ../MatrixUtils)
target_link_libraries(MatrixIO PUBLIC MatrixUtils)
