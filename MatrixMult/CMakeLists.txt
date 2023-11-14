add_library(MatrixMult STATIC MatrixMult.h MatrixMult.c)

target_include_directories(MatrixMult PUBLIC ../MatrixIO)
target_link_libraries(MatrixMult PUBLIC MatrixIO)

target_include_directories(MatrixMult PUBLIC ../MatrixUtils)
target_link_libraries(MatrixMult PUBLIC MatrixUtils)
