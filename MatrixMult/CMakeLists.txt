add_library(MatrixMult STATIC MatrixMult.h MatrixMult.cpp MatrixMultTwoPhases.h MatrixMultTwoPhases.cpp)

target_include_directories(MatrixMult PUBLIC ../MatrixUtils)
target_link_libraries(MatrixMult PUBLIC MatrixUtils)
