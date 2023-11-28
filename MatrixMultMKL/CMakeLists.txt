add_library(MatrixMultMKL STATIC MatrixMultMKL.h MatrixMultMKL.c)

target_include_directories(MatrixMultMKL PUBLIC ../MatrixUtils)
target_link_libraries(MatrixMultMKL PUBLIC MatrixUtils)
