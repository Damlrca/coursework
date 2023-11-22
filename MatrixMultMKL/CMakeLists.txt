add_executable(MatrixMultMKL main.cpp)

target_include_directories(MatrixMultMKL PUBLIC ../MatrixIO)
target_link_libraries(MatrixMultMKL PUBLIC MatrixIO)

target_include_directories(MatrixMultMKL PUBLIC ../MatrixUtils)
target_link_libraries(MatrixMultMKL PUBLIC MatrixUtils)

target_include_directories(MatrixMultMKL PUBLIC ../OtherUtils)
target_link_libraries(MatrixMultMKL PUBLIC OtherUtils)
