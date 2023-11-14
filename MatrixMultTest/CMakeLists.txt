add_executable(MatrixMultTest main.cpp)

target_include_directories(MatrixMultTest PUBLIC ../MatrixIO)
target_link_libraries(MatrixMultTest PUBLIC MatrixIO)

target_include_directories(MatrixMultTest PUBLIC ../MatrixUtils)
target_link_libraries(MatrixMultTest PUBLIC MatrixUtils)

target_include_directories(MatrixMultTest PUBLIC ../MatrixMult)
target_link_libraries(MatrixMultTest PUBLIC MatrixMult)

target_include_directories(MatrixMultTest PUBLIC ../OtherUtils)
target_link_libraries(MatrixMultTest PUBLIC OtherUtils)
