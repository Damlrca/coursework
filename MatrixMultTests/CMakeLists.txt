add_executable(MatrixMultTests main.cpp)

target_include_directories(MatrixMultTests PUBLIC ../MatrixIO)
target_link_libraries(MatrixMultTests PUBLIC MatrixIO)

target_include_directories(MatrixMultTests PUBLIC ../MatrixUtils)
target_link_libraries(MatrixMultTests PUBLIC MatrixUtils)

target_include_directories(MatrixMultTests PUBLIC ../MatrixMult)
target_link_libraries(MatrixMultTests PUBLIC MatrixMult)

target_include_directories(MatrixMultTests PUBLIC ../MatrixMultMKL)
target_link_libraries(MatrixMultTests PUBLIC MatrixMultMKL)

target_include_directories(MatrixMultTests PUBLIC ../OtherUtils)
target_link_libraries(MatrixMultTests PUBLIC OtherUtils)
