add_executable(MatrixMatrixMult main.cpp)

target_include_directories(MatrixMatrixMult PUBLIC ../MatrixIO)
target_link_libraries(MatrixMatrixMult PUBLIC MatrixIO)
