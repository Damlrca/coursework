add_executable(MatrixGenerator main.cpp)

target_include_directories(MatrixGenerator PUBLIC ../MatrixIO)
target_link_libraries(MatrixGenerator PUBLIC MatrixIO)

target_include_directories(MatrixGenerator PUBLIC ../MatrixUtils)
target_link_libraries(MatrixGenerator PUBLIC MatrixUtils)
