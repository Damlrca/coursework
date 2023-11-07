add_executable(MatrixTests main.cpp)

target_include_directories(MatrixTests PUBLIC ../MatrixIO)
target_link_libraries(MatrixTests PUBLIC MatrixIO)

target_include_directories(MatrixTests PUBLIC ../MatrixUtils)
target_link_libraries(MatrixTests PUBLIC MatrixUtils)
