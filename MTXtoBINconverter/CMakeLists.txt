add_executable(MTXtoBINconverter main.cpp)

target_include_directories(MTXtoBINconverter PUBLIC ../MatrixIO)
target_link_libraries(MTXtoBINconverter PUBLIC MatrixIO)

target_include_directories(MTXtoBINconverter PUBLIC ../MatrixUtils)
target_link_libraries(MTXtoBINconverter PUBLIC MatrixUtils)

target_include_directories(MTXtoBINconverter PUBLIC ../OtherUtils)
target_link_libraries(MTXtoBINconverter PUBLIC OtherUtils)
