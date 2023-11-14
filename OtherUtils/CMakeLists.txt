add_library(OtherUtils STATIC OtherUtils.hpp OtherUtils.cpp)

target_include_directories(OtherUtils PUBLIC ../MatrixUtils)
target_link_libraries(OtherUtils PUBLIC MatrixUtils)
