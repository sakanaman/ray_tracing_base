add_library(parallel concurrent.cpp)
target_include_directories(parallel 
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/RNG)
target_link_libraries(parallel random)