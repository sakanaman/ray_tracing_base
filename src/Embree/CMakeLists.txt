FIND_PACKAGE(embree 3.0 REQUIRED)
add_library(embree_intersector wrap_embree.cpp)

    target_include_directories(embree_intersector
                            PRIVATE ${CMAKE_SOURCE_DIR}/lib/Math
                            PUBLIC ${CMAKE_SOURCE_DIR}/ext/embree-3.13.2.x86_64.linux/include
                            PRIVATE ${CMAKE_SOURCE_DIR}/src/Embree)

target_link_libraries(embree_intersector Vec embree)