add_executable(path_tracing_test_embree main.cpp)
target_include_directories(path_tracing_test_embree
PRIVATE ${CMAKE_SOURCE_DIR}/ext/tinyobjloader
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/Camera
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/RNG
    PRIVATE ${CMAKE_SOURCE_DIR}/src/Parallel
    PRIVATE ${CMAKE_SOURCE_DIR}/src/Embree
    PRIVATE ${CMAKE_SOURCE_DIR}/src/Light
    PRIVATE ${CMAKE_SOURCE_DIR}/ext/stb
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/Math)
target_link_libraries(path_tracing_test_embree 
                                            stb_image_write 
                                            Camera 
                                            random
                                            IBL 
                                            parallel 
                                            Vec
                                            embree_intersector
                                            tinyobjloader
                                            stb_image
                                            Threads::Threads)