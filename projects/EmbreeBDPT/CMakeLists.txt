add_executable(bdpt_embree main.cpp)
target_include_directories(bdpt_embree
PRIVATE ${CMAKE_SOURCE_DIR}/ext/tinyobjloader
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/Camera
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/RNG
    PRIVATE ${CMAKE_SOURCE_DIR}/src/Parallel
    PRIVATE ${CMAKE_SOURCE_DIR}/src/Embree
    PRIVATE ${CMAKE_SOURCE_DIR}/src/Light
    PRIVATE ${CMAKE_SOURCE_DIR}/ext/stb
    PRIVATE ${CMAKE_SOURCE_DIR}/lib/Math
    PRIVATE ${CMAKE_SOURCE_DIR}/src/Texture)
target_link_libraries(bdpt_embree 
                            stb_image_write 
                            Camera 
                            random
                            IBL 
                            parallel 
                            Vec
                            embree_intersector
                            tinyobjloader
                            stb_image
                            Threads::Threads
                            texture)