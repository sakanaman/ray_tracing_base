add_library(texture texture.cpp)
target_include_directories(texture
        PRIVATE ${CMAKE_SOURCE_DIR}/ext/stb)