cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

project(concave_hull)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(BUILD_EXAMPLE FALSE)

add_library(${PROJECT_NAME}
    SHARED
    ${CMAKE_CURRENT_SOURCE_DIR}/concave_hull.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/delaunator.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/delaunator.cpp
)

target_include_directories(${PROJECT_NAME}
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(concave_hull::concave_hull
    ALIAS
    ${PROJECT_NAME}
)

if(BUILD_EXAMPLE)
    find_package(Python REQUIRED
        COMPONENTS
        Interpreter
        Development
    )
    add_executable(main main.cpp)
    target_include_directories(main
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_link_libraries(main
        PRIVATE
        concave_hull
        Python::Python
    )
endif()