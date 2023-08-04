cmake_minimum_required(VERSION 3.18)

project(heterogeneity VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE src_files "*")
add_executable(heterogeneity ${src_files})
target_include_directories(heterogeneity PRIVATE "../..")
target_link_libraries(heterogeneity PRIVATE PhysiCellCore)

if(MSVC)
    target_compile_options(heterogeneity PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/W4 /bigobj>)
else()
    target_compile_options(heterogeneity PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -pedantic>)
endif()

add_custom_command(TARGET heterogeneity POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/config $<TARGET_FILE_DIR:${PROJECT_NAME}>/config)
