cmake_minimum_required(VERSION 3.18)

project(project VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE src_files "*")
add_executable(project ${src_files})
target_include_directories(project PRIVATE "../..")
target_link_libraries(project PRIVATE PhysiCellCore)

if(MSVC)
    target_compile_options(project PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/W4 /bigobj>)
else()
    target_compile_options(project PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -pedantic>)
endif()

add_custom_command(TARGET project POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/config $<TARGET_FILE_DIR:${PROJECT_NAME}>/config)
