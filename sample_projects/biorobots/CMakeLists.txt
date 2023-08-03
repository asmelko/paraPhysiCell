cmake_minimum_required(VERSION 3.18)

project(biorobots VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE src_files "*")
add_executable(biorobots ${src_files})
target_include_directories(biorobots PRIVATE "../..")
target_link_libraries(biorobots PRIVATE PhysiCellCore)

if(MSVC)
    target_compile_options(biorobots PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/bigobj>)
else()
    target_compile_options(biorobots PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
endif()

add_custom_command(TARGET biorobots POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/config $<TARGET_FILE_DIR:${PROJECT_NAME}>/config)

# add_custom_command(TARGET biorobots POST_BUILD
#     COMMAND ${CMAKE_COMMAND} -E make_directory
#     $<TARGET_FILE_DIR:${PROJECT_NAME}>/output)
