cmake_minimum_required(VERSION 3.18)

project(
  biorobots
  VERSION 0.1
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE src_files "*")
add_executable(biorobots ${src_files})
target_include_directories(biorobots PRIVATE "../..")
target_link_libraries(biorobots PRIVATE PhysiCellCore)

add_custom_command(
  TARGET biorobots
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/config
          $<TARGET_FILE_DIR:${PROJECT_NAME}>/config)
