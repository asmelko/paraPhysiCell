cmake_minimum_required(VERSION 3.18)

project(
  sem-intrusion-ring
  VERSION 0.1
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE src_files "*")
add_executable(sem-intrusion-ring ${src_files} ${CMAKE_CURRENT_SOURCE_DIR}/../helpers.cpp)
target_include_directories(sem-intrusion-ring PRIVATE "../..")
target_link_libraries(sem-intrusion-ring PRIVATE PhysiCellCore)

add_custom_command(
  TARGET sem-intrusion-ring
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/config
          $<TARGET_FILE_DIR:${PROJECT_NAME}>/config)
