cmake_minimum_required(VERSION 3.18)

project(
  sem-intrusion
  VERSION 0.1
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE src_files "*")
add_executable(sem-intrusion ${src_files})
target_include_directories(sem-intrusion PRIVATE "../..")
target_link_libraries(sem-intrusion PRIVATE PhysiCellCore)

add_custom_command(
  TARGET sem-intrusion
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/config
          $<TARGET_FILE_DIR:${PROJECT_NAME}>/config)
