cmake_minimum_required(VERSION 3.18)

project(
  sem-contact
  VERSION 0.1
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE src_files "*")
add_executable(sem-contact ${src_files} ${CMAKE_CURRENT_SOURCE_DIR}/../helpers.cpp)
target_include_directories(sem-contact PRIVATE "../..")
target_link_libraries(sem-contact PRIVATE PhysiCellCore)

add_custom_command(
  TARGET sem-contact
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/config
          $<TARGET_FILE_DIR:${PROJECT_NAME}>/config)
