cmake_minimum_required(VERSION 3.22)

project(
  paraPhysiCell
  VERSION 0.1
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

if(BUILD_TESTS)
  # Fetch Googletest
  FetchContent_Declare(
    googletest
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    URL file://${CMAKE_CURRENT_SOURCE_DIR}/submodules/paraBioFVM/submodules/googletest
  )
  # For Windows: Prevent overriding the parent project's compiler/linker
  # settings
  set(gtest_force_shared_crt
      ON
      CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endif()

# Fetch BioFVM
FetchContent_Declare(
  paraBioFVM
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL file://${CMAKE_CURRENT_SOURCE_DIR}/submodules/paraBioFVM)
FetchContent_MakeAvailable(paraBioFVM)

# Fetch pugixml
FetchContent_Declare(
  pugixml
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL file://${CMAKE_CURRENT_SOURCE_DIR}/submodules/pugixml)
FetchContent_MakeAvailable(pugixml)

if(BUILD_FOR_DEVICE AND DEVICE_IMPL STREQUAL "cuda")
  enable_language(CUDA)
endif()

# Target PhysiCellCore
file(GLOB_RECURSE src_files "src/*")
list(FILTER src_files EXCLUDE REGEX ".*main\\..*")

add_library(PhysiCellCore ${src_files})

target_include_directories(PhysiCellCore PUBLIC ${paraBioFVM_SOURCE_DIR}/include
                                                ${pugixml_SOURCE_DIR}/src)

target_link_libraries(PhysiCellCore PUBLIC BioFVMCore pugixml-static)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(PhysiCellCore PUBLIC OpenMP::OpenMP_CXX)
endif()

# Target PhysiCell
add_executable(PhysiCell src/main.cpp)
target_link_libraries(PhysiCell PhysiCellCore)

if(BUILD_TESTS)
  # Target unit_PhysiCell
  file(GLOB_RECURSE test_files "test/*")
  add_executable(unit_PhysiCell ${test_files})

  target_link_libraries(unit_PhysiCell PRIVATE gtest_main gmock_main
                                               PhysiCellCore)

  target_include_directories(unit_PhysiCell PUBLIC "src")

  enable_testing()
  include(GoogleTest)
  gtest_discover_tests(unit_PhysiCell)
endif()

add_subdirectory(sample_projects/biorobots)
add_subdirectory(sample_projects/heterogeneity)
add_subdirectory(sample_projects/interactions)
add_subdirectory(sample_projects/template)
add_subdirectory(sample_projects/sem-test)
