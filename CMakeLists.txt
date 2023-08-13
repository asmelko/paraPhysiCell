cmake_minimum_required(VERSION 3.18)

project(paraPhysiCell VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Fetch Googletest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Fetch BioFVM
FetchContent_Declare(
  paraBioFVM
  GIT_REPOSITORY https://github.com/asmelko/paraBioFVM.git
  GIT_TAG        omp
)
FetchContent_MakeAvailable(paraBioFVM)

# Fetch pugixml
FetchContent_Declare(
  pugixml
  URL https://github.com/zeux/pugixml/archive/refs/tags/v1.13.zip
)
FetchContent_MakeAvailable(pugixml)

# Target PhysiCellCore
file(GLOB_RECURSE src_files "src/*")
list(FILTER src_files EXCLUDE REGEX ".*main\\..*")

add_library(PhysiCellCore ${src_files})

target_include_directories(PhysiCellCore PUBLIC ${noarr_structures_SOURCE_DIR}/include ${paraBioFVM_SOURCE_DIR}/include ${pugixml_SOURCE_DIR}/src)

target_link_libraries(PhysiCellCore PUBLIC BioFVMCore pugixml-static)

# Target unit_PhysiCell
file(GLOB_RECURSE test_files "test/*")
add_executable(unit_PhysiCell ${test_files})

# Target PhysiCell
add_executable(PhysiCell src/main.cpp)
target_link_libraries(PhysiCell PhysiCellCore)

target_link_libraries(unit_PhysiCell PRIVATE
	gtest_main gmock_main PhysiCellCore)

target_include_directories(unit_PhysiCell PUBLIC
	"src")

foreach(target PhysiCell unit_PhysiCell PhysiCellCore)
	if(MSVC)
		target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/W4 /bigobj>)
	else()
		target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -pedantic>)
	endif()
endforeach()

enable_testing()
include(GoogleTest)
gtest_discover_tests(unit_PhysiCell)

add_subdirectory(sample_projects/biorobots)
add_subdirectory(sample_projects/heterogeneity)
add_subdirectory(sample_projects/template)
