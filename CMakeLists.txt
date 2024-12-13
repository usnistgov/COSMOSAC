project(cCOSMO)
cmake_minimum_required(VERSION 3.17)

# Set the C++ standard to C++17
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

string(CONFIGURE "#define COSMO_SAC_HOME \"${CMAKE_CURRENT_SOURCE_DIR}\"" CONFIG_HEADER)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/HOME_PATH.h" ${CONFIG_HEADER})

set(MY_INCLUDES
    "${CMAKE_CURRENT_SOURCE_DIR}/include" 
    "${CMAKE_CURRENT_SOURCE_DIR}/externals/Eigen" 
    "${CMAKE_CURRENT_BINARY_DIR}" 
    "${CMAKE_CURRENT_SOURCE_DIR}/externals/nlohmann/single_include"
)

add_executable(main "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
target_sources(main PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/externals/Eigen/debug/msvc/eigen.natvis")
target_include_directories(main PUBLIC ${MY_INCLUDES})

add_subdirectory(externals/pybind11)
pybind11_add_module(cCOSMO "${CMAKE_CURRENT_SOURCE_DIR}/src/pybind11_interface.cxx")
set_target_properties(cCOSMO PROPERTIES COMPILE_FLAGS -DPYBIND11)
target_include_directories(cCOSMO PUBLIC ${MY_INCLUDES})

add_executable(catch_tests "${CMAKE_CURRENT_SOURCE_DIR}/tests/catch_tests.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/tests/catch_main.cpp")
target_sources(catch_tests PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/externals/Eigen/debug/msvc/eigen.natvis" "${CMAKE_CURRENT_SOURCE_DIR}/externals/nlohmann/nlohmann_json.natvis")
target_include_directories(catch_tests PUBLIC ${MY_INCLUDES} "${CMAKE_CURRENT_SOURCE_DIR}/externals/Catch/single_include" )
