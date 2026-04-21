cmake_minimum_required(VERSION 3.20)

project(Molecule-Visualization LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#edit as needed
find_package(VTK COMPONENTS 
  CommonColor
  CommonCore
  DomainsChemistry
  DomainsChemistryOpenGL2
  IOChemistry
  InteractionStyle
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "Molecule-Visualization: Unable to find the VTK build folder.")
endif()


find_package(nlohmann_json REQUIRED)
find_package(HDF5 REQUIRED)
find_package(HighFive REQUIRED)
find_package(Armadillo REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(chemfiles REQUIRED)

add_library(common INTERFACE)

target_link_libraries(common INTERFACE
    Eigen3::Eigen
    armadillo
    nlohmann_json::nlohmann_json
    HDF5::HDF5
    chemfiles
)
target_link_libraries(common INTERFACE ${VTK_LIBRARIES})

target_compile_definitions(common INTERFACE 
    ARMA_USE_HDF5
    H5_USE_110_API
)

link_libraries(common)


add_executable(
    Molecule-Visualization
    src/main.cpp
    src/molecule_functions.cpp
    src/fock_functions.cpp
    src/sto3g_functions.cpp
    src/visualize.cpp
)

vtk_module_autoinit(
  TARGETS Molecule-Visualization
  MODULES ${VTK_LIBRARIES}
)
