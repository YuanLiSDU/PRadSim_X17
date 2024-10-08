#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(PRadSim)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -luuid")

add_compile_definitions(G4VIS_USE)
add_compile_definitions(G4UI_USE)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Tell CMake where to find the ROOT installation
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED COMPONENTS MathMore Foam)

#----------------------------------------------------------------------------
# Setup include and lib directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB PRADSIM_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB PRADSIM_HEADERS ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the libraries
#
add_executable(PRadSim PRadSim.cc ${PRADSIM_SOURCES} ${PRADSIM_HEADERS})
target_link_libraries(PRadSim ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Add subdirectories to build
#
add_subdirectory(digitization)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory
# This is so that we can run the executable directly because it relies on these scripts being in the current working directory
#
set(PRADSIM_RESOURCE_FILES init_vis.mac vis.mac gui.mac run.mac)

add_custom_target(copy_files ALL COMMENT "Copy files")
foreach(_file ${PRADSIM_RESOURCE_FILES})
  add_custom_command(TARGET copy_files POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/${_file} ${PROJECT_BINARY_DIR}/${_file}
  )
endforeach()

