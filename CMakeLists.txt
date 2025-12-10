###################
##### Issue the following cmake command from ./build/ for Xerxes (GDML support). Drop all flags if GDML not needed.
## cmake -DUSEG4MT=1 -DXercesC_INCLUDE_DIR=/Users/chur558/xerces-c-3.3.0/src -DXercesC_LIBRARY=/Users/chur558/xerces-c-3.3.0/build/src/libxerces-c-3.3.dylib -DXercesC_VERSION=3.3.0 ../
###################

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(sndwrm)

# Required by MARLEY
set(CMAKE_CXX_STANDARD 14)


#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED gdml)
endif()

message(STATUS "GEANT4 include dir: ${Geant4_INCLUDE_DIR}")
message(STATUS "GEANT4 libraries dir: ${Geant4_LIBRARIES_DIR}")
message(STATUS "GEANT4 libraries libs: ${Geant4_LIBRARIES}")


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})
#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR} $ENV{MARLEY}/include/
		    $ENV{XERCESROOT}/../../src/
		    )
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
##file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cxx)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
##file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
# Add -DUSEG4MT=1 to cmake options to use multithreading
option(USEG4MT "Compile with G4MULTITHREADing " OFF)
if (USEG4MT)
   add_compile_definitions(USEG4MT)
endif()

add_executable(sndwrm rdecay02.cc ${sources} ${headers})
#target_link_libraries(sndwrm  ${Geant4_LIBRARIES}  $ENV{MARLEY}/build/libMARLEY.so $ENV{XERCESROOT}/libxerces-c-3.3.so) ## .so for linux
# Find Xerces-C library
find_package(XercesC REQUIRED)

# Link libraries
target_link_libraries(sndwrm
    ${Geant4_LIBRARIES}
    $ENV{MARLEY}/build/libMARLEY.so
    ${XercesC_LIBRARIES}
)

# Include directories for Xerces-C
target_include_directories(sndwrm PRIVATE ${XercesC_INCLUDE_DIRS})

message(STATUS "USEG4MT is ${USEG4MT}")

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build rdecay02. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

##set(rdecay02_SCRIPTS
##    debug.mac
##    plotHisto.C
##    rdecay02.in
##    vis.mac
##  )

##foreach(_script ${rdecay02_SCRIPTS})
##  configure_file(
##    ${PROJECT_SOURCE_DIR}/${_script}
##    ${PROJECT_BINARY_DIR}/${_script}
##    COPYONLY
##    )
##endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS sndwrm DESTINATION bin)

