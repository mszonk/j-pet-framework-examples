################################################################################
<<<<<<< HEAD
## Data analysis project based on J-PET Framework
## Created by J-PET Framework developers 2016-2018
##
## Description:
##   Builds template for new analysis project
=======
## Set of projects based on J-PET Framework
## Created by J-PET Framework developers 2016-2018
##
## Description:
##   Script for building J-PET Framework, J-PET MLEM and all the examples
>>>>>>> b76e8c289e91932a2e7bf30b19ab6fcaaa8d2d1e
################################################################################

cmake_minimum_required(VERSION 2.6)

<<<<<<< HEAD
######################################################################
## Project name - set this section appropriately for your project
## (please remove this comment in your template-derived analysis)
######################################################################
set(projectName pPsto3gamma)

## Auxiliary files are copied to build directory after CMake
## You can add any files to this list
set(AUXILLIARY_FILES
  README.md
)

set(ROOT_SCRIPTS
  rootlogon.C
  rootalias.C
)

################################################################################
## Binary, header and source files definitions
set(projectBinary ${projectName}.x)
project(${projectName} CXX) # using only C++

## Name your header files with .h extension and source files with .cpp
file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

######################################################################
## If you are using your own version of some module which is also present
## in LargeBarrelAnalysis, comment out the respective lines below
## so as not to use this module's source from LargeBarrelAnalysis.
## (please remove this comment in your template-derived analysis)
set(use_modules_from ../LargeBarrelAnalysis)
list(APPEND HEADERS ${use_modules_from}/TimeWindowCreator.h)
list(APPEND SOURCES ${use_modules_from}/TimeWindowCreator.cpp)
list(APPEND HEADERS ${use_modules_from}/TimeWindowCreatorTools.h)
list(APPEND SOURCES ${use_modules_from}/TimeWindowCreatorTools.cpp)
list(APPEND HEADERS ${use_modules_from}/UniversalFileLoader.h)
list(APPEND SOURCES ${use_modules_from}/UniversalFileLoader.cpp)
list(APPEND HEADERS ${use_modules_from}/SignalFinder.h)
list(APPEND SOURCES ${use_modules_from}/SignalFinder.cpp)
list(APPEND HEADERS ${use_modules_from}/SignalFinderTools.h)
list(APPEND SOURCES ${use_modules_from}/SignalFinderTools.cpp)
list(APPEND HEADERS ${use_modules_from}/SignalTransformer.h)
list(APPEND SOURCES ${use_modules_from}/SignalTransformer.cpp)
list(APPEND HEADERS ${use_modules_from}/HitFinder.h)
list(APPEND SOURCES ${use_modules_from}/HitFinder.cpp)
list(APPEND HEADERS ${use_modules_from}/HitFinderTools.h)
list(APPEND SOURCES ${use_modules_from}/HitFinderTools.cpp)
list(APPEND HEADERS ${use_modules_from}/EventFinder.h)
list(APPEND SOURCES ${use_modules_from}/EventFinder.cpp)
list(APPEND HEADERS ${use_modules_from}/EventCategorizer.h)
list(APPEND SOURCES ${use_modules_from}/EventCategorizer.cpp)
list(APPEND HEADERS ${use_modules_from}/EventCategorizerTools.h)
list(APPEND SOURCES ${use_modules_from}/EventCategorizerTools.cpp)

include_directories(${Framework_INCLUDE_DIRS})
add_definitions(${Framework_DEFINITIONS})

add_executable(${projectBinary} ${SOURCES} ${HEADERS})
target_link_libraries(${projectBinary} JPetFramework)

add_custom_target(clean_data_${projectName}
  COMMAND rm -f *.tslot.*.root *.phys.*.root *.sig.root
)

## Create variable for list with depends files path
set(copy_depends)

################################################################################
## Copy the example auxiliary files
foreach(file_i ${AUXILLIARY_FILES})
  list(APPEND copy_depends ${CMAKE_CURRENT_BINARY_DIR}/${file_i})
  if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${file_i})
    set(CP_CMD copy_directory)
  else()
    set(CP_CMD copy)
  endif()
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${file_i}
    POST_BUILD COMMAND ${CMAKE_COMMAND}
    ARGS -E ${CP_CMD} ${CMAKE_CURRENT_SOURCE_DIR}/${file_i} ${CMAKE_CURRENT_BINARY_DIR}/${file_i}
  )
endforeach(file_i)

## Copy the ROOT scripts
foreach(file_i ${ROOT_SCRIPTS})
  list(APPEND copy_depends ${CMAKE_CURRENT_BINARY_DIR}/${file_i})
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${file_i}
    POST_BUILD COMMAND ${CMAKE_COMMAND}
    ARGS -E ${CP_CMD} ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/${file_i} ${CMAKE_CURRENT_BINARY_DIR}/${file_i}
  )
endforeach(file_i)

################################################################################
## Add new target that depends on copied files
add_custom_target(copy_files_${projectName} DEPENDS ${copy_depends})
## Add dependencies between building project and copy files
add_dependencies(${projectBinary} copy_files_${projectName})
=======
################################################################################
## Add Framework submodule
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/j-pet-framework/cmake)
add_subdirectory(j-pet-framework)
include_directories(j-pet-framework ${Framework_INCLUDE_DIRS})
add_definitions(${Framework_DEFINITIONS})

################################################################################
# Add MLEM submodule
add_subdirectory(j-pet-mlem)
#SYSTEM keyword is a trick to disable all warnings from the j-pet-mlem package.
#This code invokes some horrible warnings that we deliberately turn off, to
#keep our build process clean.
include_directories(SYSTEM j-pet-mlem ${MLEM_INCLUDE_DIRS})

################################################################################
add_definitions(-std=c++11 -Wall -Wunused-parameter)

################################################################################
## Enable rpath on OS X and point it to ROOT
if(APPLE)
  set(CMAKE_MACOSX_RPATH ON)
  set(CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -std=c++11 -Wl,-rpath,${ROOT_LIBRARY_DIR}"
    )
endif()

################################################################################
## Download input and configuration files
## The script shouldn't do anything if the data is present and correct.
set(DOWNLOAD_MORE_DATA ${CMAKE_CURRENT_SOURCE_DIR}/download_reformed_data.sh ${CMAKE_CURRENT_SOURCE_DIR})
execute_process(COMMAND ${DOWNLOAD_MORE_DATA})

################################################################################
## Directories with examples
add_subdirectory(InterThresholdCalibration)
add_subdirectory(LargeBarrelAnalysis)
add_subdirectory(VelocityCalibration)
add_subdirectory(NewAnalysisTemplate)
add_subdirectory(ImageReconstruction)
add_subdirectory(ScopeLoaderExample)
add_subdirectory(TimeCalibration)
add_subdirectory(MCGeantAnalysis)
add_subdirectory(PhysicAnalysis)
add_subdirectory(CosmicAnalysis)
add_subdirectory(ScopeAnalysis)
add_subdirectory(Imaging)
add_subdirectory(modules)
>>>>>>> b76e8c289e91932a2e7bf30b19ab6fcaaa8d2d1e
