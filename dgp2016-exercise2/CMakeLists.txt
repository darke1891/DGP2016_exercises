# Set the minimum required CMake version
cmake_minimum_required(VERSION 3.0)

# Guess project name from folder name 
get_filename_component(PROJECTNAME ${CMAKE_CURRENT_LIST_DIR} NAME)
project(${PROJECTNAME}) #< note: this wipes CMAKE_MODULE_PATH

# Directory of cmake config files
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#================================
# Find required packages
#================================
set(CMAKE_FIND_QUIETLY TRUE)
include(cmake/ConfigureEigen.cmake)
include(cmake/ConfigureOpenGP.cmake)
include(cmake/ConfigureCompiler.cmake)
include(cmake/ConfigureGLEW.cmake)
include(cmake/ConfigureOpenGL.cmake)
include(cmake/ConfigureGLFW3.cmake)

#================================
# Subdirectories
#================================
add_subdirectory(curves)