cmake_minimum_required(VERSION "3.10")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#Enable all standard build types
set(CMAKE_CONFIGURATION_TYPES "Debug;Release;RelWithDebInfo;MinSizeRel" CACHE STRING "" FORCE)

# Default to Release build if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build (Debug, Release, RelWithDebInfo, MinSizeRel)." FORCE)
endif()

project(GFTEvolver)

# MPI support is optional so that the original serial build remains available
# without requiring an MPI installation.
option(GFT_ENABLE_MPI "Enable MPI parallelisation" OFF)

if (GFT_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
endif()

# Include path
set(PROJECT_HEADERS_DIR "${CMAKE_SOURCE_DIR}/Header")


# -----------------------------
# Sources
# -----------------------------
#set(
  #SOURCES
  #"Source/main.cpp"
  #"Source/Lattice.cpp"
  #"Source/Model.cpp"
#)

add_executable(GFTEvolver main.cpp)
if (GFT_ENABLE_MPI)
    target_compile_definitions(GFTEvolver PRIVATE GFT_ENABLE_MPI)
endif()

if (MSVC)
    target_compile_options(GFTEvolver PRIVATE
        $<$<CONFIG:Debug>:/W4 /permissive- /Od>
        $<$<CONFIG:Release>:/W4 /permissive- /O3>
    )
else()
    target_compile_options(GFTEvolver PRIVATE
        $<$<CONFIG:Debug>: -Wall -Wextra -Wpedantic -Wno-ignored-qualifiers -Wfloat-conversion -Wsign-compare -Wconversion -O0 -g>
        $<$<CONFIG:Release>: -O3>
        $<$<CONFIG:RelWithDebInfo>: -Wall -Wextra -Wpedantic -Wno-ignored-qualifiers -Wsign-compare -Wconversion -O3>
    )
endif()

# Create macro for the source directory
target_compile_definitions(GFTEvolver PRIVATE
    SOURCE_DIR="${CMAKE_SOURCE_DIR}"
    DATA_DIR="${CMAKE_SOURCE_DIR}/Data"
)


# -----------------------------
# Subdirectories (modules)
# -----------------------------

add_subdirectory(Core)
add_subdirectory(Potentials)
add_subdirectory(Gradients)
add_subdirectory(WilsonLoops)
add_subdirectory(InitialConditions)
add_subdirectory(BoundaryConditions)
add_subdirectory(Analysers)

# -----------------------------
# Header files set-up
# -----------------------------

target_include_directories(GFTEvolver PRIVATE ${PROJECT_HEADERS_DIR})

# -----------------------------
# Linking
# -----------------------------

target_link_libraries(GFTEvolver PRIVATE
    Core
    Potentials
    Gradients
    WilsonLoops
    InitialConditions
    BoundaryConditions
    Analysers)

if (GFT_ENABLE_MPI)
    target_link_libraries(GFTEvolver PRIVATE MPI::MPI_CXX)
endif()