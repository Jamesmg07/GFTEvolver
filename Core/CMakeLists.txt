# Include path
set(HEADERS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Header")

# -----------------------------
# Sources
# -----------------------------
set(
  SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/Lattice.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/Model.cpp"
)

# Create a library
add_library(Core ${SOURCES})

# Expose this folder for includes (PUBLIC so dependents see it)
target_include_directories(Core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Create macro for the source directory
target_compile_definitions(Core PRIVATE
    SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
)

if (GFT_ENABLE_MPI)
    target_compile_definitions(Core PRIVATE GFT_ENABLE_MPI)
    target_link_libraries(Core PUBLIC MPI::MPI_CXX)
endif()

# -----------------------------
# Header files set-up
# -----------------------------

target_include_directories(Core PUBLIC ${HEADERS_DIR})

# (Optional) If module depends on other libraries
target_link_libraries(Core PUBLIC 
    Potentials 
    Gradients 
    WilsonLoops 
    InitialConditions 
    BoundaryConditions 
    Analysers
)
