# Include path
set(HEADERS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Header")
set(ROOT_HEADERS_DIR "${CMAKE_SOURCE_DIR}/Header")

# -----------------------------
# Sources
# -----------------------------
set(
  SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/BoundaryCondition.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/Periodic.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/Interior.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/Fixed.cpp"
)

# Create a library
add_library(BoundaryConditions ${SOURCES})

# Expose this folder for includes (PUBLIC so dependents see it)
target_include_directories(BoundaryConditions
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Create macro for the source directory
target_compile_definitions(BoundaryConditions PRIVATE
    SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
)

# -----------------------------
# Header files set-up
# -----------------------------

target_include_directories(BoundaryConditions PUBLIC ${HEADERS_DIR} ${ROOT_HEADERS_DIR})

# (Optional) If module depends on other libraries
target_link_libraries(BoundaryConditions PUBLIC Potentials Gradients WilsonLoops Analysers)
