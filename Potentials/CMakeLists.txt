# Include path
set(HEADERS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Header")

# -----------------------------
# Sources
# -----------------------------
set(
  SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/Potential.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/SO_N.cpp"
)

# Create a library
add_library(Potentials ${SOURCES})

# Expose this folder for includes (PUBLIC so dependents see it)
target_include_directories(Potentials
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Create macro for the source directory
target_compile_definitions(Potentials PRIVATE
    SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
)

# -----------------------------
# Header files set-up
# -----------------------------

target_include_directories(Potentials PUBLIC ${HEADERS_DIR})

# (Optional) If module depends on other libraries
#target_link_libraries(this_library PUBLIC other libraries)
