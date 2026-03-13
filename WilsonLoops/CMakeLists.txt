# Include path
set(HEADERS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Header")

# -----------------------------
# Sources
# -----------------------------
set(
  SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/WilsonLoop.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/StandardModel.cpp"
)

# Create a library
add_library(WilsonLoops ${SOURCES})

# Expose this folder for includes (PUBLIC so dependents see it)
target_include_directories(WilsonLoops
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Create macro for the source directory
target_compile_definitions(WilsonLoops PRIVATE
    SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
)

# -----------------------------
# Header files set-up
# -----------------------------

target_include_directories(WilsonLoops PUBLIC ${HEADERS_DIR})

# (Optional) If module depends on other libraries
#target_link_libraries(this_library PUBLIC other libraries)
