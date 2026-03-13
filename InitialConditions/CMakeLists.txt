# Include path
set(HEADERS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Header")

# -----------------------------
# Sources
# -----------------------------
set(
  SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/InitialCondition.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/RandomUniform.cpp"
)

# Create a library
add_library(InitialConditions ${SOURCES})

# Expose this folder for includes (PUBLIC so dependents see it)
target_include_directories(InitialConditions
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Create macro for the source directory
target_compile_definitions(InitialConditions PRIVATE
    SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
)

# -----------------------------
# Header files set-up
# -----------------------------

target_include_directories(InitialConditions PUBLIC ${HEADERS_DIR})

# (Optional) If module depends on other libraries
#target_link_libraries(this_library PUBLIC other libraries)
