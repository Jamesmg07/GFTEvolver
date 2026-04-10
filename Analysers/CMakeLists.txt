# Include path
set(HEADERS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Header")
set(ROOT_HEADERS_DIR "${CMAKE_SOURCE_DIR}/Core/Header")

# -----------------------------
# Sources
# -----------------------------
set(
  SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/Analyser.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/OutputFields.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/Energy.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/GaugeCondition.cpp"
)

# Create a library
add_library(Analysers ${SOURCES})

# Expose this folder for includes (PUBLIC so dependents see it)
target_include_directories(Analysers
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Create macro for the source directory
target_compile_definitions(Analysers PRIVATE
    SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
    DATA_DIR="${CMAKE_SOURCE_DIR}/Data"
)

# -----------------------------
# Header files set-up
# -----------------------------

target_include_directories(Analysers PUBLIC ${HEADERS_DIR} ${ROOT_HEADERS_DIR})

# (Optional) If module depends on other libraries
target_link_libraries(Analysers PUBLIC Potentials Gradients WilsonLoops)
