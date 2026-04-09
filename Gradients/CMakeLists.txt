# Include path
set(HEADERS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Header")

# -----------------------------
# Sources
# -----------------------------
set(
  SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/Gradient.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/Global.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/StandardModel.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Source/TwoHDM.cpp"
)

# Create a library
add_library(Gradients ${SOURCES})

# Expose this folder for includes (PUBLIC so dependents see it)
target_include_directories(Gradients
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Create macro for the source directory
target_compile_definitions(Gradients PRIVATE
    SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
)

# -----------------------------
# Header files set-up
# -----------------------------

target_include_directories(Gradients PUBLIC ${HEADERS_DIR})

# (Optional) If module depends on other libraries
#target_link_libraries(this_library PUBLIC other libraries)
