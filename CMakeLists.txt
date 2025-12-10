cmake_minimum_required(VERSION 3.16)

# Sets project name
project(projects_relativistic LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

# Sets C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specifies required Qt components
find_package(Qt6 REQUIRED COMPONENTS Concurrent)
find_package(Qt6 REQUIRED COMPONENTS Core)
find_package(Qt6 REQUIRED COMPONENTS Gui)
find_package(Qt6 REQUIRED COMPONENTS Xml)

# Allows you to include files from within those directories, without prefixing their filepaths
include_directories(src)

add_definitions(-DGLM_FORCE_SWIZZLE)

# Specifies .cpp and .h files to be passed to the compiler
add_executable(${PROJECT_NAME}
  src/main.cpp
  src/raytrace.h
  src/utils/rgba.h
  src/raytrace.cpp
  src/utils/scenefileparser.cpp
  src/utils/scenefileparser.h
  src/camerapath.cpp
  src/camerapath.h
  src/renderer.h src/renderer.cpp
  src/lighting.h src/lighting.cpp
  src/textures.h src/textures.cpp
  src/framedata.h src/framedata.cpp
  src/utils/structs.h
  src/utils/imagereader.h src/utils/imagereader.cpp
)

# GLM: this creates its library and allows you to `#include "glm/..."`
add_subdirectory(glm)

target_link_libraries(${PROJECT_NAME} PRIVATE
    Qt::Concurrent
    Qt::Core
    Qt::Gui
    Qt::Xml
)

# Set this flag to silence warnings on Windows
if (MSVC OR MSYS OR MINGW)
  set(CMAKE_CXX_FLAGS "-Wno-volatile")
endif()

# Set this flag to silence warnings on MacOS
if (APPLE)
  set(CMAKE_CXX_FLAGS "-Wno-deprecated-volatile")
endif()
