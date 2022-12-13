cmake_minimum_required(VERSION 3.14)

project(Sandbox LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 COMPONENTS Core REQUIRED)
find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Core REQUIRED)

add_executable(Sandbox
  mytimer.cpp
  mytimer.hpp

  main.cpp
  simulation.cpp
  simulation.hpp
  ipotential.hpp
  lj.cpp
  lj.hpp
  dz.cpp
  dz.hpp
)
target_link_libraries(Sandbox Qt${QT_VERSION_MAJOR}::Core)
