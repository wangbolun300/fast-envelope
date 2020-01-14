# Install script for directory: /home/bw1760/rational/fast-envelope/3rdparty/geogram/doc

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xruntimex" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/doc/geogram" TYPE FILE OPTIONAL FILES "/home/bw1760/rational/fast-envelope/debug/geogram/doc/VERSION.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdoc-devkitx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/doc/devkit" TYPE DIRECTORY OPTIONAL FILES "/home/bw1760/rational/fast-envelope/debug/geogram/doc/devkit/html")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdoc-devkit-fullx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/doc/devkit" TYPE DIRECTORY OPTIONAL FILES "/home/bw1760/rational/fast-envelope/debug/geogram/doc/devkit-full/html")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdoc-devkit-internalx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/doc/devkit" TYPE DIRECTORY OPTIONAL FILES "/home/bw1760/rational/fast-envelope/debug/geogram/doc/devkit-internal/html")
endif()

