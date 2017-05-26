cmake_minimum_required(VERSION 2.6)

project("cshell")

set(PACKAGE_VERSION "0.1.0")

####### ADD PACKAGES #######
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(CGAL REQUIRED)

####### SET COMPILER FLAGS #######
if(NOT WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
endif()

if(WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /d _CRT_SECURE_NO_WARNINGS")
endif()

####### ADD SUBDIRECTORIES #######
add_subdirectory(src)

