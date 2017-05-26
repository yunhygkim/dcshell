######## common Libraries ########

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(CGAL REQUIRED)


include_directories(
	${OPENGL_INCLUDE_DIRS}
	${GLUT_INCLUDE_DIRS}
	${CGAL_INCLUDE_DIRS}
	)

link_directories(
	"${CMAKE_SOURCE_DIR}/lib"
	"${CMAKE_SOURCE_DIR}/lib/Eigen"
	"${CMAKE_SOURCE_DIR}/lib/FIST"
	"${CMAKE_SOURCE_DIR}/lib/GL"
	"${CMAKE_SOURCE_DIR}/lib/glew"
	"${CMAKE_SOURCE_DIR}/lib/glew/lib"
	"${CMAKE_SOURCE_DIR}/lib/glpk"
	"${CMAKE_SOURCE_DIR}/lib/gsl"
	"${CMAKE_SOURCE_DIR}/lib/gsl/lib"
	"${CMAKE_SOURCE_DIR}/lib/qhull"
	"${CMAKE_SOURCE_DIR}/lib/RAPID"
	)

set(LIBCSHELL_DEPS
	mathtool
	io
	util
	gui
	GL
	measure
	${OpenCV_LIBS}
	${GLPK_LIBRARY}
	${CGAL_LIBRARY}
	${GMP_LIBRARY}	
	)

if(APPLE)
	set(LIBCSHELL_DEPS
		${LIBCSHELL_DEPS}
		libGLEW.osx.a
		libFIST.osx.a
		libglpk.osx.a
		libqhull.osx.a
		libRAPID.osx.a
		)
endif()

if(WIN32)
	set(LIB_DEPS
		${LIB_DEPS}
		RAPID.lib
		glut32.lib
		glpk.lib
	)
endif()

######## Add subdirectories########

add_subdirectory(gui)
add_subdirectory(io)
add_subdirectory(measure)
add_subdirectory(util)

######## cshell libraries ########

set(LIBCSHELL_SRC
	${LIBCSHELL_SRC}
	CD3d_model.cpp
	cshell.cpp
	)

add_library(libcshell STATIC 
	${LIBCSHELL_SRC}
)
target_link_libraries(libcshell ${LIBCSHELL_DEPS})

######## cshell GUI ########
set(APP_NAME cshell)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})


include_directories(
	"${CMAKE_SOURCE_DIR}/src"
	"${CMAKE_SOURCE_DIR}/src/gui"
	"${CMAKE_SOURCE_DIR}/src/io"
	"${CMAKE_SOURCE_DIR}/src/measure"
	"${CMAKE_SOURCE_DIR}/src/util"
	"${CMAKE_SOURCE_DIR}/lib"
	"${CMAKE_SOURCE_DIR}/lib/qhull"
	)

link_directories(
	"${CMAKE_SOURCE_DIR}/lib"
	)

file(GLOB LIBRARIES libCGAL_Core.12.0.0.dylib
	libCGAL.12.0.0.dylib
	libboost_system-mt.dylib
	libboost_thread-mt.dylib
	"/opt/local/lib/libgmp*.dylib"
	libmpfr.4.dylib
)

file(GLOB LIBRARIES2 "${CMAKE_SOURCE_DIR}/lib/GL/*.a")

set(CSHELL_SRC
	${CSHELL_SRC}
	cshell.cpp
	)

add_executable(${APP_NAME} ${CSHELL_SRC})
add_dependencies(${APP_NAME} libcshell)


target_link_libraries(${APP_NAME}
	${LIBCSHELL_DEPS}
	libcshell
	${OPENGL_LIBRARIES}
	${GLUT_LIBRARY}
	${CGAL_LIBRARY}
	${GMP_LIBRARY}
	${LIBRARIES}
	${LIBRARIES2}
	)
