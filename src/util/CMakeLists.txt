######## common Libraries ########



set(LIB_SRC
	${LIBCSHELL_SRC}
	CD.cpp
	CD3d_convexity.cpp
	CD3d_hspace_intersection.cpp
	CD3d_progressive_hull.cpp
	CD3d_trim_hull.cpp
	CD3d_util.cpp
	CD3d_vsa.cpp
	svm.cpp
	)

include_directories(
	"${CMAKE_SOURCE_DIR}/src"
	"${CMAKE_SOURCE_DIR}/util"
	"${CMAKE_SOURCE_DIR}/measure"
	"${CMAKE_SOURCE_DIR}/src/measure"
	"${CMAKE_SOURCE_DIR}/lib/qhull"	
	"${CMAKE_SOURCE_DIR}/lib"	
	"${CMAKE_SOURCE_DIR}/lib/RAPID"
	"${CMAKE_SOURCE_DIR}/lib/glpk"
	)

add_library(util ${LIB_SRC})

add_subdirectory(mathtool)
