######## common Libraries ########

set(LIB_SRC
	CD3d_hull.cpp
	hull_vol.cpp
)

include_directories(
	"${CMAKE_SOURCE_DIR}/src"
	"${CMAKE_SOURCE_DIR}/src/measure"
	"${CMAKE_SOURCE_DIR}/src/util"
	"${CMAKE_SOURCE_DIR}/lib/qhull"
	"${CMAKE_SOURCE_DIR}/lib/qhull"
	"${CMAKE_SOURCE_DIR}/lib"
	)


add_library(measure ${LIB_SRC})

