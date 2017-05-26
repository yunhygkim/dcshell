######## common Libraries ########


set(LIB_SRC
	CD3d_draw_GL.cpp
	CD3d_stat.cpp
	dude3d_draw.cpp
)

include_directories(
	"${CMAKE_SOURCE_DIR}/src/gui"
	"${CMAKE_SOURCE_DIR}/lib"
	"${CMAKE_SOURCE_DIR}/src"
	"${CMAKE_SOURCE_DIR}/src/measure"
	"${CMAKE_SOURCE_DIR}/lib/qhull"
	"${CMAKE_SOURCE_DIR}/lib/RAPID"
	"${CMAKE_SOURCE_DIR}/lib/GL"
	)



add_library(gui ${LIB_SRC})

add_subdirectory(GL)
