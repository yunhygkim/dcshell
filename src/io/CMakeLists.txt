######## common Libraries ########


set(LIB_SRC
	CD3d_msg.cpp
	CD3dexport.cpp
)

include_directories(
	"${CMAKE_SOURCE_DIR}/src"
	"${CMAKE_SOURCE_DIR}/io"
	"${CMAKE_SOURCE_DIR}/lib/qhull"
	)

add_library(io ${LIB_SRC})

