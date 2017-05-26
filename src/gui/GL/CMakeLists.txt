######## common Libraries ########


set(LIB_SRC
	gli.cpp
	gliDump.cpp
	gliFont.cpp
	shadowMap_glsl.cpp
)

include_directories(
	"${CMAKE_SOURCE_DIR}/lib"
	"${CMAKE_SOURCE_DIR}/src/gui"
	)



add_library(GL ${LIB_SRC})
