######## common Libraries ########

set(LIB_SRC
	Basic.cpp
	Gauss.cpp
	Quaternion.cpp
)

include_directories(
	"${CMAKE_SOURCE_DIR}/src"
	)

add_library(mathtool ${LIB_SRC})

