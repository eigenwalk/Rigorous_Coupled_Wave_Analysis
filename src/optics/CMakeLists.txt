SET(rcwa_source
    rcwaRun.cpp
	light.cpp
    scatterMatrix.cpp
)

add_library(rcwa_lib 
    ${rcwa_source}
)
