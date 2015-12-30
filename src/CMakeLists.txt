set(ekats_source_files
        GramMatrix.cpp
	DistanceL1.cpp
        UniformTimeseries.cpp
	DistanceL2.cpp
	DTW.cpp
	KDTW.cpp
	kmedoids.cpp
	kcentroids.cpp
	cluster.cpp
	DBA.cpp
	utils.cpp
	matrix.h

        DistanceL2.h
  	DTW.h
	KDTW.h
        Point.h
	DistanceL1.h
      	GramMatrix.h
        UniformTimeseries.h
	kmedoids.h
	kcentroids.h
	cluster.h
	DBA.h
	utils.h
)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${Project_HOME_DIR}/lib")

add_library(ekats ${ekats_source_files})
target_link_libraries(ekats)
