include (plugin)
SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} /usr/local/include)

SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/local/lib)
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/local/qwt-6.1.1/lib)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(Boost 1.55 REQUIRED serialization)

openflipper_plugin (DEPS GMP CGAL EIGEN3 Boost
					DIRS SuperDeform
					INCDIRS SuperDeform)