include (plugin)
set(CGAL_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/../CGAL-4.5/include)
set(CGAL_LIBRARY_DIR ${CMAKE_SOURCE_DIR}/../CGAL-4.5/build/lib)
openflipper_plugin (DEPS GMP CGAL EIGEN3
					DIRS SuperDeform
					INCDIRS SuperDeform)