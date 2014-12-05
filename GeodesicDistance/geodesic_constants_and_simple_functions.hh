//Copyright (C) 2008 Danil Kirsanov, MIT License
#ifndef GEODESIC_CONSTANTS_20071231
#define GEODESIC_CONSTANTS_20071231

// some constants and simple math functions

#include <assert.h>
#include <math.h>
#include <limits>
#include <fstream>

namespace geodesic{

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//double const GEODESIC_INF = std::numeric_limits<double>::max();
double const GEODESIC_INF = 1e100;

//in order to avoid numerical problems with "infinitely small" intervals,
//we drop all the intervals smaller than SMALLEST_INTERVAL_RATIO*edge_length
double const SMALLEST_INTERVAL_RATIO = 1e-6;		
//double const SMALL_EPSILON = 1e-10;


inline double cos_from_edges(double const a,			//compute the cosine of the angle given the lengths of the edges
							 double const b,
							 double const c)
{
	assert(a>1e-50);
	assert(b>1e-50);
	assert(c>1e-50);

	double result = (b*b + c*c - a*a)/(2.0*b*c);
	result = std::max(result, -1.0);
	return std::min(result, 1.0);
}

inline double angle_from_edges(double const a,			//compute the cosine of the angle given the lengths of the edges
							   double const b,
							   double const c)
{
	return acos(cos_from_edges(a,b,c));
}

inline int get_first_integer( const char *v ){
    int ival;
    std::string s( v );
    std::replace( s.begin(), s.end(), '/', ' ' );
    sscanf( s.c_str(), "%d", &ival );
    return ival;
}

template<class Points, class Faces>
inline bool read_mesh_from_file(QString filename,
								Points& points,
								Faces& faces)
{
	std::ifstream file(filename.toUtf8().constData());
	assert(file.is_open());
	if(!file.is_open()) return false;
	
	unsigned num_points = 0;

    // read lines from the file, if the first character of the
    // line is 'v', we are reading a vertex, otherwise, if the
    // first character is a 'f' we are reading a facet
    char line[1024], v0[1024], v1[1024], v2[1024];
    while( file.getline(line, 1024) )
    {
        if( line[0] == 'v' && line[1] == ' '){
            double x, y, z;
            sscanf( line, "%*s%lf%lf%lf", &x, &y, &z );
            points.push_back(x);
            points.push_back(y);
            points.push_back(z);
            num_points += 3;
        } else if( line[0] == 'f' ){
            sscanf( line, "%*s%s%s%s", v0, v1, v2 );
            faces.push_back( get_first_integer( v0 )-1 );
            faces.push_back( get_first_integer( v1 )-1 );
            faces.push_back( get_first_integer( v2 )-1 );
        }
    }
	file.close();

    assert(num_points>=3);
    
	return true;
}

} //geodesic

#endif	//GEODESIC_CONSTANTS_20071231
