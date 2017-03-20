// Include files for R3 package
#ifndef R3_INCLUDED
#define R3_INCLUDED


// Windows include files 

#ifdef _WIN32
#include <windows.h>
#pragma warning(disable:4996)
#endif

#ifndef M_PI
#define M_PI 3.1415925
#endif



// Dependency include files 

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif
#include "../R2/R2.h"



// Constant declarations

#define R3_X 0
#define R3_Y 1
#define R3_Z 2



// Class declarations

class R3Point;
class R3Vector;
class R3Line;
class R3Ray;
class R3Segment;
class R3Plane;
class R3Box;
class R3Sphere;
class R3Matrix;



// Class include files

#include "R3Point.h"
#include "R3Vector.h"
#include "R3Line.h"
#include "R3Ray.h"
#include "R3Segment.h"
#include "R3Plane.h"
#include "R3Matrix.h"
#include "R3Box.h"
#include "R3Sphere.h"
#include "R3Cylinder.h"
#include "R3Cone.h"
#include "R3Mesh.h"


// Utility include files

#include "R3Distance.h"



#endif
