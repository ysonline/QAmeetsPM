#define HAVE_SINGLEPRECISION_MATH //Eigen package does #include <complex> which conflicts with Inventor/C/basic.h definitions; so prevent them with this in files using include Inventor
#define HAVE_INT8_T //similarly Eigen & Inventor/Inventor/system/inttypes.h has a redefinition problem: int8_t; so prevent them with this in files using include Inventor

#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoIndexedLineSet.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoCoordinate3.h>

#include "Mesh.h"

class Painter
{
public:
	int count = 0;
	SoSeparator* shapeSep(Mesh* mesh, double deltaX = 0.0);
	SoSeparator* spheresSep(Mesh* mesh, double deltaX = 0.0);
	SoSeparator* spheresSep(Mesh* mesh1, Mesh* mesh2, double deltaX);
	SoSeparator* spheresSep2(Mesh* mesh, double deltaX = 0.0);
	SoSeparator* spheresSep3(Mesh* mesh, bool pathVert, double deltaX = 0.0);
	SoSeparator* spheresSep4(Mesh* mesh1, Mesh* mesh2, int worst, int best, double deltaX);
	SoSeparator* sphereSep(Mesh* mesh, int v, double deltaX = 0.0);
	SoSeparator* matchingLinesSep(Mesh* mesh1, Mesh* mesh2, double deltaX = 0.0);	
};
