
//#include "R2/R2.h"
//#include "R3/R3.h"
#include "SurfaceTextureMap.h"
/*#include "R3Scene.h"
//#include "raytrace.h"

	R2Pixel surfaceTextureMap(R2Image *image, R3Sphere *sphere, const R3Point &pos) {
		
		//obtain texture map coordinates from surface coordinates
		R3Vector V0 = pos - sphere->Center();
		double r = sphere->Radius();
		
		R3Vector V0xy = V0;
		V0xy.Project(R3xy_plane);
		double cphi = V0.Dot(V0xy) / r / V0xy.Length();
		R3Vector V0xy_xz = V0xy;
		V0xy_xz.Project(R3xz_plane);
		double ctheta = V0xy.Dot(V0xy_xz) / V0xy.Length() / V0xy_xz.Length();

		
		double phi = acos(cphi); 
		double theta = acos(ctheta);
		if (V0.Z()<0) phi *= -1;
		if (V0.Y()>0) theta *= -1;
		phi += M_PI;
		theta += M_PI;
		phi /= 2.0*M_PI;
		theta /= 2.0*M_PI;

		//printf("%g--%g--",phi,theta);
		//obtain pixel given texture coordinates with range 0,1
		return image->Sample(phi,theta,1);

	}

	R2Pixel GetTexturePixel(R2Image *image, const R3Node *node, const R3Point &pos) {
		switch(node->shape->type) {
		case R3_BOX_SHAPE:
		case R3_SPHERE_SHAPE:
			return surfaceTextureMap(image,node->shape->sphere,pos);
		case R3_CYLINDER_SHAPE:
				break;
		case R3_MESH_SHAPE:
				break;
		case R3_CONE_SHAPE:
				break;
		default:
			assert(false);
			return R2black_pixel;
	}
		
	}
*/