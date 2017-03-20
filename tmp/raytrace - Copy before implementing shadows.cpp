// Source file for raytracing code


// Include files

#ifdef _WIN32
#include <windows.h>
#endif

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "raytrace.h"

R3Point O;
R3Vector bottomframe, leftframe;

void GenEndpoints(R3Point eyept, double d, R3Vector up, R3Vector towards, R3Vector right, double xfov,double yfov) {

	assert(up.IsNormalized());
	assert(towards.IsNormalized());

	R3Point P1v, P2v, P1h, P2h;

	P1v = eyept + d*towards - d*tan(yfov)*up;
	P2v = eyept + d*towards + d*tan(yfov)*up;

	P1h = eyept + d*towards - d*tan(xfov)*right;
	P2h = eyept + d*towards + d*tan(xfov)*right;

	//move origin to lower left
	P1v -= d*tan(xfov)*right;
	P2v -= d*tan(xfov)*right;
	P1h -= d*tan(yfov)*up;
	P2h -= d*tan(yfov)*up;
	O = P1h; //should assert that P1h and P1v are apprx equal
	leftframe = P2v - P1v;
	bottomframe = P2h - P1h;

}

R3Ray GenRay(int i,int j, R3Point &eyept, R3Vector &up, R3Vector &down, double xfov, double yfov, 
			 int width, int height) {

				 assert(up.IsNormalized()); 
				 assert(down.IsNormalized());

				 R3Vector dx = ((double(i)-.5)/width)*bottomframe;
				 R3Vector dy = ((double(j)-.5)/height)*leftframe;

				 R3Ray ray(eyept,O+dx + dy);

				 return(ray);

}


bool R3Intersects(const R3Sphere &sphere, const R3Ray &ray, R3Point &pos, double &t, R3Vector &normal) {
	pos = R3Point(0,0,0);
	t=0;
	normal = R3Vector(0,0,0);
/*
	R3Vector L = (sphere.Center() - ray.Point(0));
	L += sphere.Center().Vector();

	double t_ca = L.Dot(ray.Vector());

	double d = sqrt(L.Dot(L) - t_ca*t_ca);
	if (d>sphere.Radius()) return false;

	double t_hc = sqrt(sphere.Radius()*sphere.Radius() - d*d);
	double t1 = t_ca - t_hc;
	double t2 = t_ca + t_hc;
*/
	R3Point P0(ray.Point(0));
	R3Point O(sphere.Center());
	R3Vector V(ray.Vector());
	V.Normalize();
	double r = sphere.Radius();

	double a = 1.0;
	double b = 2.0*V.Dot(P0 - O);
	double c = (P0 - O).Dot(P0 - O) - r*r;
	double discr = b*b - 4.0*a*c;

	if (discr<0) return false;
	double t1 = (-b + sqrt(discr))/(2.0*a);
	double t2 = (-b - sqrt(discr))/(2.0*a);

	if (abs(t1) < abs(t2)) t = t1; else t = t2;
	pos = ray.Point(t);
	normal = pos-sphere.Center();
	normal.Normalize();
	return true;
}


bool R3Intersects(const R3Plane plane, const R3Ray &ray, R3Point &pos, double &t, R3Vector &normal) {
	pos = R3Point(0,0,0);
	t=0;
	normal = plane.Normal();
	normal.Normalize();

	if (ray.Vector().Dot(normal)==0) return false;

	t = -(ray.Point(0).Vector().Dot(normal) + plane.D())/ray.Vector().Dot(normal);;

	pos = ray.Point(t);

	return true;
}

//whether the ray intersects the axes-parallel box on face specified by dim1
//helper for box version of r3intersects
bool IntersectsFace(R3Ray ray, R3Box box, int dim1, int min, double &t) {

	if (ray.Vector()[dim1]==0) return(false);

	int dim2,dim3;
	dim2 = 1-dim1;
	if (dim2<0) {dim2=0; dim3=1;} else dim3=2;

	double facecoord = box[min][dim1];

	t = (facecoord - ray.Point(0)[dim1])/ray.Vector()[dim1];

	return ((ray.Point(t)[dim2] < box[1][dim2] && ray.Point(t)[dim2] > box[0][dim2])&&
		(ray.Point(t)[dim3] < box[1][dim3] && ray.Point(t)[dim3] > box[0][dim3]));

}


bool R3Intersects(const R3Box box, const R3Ray &ray, R3Point &pos, double &t, R3Vector &normal) {


	pos = R3Point(0,0,0);
	t=MAXDWORD;
	normal = R3Vector(0,0,0);

	vector<double> ts(2);
	R3Point P0(ray.Point(0));
	R3Vector V(ray.Vector());

	double tx_min = V.X() != 0 ? (box.XMin() - P0.X()) / V.X() : MAXDWORD;
	double tx_max = V.X() != 0 ? (box.XMax() - P0.X()) / V.X() : MAXDWORD;
	double ty_min = V.Y() != 0 ? (box.YMin() - P0.Y()) / V.Y() : MAXDWORD;
	double ty_max = V.Y() != 0 ? (box.YMax() - P0.Y()) / V.Y() : MAXDWORD;
	double tz_min = V.Z() != 0 ? (box.ZMin() - P0.Z()) / V.Z() : MAXDWORD;
	double tz_max = V.Z() != 0 ? (box.ZMax() - P0.Z()) / V.Z() : MAXDWORD;

#define between(x,a,b) (((a<=x) && (x<=b)) || ((b<=x) && (x<=a)))

	if (tx_min>0 && between(ray.Point(tx_min).Y(),box.YMin(),box.YMax()) && between(ray.Point(tx_min).Z(),box.ZMin(),box.ZMax()))
		if (abs(tx_min) < abs(t)) {
			t = tx_min;
			normal = R3Vector(box.XMin() - box.XMax(),0,0);
		}
	if (tx_max>0 && between(ray.Point(tx_max).Y(),box.YMin(),box.YMax()) && between(ray.Point(tx_max).Z(),box.ZMin(),box.ZMax()))
		if (abs(tx_max) < abs(t)) {
			t = tx_max;
			normal = R3Vector(box.XMax() - box.XMin(),0,0);
		}
	if (ty_min>0 && between(ray.Point(ty_min).X(),box.XMin(),box.XMax()) && between(ray.Point(ty_min).Z(),box.ZMin(),box.ZMax()))
		if (abs(ty_min) < abs(t)) {
			t = ty_min;
			normal = R3Vector(0,box.YMin() - box.YMax(),0);
		}
	if (ty_max>0 && between(ray.Point(ty_max).X(),box.XMin(),box.XMax()) && between(ray.Point(ty_max).Z(),box.ZMin(),box.ZMax()))
		if (abs(ty_max) < abs(t)) {
			t = ty_max;
			normal = R3Vector(0,box.YMax() - box.YMin(),0);
		}
	if (tz_min>0 && between(ray.Point(tz_min).X(),box.XMin(),box.XMax()) && between(ray.Point(tz_min).Y(),box.YMin(),box.YMax()))
		if (abs(tz_min) < abs(t)) {
			t = tz_min;
			normal = R3Vector(0,0,box.ZMin() - box.ZMax());
		}
	if (tz_max>0 && between(ray.Point(tz_max).X(),box.XMin(),box.XMax()) && between(ray.Point(tz_max).Y(),box.YMin(),box.YMax()))
		if (abs(tz_max) < abs(t)) {
			t = tz_max;
			normal = R3Vector(0,0,box.ZMax() - box.ZMin());
		}

#undef between

/*	int d;

	double min_t = MAXDWORD;

	for (int dim1=0; dim1<3; dim1++)
		for (int max=0; max<=1; max++) {
			if (IntersectsFace(ray,box,dim1,max,t)) {
				if (abs(min_t) > abs(t)) {
					min_t = t;
					normal = R3Vector(0,0,0);
					normal[dim1] = -1*(1-2*max);
				}
			}
		}


		//	if (abs(t1) < abs(t2)) t = t1; else t = t2;
		if (min_t == MAXDWORD) return false;
		t = min_t;
*/
		if (t == MAXDWORD) return false;
		pos = ray.Point(t);
		normal.Normalize();
		return true;
		//	normal = pos-sphere.Center();
		//normal.Normalize();
}

bool R3Intersects(const R3MeshFace &meshFace, const R3Ray &ray, R3Point &pos, double &t, R3Vector &normal) {

	pos = R3Point(0,0,0);
	t=0;
	normal = R3Vector(0,0,0);

	int d;

	if (!R3Intersects(meshFace.plane,ray,pos,t,normal)) return false;

	R3Point P0 = ray.Point(0);
	R3Vector v1,v2,N;
	double res=0;
	for (int side_i = 0; side_i < 3; side_i++) {
		v1 = meshFace.vertices[side_i]->position - P0;
		v2 = meshFace.vertices[(side_i+1)%3]->position - P0;
		N = v1;
		N.Cross(v2);

		if (res==0 || (N.Dot(pos-P0)*res) > 0) res = N.Dot(pos-P0);
		else return false;
	}

	return true;


	//	if (abs(t1) < abs(t2)) t = t1; else t = t2;
	//		t = (min_t == MAXDWORD) ? 0 : min_t;
	//	pos = ray.Point(t);
	//	normal = pos-sphere.Center();
	//normal.Normalize();
}

bool R3Intersects(const R3Mesh &mesh, const R3Ray &ray, R3Point &pos, double &t, R3Vector &normal) {

	pos = R3Point(0,0,0);
	t=0;
	normal = R3Vector(0,0,0);

	t = MAXDWORD;
	double t_i;
	R3Point pos_i;
	R3Vector normal_i;

	for (int i=0; i<mesh.NFaces(); i++) 
		if (R3Intersects(*mesh.Face(i),ray,pos_i,t_i,normal_i))
			if (abs(t_i) < abs(t)) {
				t = t_i;
				pos = pos_i;
				normal = normal_i;
			}

			if (t==MAXDWORD) return false;
			else return true;

}
/*
//project point onto a line
//helper used for R3Intersects(cylidner)
double projectPointLine(R3Point p0, R3Ray ray) {
R3Point p1 = ray.Point(0);
R3Vector v = ray.Vector();
v.Normalize();
double t = (p0-p1).Dot(v);
return abs(t);
}
*/
bool R3Intersects(const R3Cylinder &cylinder, const R3Ray &ray, R3Point &pos, double &t, R3Vector &normal) {

	pos = R3Point(0,0,0);
	t=0;
	normal = R3Vector(0,0,0);

	double t1,t2,t3;

	R3Vector C = cylinder.Center().Vector();
	R3Plane plane(R3Point(0,0,0),cylinder.Axis().Vector());
	C.Project(plane);
	double r = cylinder.Radius();
	R3Vector P = ray.Point(0).Vector();
	P.Project(plane);
	R3Vector V = ray.Vector();
	V.Project(plane);
	V.Normalize();

	R3Vector N = cylinder.Axis().Vector();
	R3Plane cap_plane1(cylinder.Axis().Start(),N);
	R3Plane cap_plane2(cylinder.Axis().End(),N);
	double d1 = cap_plane1.D();
	double d2 = cap_plane2.D();
	R3Vector N1 = cap_plane1.Normal();
	R3Vector N2 = cap_plane2.Normal();

	if (V.Length() == 0) //ray parallel to axis
		t1 = MAXDWORD;
	else {
		double discr = 4.0*pow(P.Dot(V) - C.Dot(V),2) - 4.0*(P.Dot(P) + C.Dot(C) - 2.0*P.Dot(C) - r*r);

		if (discr<0) t1 = MAXDWORD;
		//else t1 = min(abs((-2*(P.Dot(V)-C.Dot(V)) + sqrt(discr))/(2*V.Dot(V))),abs((-2*(P.Dot(V)-C.Dot(V)) - sqrt(discr))/(2*V.Dot(V))));
		else {
			double r1 = (-2.0*(P.Dot(V)-C.Dot(V)) + sqrt(discr))/2.0;
			double r2 = (-2.0*(P.Dot(V)-C.Dot(V)) - sqrt(discr))/2.0;
			t1 = (abs(r1)<abs(r2)) ? r2 : r2;
		}

		//is intersection with side outside segment of cylinder?
		R3Point proj1 = ray.Point(t1); R3Point proj2 = proj1;
		proj1.Project(cap_plane1); proj2.Project(cap_plane2);
		if (t1 != MAXDWORD)
			if (((ray.Point(t1) - proj1).Length() > cylinder.Height()) ||
				((ray.Point(t1) - proj2).Length() > cylinder.Height()))
				t1 = MAXDWORD;
	}


	if (N.Dot(ray.Vector())==0) {//if ray orthogonal to cylinder axis 
		if (t1==MAXDWORD) return false; //might as well return immediately if no intersection possible
	}
	else {

		//check for intersection with caps
		//if ray vector parallel to cylidner vector, no intersection
		//(ignores case when camera eye is inside cylinder)

		//get intersections with plances of the caps
		V = ray.Vector(); V.Normalize();
		P = ray.Point(0).Vector();
		t2 = (d1 - N1.Dot(P))/(N1.Dot(V));
		t3 = (d2 - N2.Dot(P))/(N2.Dot(V));

		//is intersection with cap plane outside the cap circle?
		if ((ray.Point(t2) - cylinder.Axis().Start()).Length() > r) t2 = MAXDWORD;
		if ((ray.Point(t3) - cylinder.Axis().End()).Length() > r) t3 = MAXDWORD;
	}

	if (t1==MAXDWORD && t2==MAXDWORD && t3==MAXDWORD) return false;

	if (abs(t1) < abs(t2) && abs(t1) < abs(t3)) {
		t = t1;
		R3Point proj = ray.Point(t1);           
		proj.Project(cylinder.Axis().Line());
		normal = ray.Point(t1)-proj;
		normal.Normalize();
		return true;
	}
	normal = N;
	if (abs(t3)<abs(t2)) {
		t = t3;
		normal *= -1.0;
	} 
	else t = t3;
	return true;

}

R3Vector transformCol(R3Vector A0, R3Vector A1, R3Vector A2, R3Vector X) {
	return R3Vector(X.Dot(A0),X.Dot(A1),X.Dot(A2));
}

bool R3Intersects(const R3Cone cone, const R3Ray &ray, R3Point &pos, double &t, R3Vector &normal) {

	pos = R3Point(0,0,0);
	t=0;
	normal = R3Vector(0,0,0);
	//printf("cone intersection not implemented\n");

	R3Vector P = ray.Point(0).Vector();
	R3Vector D = ray.Vector();
	D.Normalize();
	R3Vector V = cone.Axis().Start().Vector();
	R3Vector A = cone.Axis().Vector();
	A.Normalize();
	double theta = atan2(cone.Radius(),cone.Height());
	double gamma = cos(theta);

	double a0=A.X(),a1=A.Y(),a2=A.Z();
	R3Vector M0 = a0*R3Vector(a0,a1,a2);
	R3Vector M1 = a1*R3Vector(a0,a1,a2);
	R3Vector M2 = a2*R3Vector(a0,a1,a2);
	M0 -= R3Vector(gamma*gamma,0,0);
	M1 -= R3Vector(0,gamma*gamma,0);
	M2 -= R3Vector(0,0,gamma*gamma);

	//double M = A.Dot(A) - pow(cos(theta),2);
	R3Vector Delta = P - V;
	double c2 = D.Dot(transformCol(M0,M1,M2,D));
	double c1 = D.Dot(transformCol(M0,M1,M2,Delta));
	double c0 = Delta.Dot(transformCol(M0,M1,M2,Delta));

	double delta = c1*c1-c0*c2;
	if (delta < 0) return false;
	t = min(abs((-c1+sqrt(delta))/c2),abs((-c1-sqrt(delta))/c2));
	if (A.Dot(ray.Point(t).Vector()-V) < 0) return false;
	R3Point proj = ray.Point(t);
	proj.Project(cone.Axis().Line());
	if ((proj-cone.Axis().Start()).Length() > cone.Height()) return false;
	pos = ray.Point(t);
	//normal is normal to normal section through axis and line containing pos and head, and normal to pos itself
	//R3Line surface_line(cone.Axis().Start(),pos);
	R3Vector n1(cone.Axis().Start()-pos);
	R3Vector n2(n1);
	n2.Cross(cone.Axis().Vector());
	normal = n1;
	normal.Cross(n2);
	normal.Normalize();
	return true;

}

void traverseTree(R3Node *node,const R3Ray &ray, double &t, R3Vector &normal, R3Node *&min_node) {

	double min_t = MAXDWORD;
	R3Vector min_normal(0,0,0);
	min_node = 0;
	//double t;

	//check if this node's shapes is intersected by ray
	//if so, set t to the minimum of interseting shapes
	if (node->shape)
		if (node->shape->type==R3_BOX_SHAPE) {
			if (R3Intersects(*(node->shape->box),ray,R3Point(),t,normal))
				if (abs(t)<abs(min_t)) {
					min_t = t;
					min_normal = normal;
					min_node = node;
				}
		} else
			if (node->shape->type==R3_SPHERE_SHAPE) {
				if (R3Intersects(*(node->shape->sphere),ray,R3Point(),t,normal))
					if (abs(t)<abs(min_t)) {
						min_t = t;
						min_normal = normal;
						min_node = node;
					}
			} else
				if (node->shape->type==R3_CYLINDER_SHAPE) {
					if (R3Intersects(*(node->shape->cylinder),ray,R3Point(),t,normal))
						if (abs(t)<abs(min_t)) {
							min_t = t;
							min_normal = normal;
							min_node = node;
						}
				} else
					if (node->shape->type==R3_MESH_SHAPE) {
						if (R3Intersects(*(node->shape->mesh),ray,R3Point(),t,normal))
							if (abs(t)<abs(min_t)) {
								min_t = t;
								min_normal = normal;
								min_node = node;
							}
					} else
						if (node->shape->type==R3_CONE_SHAPE) {
							if (R3Intersects(*(node->shape->cone),ray,R3Point(),t,normal))
								if (abs(t)<abs(min_t)) {
									min_t = t;
									min_normal = normal;
									min_node = node;
								}
						}


						//for each child of this node, check if child is intersected by this node
						//if so, set t to minimum of (children's point of intersection, t)
						for (int i=0; i<node->children.size(); i++) {
							R3Node *min_node_tmp;
							traverseTree(node->children[i],ray,t,normal,min_node_tmp);
							if (abs(t)<abs(min_t)) {
								min_t = t;
								min_normal = normal;
								min_node = min_node_tmp;
							}
						}

						t = min_t;
						normal = min_normal;
						//node = min_node;
}

bool R3Intersects(const R3Scene *scene, const R3Ray &ray, R3Node *&node, R3Point &pos, double &t, R3Vector &normal) {

	pos = R3Point(0,0,0);
	t=MAXDWORD;
	normal = R3Vector(0,0,0);
	node = 0;

	//	R3Mesh *s = scene->root->children[0]->shape->mesh;

	//R3Node *node_i = scene->root;
	traverseTree(scene->root,ray,t,normal,node);

	if (t==MAXDWORD) return false;
	else {
		pos = ray.Point(t);
		return true;
	}

}

//helper for getting specular reflection; reflects vector through another vector, in their common plane
R3Vector vecVecReflect(const R3Vector &vec, const R3Vector &axis) {
	assert(axis.IsNormalized());
	double t = vec.Dot(axis);
	return vec - 2*(vec-t*axis);
}

R2Pixel GetColor(R3Scene *scene, R3Ray &ray, R3Node *node, R3Point &pos, const double &t, const R3Vector &normal) {

	vector<R2Pixel> I_L(scene->NLights());
	vector<boolean> S_L(scene->NLights());
	R3Vector L,R,V;

	//calculate intensity of each light source at given point of object
	for (int i=0; i<scene->NLights(); i++) {
		double d,k_c,k_l,k_q,alpha,gamma,ctheta;
		R3Vector L;
		R3Light *light(scene->Light(i));
		R2Pixel I_0 = light->color;
		switch (light->type) {
			case R3_DIRECTIONAL_LIGHT:
				I_L[i] = I_0;
				break;
			case R3_POINT_LIGHT:
				d = (pos-light->position).Length();
				k_c = light->constant_attenuation;
				k_l = light->linear_attenuation;
				k_q = light->quadratic_attenuation;
				I_L[i] = I_0 / (k_c + k_l*d + k_q*d*d);
				break;
			case R3_SPOT_LIGHT:
				d = (pos-light->position).Length();
				k_c = light->constant_attenuation;
				k_l = light->linear_attenuation;
				k_q = light->quadratic_attenuation;
				alpha = light->angle_attenuation;
				gamma = light->angle_cutoff;
				L = pos-light->position;
				L.Normalize();
				assert(normal.IsNormalized());
				ctheta = normal.Dot(L);
				if (ctheta < cos(gamma))
				I_L[i] = I_0 * pow(ctheta,alpha) / (k_c + k_l*d + k_q*d*d);
				else
					I_L[i] = 0;
				break;
			default:
				assert(false);
		}
	}

	//obtain shadow switch for each light source
	for (int i=0; i<scene->NLights(); i++) {

	}


	//material properties
	R3Material *material = node->material;
	int n = material->shininess;
	R2Pixel K_A=material->ka;
	R2Pixel K_D=material->kd;
	R2Pixel K_S=material->ks;

	R2Pixel I_D(R2black_pixel);
	R2Pixel I_S(R2black_pixel);
	R2Pixel I_A = K_A*scene->ambient;
	R2Pixel I_E(material->emission);

	//calculate diffuse/specular reflection from each light source
	for (int i=0; i<scene->NLights(); i++) {
		L = -scene->Light(i)->direction;
		L.Normalize();
		//obtain I_D (diffuse reflection)
		if (L.Dot(normal)>0) {
		assert(normal.IsNormalized());
		I_D += K_D * (normal.Dot(L)) * I_L[i];
		}
		//R3Vector ntmp(pos-R3Point(-.5,0,0)); //ntmp.Normalize();
		//I_D += K_D * (-R3Vector(-1,0,0).Dot(R3Vector(0,normal.Y(),0))) * R2Pixel(1,1,1,0);
		//obtain I_S (specular)	
		if (L.Dot(normal)>0) {
		R = vecVecReflect(L,normal); //mirror direction
		V = -ray.Vector();
		assert(V.IsNormalized());
		I_S += K_S * pow(V.Dot(R),double(n)) * I_L[i];
		}
	}
	//printf("%g:%g\n",L.Dot(normal),I_D.Red());
	//I_D /= scene->NLights();
	//I_S /= scene->NLights();
	R2Pixel I = I_A+I_E+I_D+I_S; 
	I.Clamp();
	return I;

}	


// Create image from scene
// This is the main ray tracing function 

R2Image *RenderImage(R3Scene *scene, int width, int height, int max_reflections, double min_luminance)
{
	// Allocate  image
	R2Image *image = new R2Image(width, height);
	if (!image) {
		fprintf(stderr, "Unable to allocate image\n");
		return NULL;
	}
	//scene->lights.
	//SCENE PARAMETERS
	/*	R3Point eyept(0,0,5);
	double d=2;
	R3Vector up(0,1,0);
	double xfov=.5;
	double yfov=.5;
	*/

	R3Point eyept = scene->camera.eye;
	double d = scene->camera.fardist;
	double xfov = scene->camera.xfov;
	double yfov = scene->camera.yfov;
	R3Vector towards = scene->camera.towards;
	R3Vector up = scene->camera.up;
	R3Vector right = scene->camera.right;


	//up.Normalize();
	/*	R3Vector towards = R3Point(0,0,0)-eyept;
	towards.Normalize();
	*/
	double t;
	R3Point pos;
	R3Vector normal;
	R3Node *node;


	//R3Box box(R3Point(-1,-1,-1),R3Point(1,1,1));
	//R3Cylinder cylinder(R3Point(0,0,0),1,1);

	//generage frame for image plane
	GenEndpoints(eyept,d,up,towards,right,xfov,yfov);

	R3Ray r(eyept,R3Vector(0,0,1));
	R3Intersects(scene,r,node,pos,t,normal);


	for (int i=0; i<image->Width(); i++) {
		for (int j=0; j<image->Height(); j++) {
			//printf("%d, %d: ",i,j);
			R3Ray ray = GenRay(i,j,eyept,up,towards,xfov,yfov,width,height);
			if (R3Intersects(scene,ray,node,pos,t,normal)) {			
				image->Pixel(i,j) = GetColor(scene,ray,node,pos,t,normal);
//				pos.Print();printf("--");normal.Print();printf("--");(pos-R3Vector(-.5,0,0)).Print();printf("\n");
			}
			else image->Pixel(i,j) = R2black_pixel;
			//if (((i==image->Width()-1) || i==0) && (j==0 || j==image->Height()-1)) {
			//if (i==255/2 && j==255/2) {
			//	printf("%d, %d: ",i,j); ray.Print();printf("\n");}
		}
		//printf("\n");
	}

	// FILL IN YOUR CODE HERE
	//fprintf(stderr, "Not implemented yet\n"); return 0;

	// Return image
	return image;
}
