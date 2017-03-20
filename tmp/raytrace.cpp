// Source file for raytracing code


// Include files

#ifdef _WIN32
#include <windows.h>
#endif

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "raytrace.h"



bool ANTI_ALIASING=true;

#define rand_unit() (double(rand())/(RAND_MAX+1))
#define rand_int(max_int) (floor((rand_unit())*(max_int))+1)
#define rand_interp(a,b) (min(a,b)+rand_unit()*abs(b-a))

R3Point O;
R3Vector bottomframe, leftframe;

void GenEndpoints(R3Point eye, double d, R3Vector up, R3Vector towards, R3Vector right, double xfov,double yfov) {

	assert(up.IsNormalized());
	assert(towards.IsNormalized());

	R3Point P1v, P2v, P1h, P2h;

	P1v = eye + d*towards - d*tan(yfov)*up;
	P2v = eye + d*towards + d*tan(yfov)*up;

	P1h = eye + d*towards - d*tan(xfov)*right;
	P2h = eye + d*towards + d*tan(xfov)*right;

	//move origin to lower left
	P1v -= d*tan(xfov)*right;
	P2v -= d*tan(xfov)*right;
	P1h -= d*tan(yfov)*up;
	P2h -= d*tan(yfov)*up;
	O = P1h; //should assert that P1h and P1v are apprx equal
	leftframe = P2v - P1v;
	bottomframe = P2h - P1h;

}

R3Ray  GenRay(int i,int j, int sample_no, int pixels_per_side,
			  R3Point &eye, R3Vector &up, R3Vector &down, double xfov, double yfov, 
			  int width, int height) 
{

	//assert(up.IsNormalized()); 
	//assert(down.IsNormalized());

	//double x_offset = ANTI_ALIASING ? (double(rand())/(RAND_MAX+1)) - 0.5 : 0.5;
	//double y_offset = ANTI_ALIASING ? (double(rand())/(RAND_MAX+1)) - 0.5 : 0.5;

	double stratum_side = 1.f/pixels_per_side;
	double x_offset_min = floor(double(sample_no) / pixels_per_side) * stratum_side;
	double y_offset_min = (sample_no % pixels_per_side) * stratum_side;
	double x_offset = x_offset_min + rand_unit()*stratum_side;
	x_offset -= .5;
	double y_offset = y_offset_min + rand_unit()*stratum_side;
	y_offset -= .5;

	R3Vector dx = ((double(i)-x_offset)/width)*bottomframe;
	R3Vector dy = ((double(j)-y_offset)/height)*leftframe;

	R3Ray ray(eye,O+dx + dy);

	return(ray);

}


bool R3Intersects(const R3Sphere &sphere, const R3Ray &ray, R3Point &pos, double &t, R3Vector &normal) {
	pos = R3Point(0,0,0);
	t=-1;
	normal = R3Vector(0,0,0);

//	if ((ray.Point(0) - sphere.Center()).Length() < sphere.Radius()) return false;
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

	if (t1*t2<0) t = max(t1,t2);
	else
		if (t1<0) return false; //both negative
		else t = min(t1,t2);
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

	if (t<0) return false;

	pos = ray.Point(t);

	return true;
}

/*
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
*/

bool R3Intersects(const R3Box box, const R3Ray &ray, R3Point &pos, double &t, R3Vector &normal) {


	pos = R3Point(0,0,0);
	t=-1;
	normal = R3Vector(0,0,0);

	vector<double> ts(2);
	R3Point P0(ray.Point(0));
	R3Vector V(ray.Vector());

	double tx_min = V.X() != 0 ? (box.XMin() - P0.X()) / V.X() : -1;
	double tx_max = V.X() != 0 ? (box.XMax() - P0.X()) / V.X() : -1;
	double ty_min = V.Y() != 0 ? (box.YMin() - P0.Y()) / V.Y() : -1;
	double ty_max = V.Y() != 0 ? (box.YMax() - P0.Y()) / V.Y() : -1;
	double tz_min = V.Z() != 0 ? (box.ZMin() - P0.Z()) / V.Z() : -1;
	double tz_max = V.Z() != 0 ? (box.ZMax() - P0.Z()) / V.Z() : -1;

#define between(x,a,b) (((a<=x) && (x<=b)) || ((b<=x) && (x<=a)))

	if (tx_min>0 && between(ray.Point(tx_min).Y(),box.YMin(),box.YMax()) && between(ray.Point(tx_min).Z(),box.ZMin(),box.ZMax()))
		if (tx_min<t || t<0) {
			t = tx_min;
			normal = R3Vector(box.XMin() - box.XMax(),0,0);
		}
		if (tx_max>0 && between(ray.Point(tx_max).Y(),box.YMin(),box.YMax()) && between(ray.Point(tx_max).Z(),box.ZMin(),box.ZMax()))
			if (tx_max<t || t<0) {
				t = tx_max;
				normal = R3Vector(box.XMax() - box.XMin(),0,0);
			}
			if (ty_min>0 && between(ray.Point(ty_min).X(),box.XMin(),box.XMax()) && between(ray.Point(ty_min).Z(),box.ZMin(),box.ZMax()))
				if (ty_min<t || t<0) {
					t = ty_min;
					normal = R3Vector(0,box.YMin() - box.YMax(),0);
				}
				if (ty_max>0 && between(ray.Point(ty_max).X(),box.XMin(),box.XMax()) && between(ray.Point(ty_max).Z(),box.ZMin(),box.ZMax()))
					if (ty_max<t || t<0) {
						t = ty_max;
						normal = R3Vector(0,box.YMax() - box.YMin(),0);
					}
					if (tz_min>0 && between(ray.Point(tz_min).X(),box.XMin(),box.XMax()) && between(ray.Point(tz_min).Y(),box.YMin(),box.YMax()))
						if (tz_min<t || t<0) {
							t = tz_min;
							normal = R3Vector(0,0,box.ZMin() - box.ZMax());
						}
						if (tz_max>0 && between(ray.Point(tz_max).X(),box.XMin(),box.XMax()) && between(ray.Point(tz_max).Y(),box.YMin(),box.YMax()))
							if ((tz_max<t || t<0)) {
								t = tz_max;
								normal = R3Vector(0,0,box.ZMax() - box.ZMin());
							}

#undef between


							if (t < 0) return false;
							pos = ray.Point(t);
							normal.Normalize();
							return true;
							//	normal = pos-sphere.Center();
							//normal.Normalize();
}

bool R3Intersects(const R3MeshFace &meshFace, const R3Ray &ray, R3Point &pos, double &t, R3Vector &normal) {

	pos = R3Point(0,0,0);
	t=-1;
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
	normal = R3Vector(0,0,0);

	t = -1;
	double t_i;
	R3Point pos_i;
	R3Vector normal_i;

	for (int i=0; i<mesh.NFaces(); i++) 
		if (R3Intersects(*mesh.Face(i),ray,pos_i,t_i,normal_i))
			if (t_i>0 && (t_i<t || t<0)) {
				t = t_i;
				pos = pos_i;
				normal = normal_i;
			}

			if (t<0) return false;
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
	t=-1;
	normal = R3Vector(0,0,0);

	double t1=-1,t2=-1;

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
	R3Plane cap_plane1(cylinder.Axis().Start(),-N);
	R3Plane cap_plane2(cylinder.Axis().End(),N);
	/*	double d1 = cap_plane1.D();
	double d2 = cap_plane2.D();
	R3Vector N1 = cap_plane1.Normal();
	R3Vector N2 = cap_plane2.Normal();
	*/
	if (V.Length() == 0) //ray parallel to axis
		return false;

	double discr = 4.0*pow(P.Dot(V) - C.Dot(V),2) - 4.0*(P.Dot(P) + C.Dot(C) - 2.0*P.Dot(C) - r*r);

	if (discr<0) 
		return false;
	//else t1 = min(abs((-2*(P.Dot(V)-C.Dot(V)) + sqrt(discr))/(2*V.Dot(V))),abs((-2*(P.Dot(V)-C.Dot(V)) - sqrt(discr))/(2*V.Dot(V))));

	t1 = (-2.0*(P.Dot(V)-C.Dot(V)) + sqrt(discr))/2.0;
	t2 = (-2.0*(P.Dot(V)-C.Dot(V)) - sqrt(discr))/2.0;
	/*if (r1*r2<0) t1 = max(r1,r2);
	else
	if (r1<0) t1 = -1; //both negative
	else t1 = min(r1,r2); //both positive
	//t1 = (abs(r1)<abs(r2)) ? r2 : r2;*/

	//is intersection with side outside segment of cylinder?
	if (t1>=0) {
	R3Point proj1 = ray.Point(t1); R3Point proj2 = proj1;
	proj1.Project(cap_plane1); proj2.Project(cap_plane2);
		if (((ray.Point(t1) - proj1).Length() > cylinder.Height()) ||
			((ray.Point(t1) - proj2).Length() > cylinder.Height()))
			t1 = -1;
	}
if (t2>=0) {
	R3Point proj1 = ray.Point(t2); R3Point proj2 = proj1;
	proj1.Project(cap_plane1); proj2.Project(cap_plane2);
		if (((ray.Point(t2) - proj1).Length() > cylinder.Height()) ||
			((ray.Point(t2) - proj2).Length() > cylinder.Height()))
			t2 = -1;
	}

	if (t1<0 && t2<0)
		return false;
	else if (t1*t2<0)
	t = max(t1,t2);
	else t = min(t1,t2);

	pos = ray.Point(t);
	R3Point proj(pos);
	proj.Project(plane);
	normal = pos-proj;
	normal.Normalize();
	return true;

/*
	if (t1>0 && (t1<t2 || t2<0) && (t1<t3 || t3<0)) {
		t = t1;
		R3Point proj = ray.Point(t1);           
		proj.Project(cylinder.Axis().Line());
		normal = ray.Point(t1)-proj;
		normal.Normalize();
		return true;
	}
	normal = N;
	if (t3>0 && (t3<t2 || t2<0)) {
		normal *= -1.0;
		t = t3;
	} 
	else {
		t = t2;
	}
	return true;
*/
}



//helper for cone implementation of r3intersects
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
	//return true;

	//R3Point ddd = cone.Axis().Start();
	/*if (normal == R3Vector(0,0,0)) { //if normal not defined, use axis
	normal = cone.Axis().Vector();
	return true;
	}*/
	//double tt=normal.X();double dd=normal.Y();

	return true;

}

void traverseTree(R3Node *node, const R3Ray &ray, R3Point &pos, double &t, R3Vector &normal, R3Node *&min_node) {

	double min_t = -1;//MAXDWORD;
	R3Vector min_normal(0,0,0);
	min_node = 0;
	//double t;

	//if non-identity transformation, transform ray by inverse
	R3Ray rayT(ray);
	if (!node->transformation.IsIdentity()) {
		rayT.InverseTransform(node->transformation);
	}

	//check if this node's shapes is intersected by ray
	//if so, set t to the minimum of interseting shapes
	if (node->shape)
		switch(node->shape->type) {
		case R3_BOX_SHAPE:
			if (R3Intersects(*(node->shape->box),rayT,R3Point(),t,normal))
				if (t>0 && (t<min_t || min_t<0)) {
					min_t = t;
					min_normal = normal;
					min_node = node;
				}
				break;
		case R3_SPHERE_SHAPE:
			if (R3Intersects(*(node->shape->sphere),rayT,R3Point(),t,normal))
				if (t>0 && (t<min_t || min_t<0)) {
					min_t = t;
					min_normal = normal;
					min_node = node;
				}
				break;
		case R3_CYLINDER_SHAPE:
			if (R3Intersects(*(node->shape->cylinder),rayT,R3Point(),t,normal))
				if (t>0 && (t<min_t || min_t<0)) {
					min_t = t;
					min_normal = normal;
					min_node = node;
				}
				break;
		case R3_MESH_SHAPE:
			if (R3Intersects(*(node->shape->mesh),rayT,R3Point(),t,normal))
				if (t>0 && (t<min_t || min_t<0)) {
					min_t = t;
					min_normal = normal;
					min_node = node;
				}
				break;
		case R3_CONE_SHAPE:
			if (R3Intersects(*(node->shape->cone),rayT,R3Point(),t,normal))
				if (t>0 && (t<min_t || min_t<0)) {
					min_t = t;
					min_normal = normal;
					min_node = node;
				}
				break;
		default:
			assert(false);
	}

	//unaccelerated intersection

	//for each child of this node, check if child is intersected by this node
	//if so, set t to minimum of (children's point of intersection, t)

	for (int i=0; i<node->children.size(); i++) {
		R3Node *min_node_tmp;//used to return the node found in the recursion
		if (/*R3Intersects(node->children[i]->bbox,rayT,pos,t,normal) && (t<min_t || min_t<0)*/true) {
		traverseTree(node->children[i],rayT,pos,t,normal,min_node_tmp);
		if (t>0 && (t<min_t || min_t<0)) {
			min_t = t;
			min_normal = normal;
			min_node = min_node_tmp;

		}
		}
	}


	t = min_t;
	pos = ray.Point(t);
	normal = min_normal;
	//cached_intersection = node;
	//node = min_node;

	//if transformation is not identity, transform hit information
	if ( t > 0 && !node->transformation.IsIdentity()) {
		pos.Transform(node->transformation);
		normal.InverseTransform(node->transformation.Transpose());
		normal.Normalize();
		t = min((pos.X() - ray.Point(0).X())/ray.Vector().X(),
			min((pos.Y() - ray.Point(0).Y())/ray.Vector().Y(),(pos.Z() - ray.Point(0).Z())/ray.Vector().Z()));
	}
}


bool R3Intersects(const R3Scene *scene, const R3Ray &ray, R3Node *&node, R3Point &pos, double &t, R3Vector &normal) {

	pos = R3Point(0,0,0);
	t = -1;
	normal = R3Vector(0,0,0);
	node = 0;


	traverseTree(scene->root,ray,pos,t,normal,node);

	if (t<0) {
		return false;
	}
	else {
		return true;
	}

}


//helper for getting specular reflection direction; reflects vector through another vector, in their common plane
R3Vector reflect_vector(const R3Vector &vec, const R3Vector &axis) {
	/*if (!axis.IsNormalized()) {
		double d=3;
	}*/
	assert(axis.IsNormalized());
	double t = vec.Dot(axis);
	return vec - 2*(vec-t*axis);
}

//helper for getting direction of refraction
R3Vector refract_vector(const R3Vector &L, const R3Vector normal, double ir) {
	ir = 1.0/ir;
	assert(L.IsNormalized());
	assert(normal.IsNormalized());
	double ctheta_i = normal.Dot(-L);
	double theta_i = acos(ctheta_i);
	double theta_r = asin(ir*sin(theta_i));
	return (ir*ctheta_i - cos(theta_r))*normal + ir*L;
}

/*
bool sample(R3Scene *scene, R3Light *light, const R3Point &pos, R3Vector &wi_L, double &t_L, double &L) {

	R3Node *node_tmp; R3Vector normal_tmp; R3Point pos_tmp(pos);
	wi_L = light->position - pos; 
	R3Ray wi_L_ray = R3Ray(pos,wi_L);
	R3Intersects(scene,wi_L_ray,node_tmp,pos_tmp,t_L,normal_tmp);
	if (wi_L.Length() < t_L)
		return false;
	else {
		wi_L.Normalize();
		return true;
	}

}
*/

R3Point sample_quadrilateral(R3Point v0, R3Point v1, R3Point v2) {
	//get a sample point on the quadrilateral using uniform sample along two axes
//			srand(time(NULL));

	R3Vector a = v2 - v0;
	R3Vector b = v1 - v0;
	double rand_a = rand_unit(), rand_b = rand_unit();
	return (rand_a*a + rand_b*b).Point();
}


class quadrilateral {
public:
	R3Point min,max;
	int fixed_dim,other_dim_1,other_dim_2;
	quadrilateral(R3Point min, R3Point max) {
		this->min = min; this->max = max;
		fixed_dim = 2; other_dim_1 = 0; other_dim_2 = 1;
		if (min.X()==max.X()) {
			fixed_dim = 0; other_dim_1 = 1; other_dim_2 = 2;
		} else if (min.Y()==max.Y()) {
			fixed_dim = 1; other_dim_1 = 0; other_dim_2 = 2;
		}
	}	

	double visibility(R3Point pos) {
		R3Point corner1;
		corner1[fixed_dim] = min[fixed_dim];
		corner1[other_dim_1] = min[other_dim_1];
		corner1[other_dim_2] = min[other_dim_2];

		R3Point corner2;
		corner2[fixed_dim] = min[fixed_dim];
		corner2[other_dim_1] = min[other_dim_1];
		corner2[other_dim_2] = max[other_dim_2];

		R3Point corner3;
		corner3[fixed_dim] = min[fixed_dim];
		corner3[other_dim_1] = max[other_dim_1];
		corner3[other_dim_2] = min[other_dim_2];
		
		R3Point corner4;
		corner4[fixed_dim] = min[fixed_dim];
		corner4[other_dim_1] = max[other_dim_1];
		corner4[other_dim_2] = max[other_dim_2];

		R3Vector e1(0,0,0),e2(0,0,0);
		e1[other_dim_1] = 1;
		e2[other_dim_2] = 1;
		R3Vector normal = e1; normal.Cross(e2);

		double a = abs(normal.Dot(corner1-pos)/(corner1-pos).Length());
		double b = abs(normal.Dot(corner2-pos)/(corner2-pos).Length());
		double c = abs(normal.Dot(corner3-pos)/(corner3-pos).Length());
		double d = abs(normal.Dot(corner4-pos)/(corner4-pos).Length());

		return (a+b+c+d)/4.f;

	}
	double area() {
		return abs(min[other_dim_1]-max[other_dim_1])*abs(min[other_dim_2]-max[other_dim_2]);
	}
	R3Point sample() {
		R3Point res;
		res[fixed_dim] = min[fixed_dim]; //==max[fixed_dim]
		res[other_dim_1] = rand_interp(min[other_dim_1],max[other_dim_1]);
		res[other_dim_2] = rand_interp(min[other_dim_2],max[other_dim_2]);
		return res;
	}
};

bool sample(const R3Box *box, const R3Point &pos, const R3Vector &normal,
			R3Point &pos2,  R3Vector &normal2, double &pdf) {

	double xmin = min(box->XMin(),box->XMax());
	double xmax = max(box->XMin(),box->XMax());
	double ymin = min(box->YMin(),box->YMax());
	double ymax = max(box->YMin(),box->YMax());
	double zmin = min(box->ZMin(),box->ZMax());
	double zmax = max(box->ZMin(),box->ZMax());

	double area = 0;

	vector<quadrilateral> box_sides;

	if (pos.X() <= xmin) 
		box_sides.push_back(quadrilateral(R3Point(xmin,ymin,zmin),R3Point(xmin,ymax,zmax)));
	if (pos.X() >= xmax) 
		box_sides.push_back(quadrilateral(R3Point(xmax,ymin,zmin),R3Point(xmax,ymax,zmax)));
	if (pos.Y() <= ymin) 
		box_sides.push_back(quadrilateral(R3Point(xmin,ymin,zmin),R3Point(xmax,ymin,zmax)));
	if (pos.Y() >= ymax)
		box_sides.push_back(quadrilateral(R3Point(xmin,ymax,zmin),R3Point(xmax,ymax,zmax)));
	if (pos.Z() <= zmin) 
		box_sides.push_back(quadrilateral(R3Point(xmin,ymin,zmin),R3Point(xmax,ymax,zmin)));
	if (pos.Z() >= zmax) 
		box_sides.push_back(quadrilateral(R3Point(xmin,ymin,zmax),R3Point(xmax,ymax,zmax)));

	//if (box_sides.size() < 1) return R2black_pixel;
	quadrilateral side = box_sides[rand_int(box_sides.size())-1];
	pos2 = side.sample();
	pdf = 0;
	for (int i=0; i<box_sides.size(); i++)
		pdf += box_sides[i].visibility(pos)*box_sides[i].area();
	pdf = 1./pdf;
	normal2 = R3Vector(0,0,0);
	normal2[side.fixed_dim] = 1;
	if (side.min[side.fixed_dim] == min(box->Min()[side.fixed_dim],box->Max()[side.fixed_dim])) normal2 *= -1;



/*
#define nearer(pt,min,max) (min-pt > pt-max ? max : min)
#define farther(pt,min,max) (min-pt > pt-max ? min : max)
	double nearer_coord, farther_coord;
	if (pos.X() < xmin || pos.X() > xmax) {
		nearer_coord = nearer(pos.X(),xmin,xmax);
		farther_coord = farther(pos.X(),xmin,xmax);
		pos2 = R3Point(nearer_coord,rand_interp(ymin,ymax),rand_interp(zmin,zmax));
		normal2 = 
	}

*/

/*
	//if point is inside box (all sides visible)
	int dim_fixed = rand_int(3)-1;
	//printf("%d\n",dim_fixed);
	int other_dim_1 = 0, other_dim_2 = 1;
	if (dim_fixed==0) {
		other_dim_1 = 1; other_dim_2 = 2;
	} else if (dim_fixed==1) {
				other_dim_1 = 0; other_dim_2 = 2;
	}

	double max_a_coord = max(box->Max()[other_dim_1],box->Min()[other_dim_1]);
	double min_a_coord = min(box->Max()[other_dim_1],box->Min()[other_dim_1]);
	double max_b_coord = max(box->Max()[other_dim_2],box->Min()[other_dim_2]);
	double min_b_coord = min(box->Max()[other_dim_2],box->Min()[other_dim_2]);

	double fixed_coord = box->Max()[dim_fixed];
	normal2 = R3Vector(0,0,0);
	normal2[dim_fixed] = box->Max()[dim_fixed] - box->Min()[dim_fixed];
	normal2.Normalize();
	if (rand_unit() > .5)
	{
		fixed_coord = box->Min()[dim_fixed];
		normal2 *= -1;
	}
	double a_coord = rand_interp(min_a_coord,max_a_coord);
	double b_coord = rand_interp(min_b_coord,max_b_coord);

	pos2[dim_fixed] = fixed_coord;
	pos2[other_dim_1] = a_coord;
	pos2[other_dim_2] = b_coord;

	pdf = abs(box->XMax()-box->XMin())*abs(box->YMax()-box->YMin())*abs(R3zaxis_vector.Dot(normal)) +
				abs(box->XMax()-box->XMin())*abs(box->ZMax()-box->ZMin())*abs(R3yaxis_vector.Dot(normal)) +
				abs(box->ZMax()-box->ZMin())*abs(box->YMax()-box->YMin())*abs(R3xaxis_vector.Dot(normal));
	pdf = 1.f/(2.f*pdf);
	//normal2.Print(); printf("\n");
*/
	return true;
}


//change frame
R3Vector trans_frame(R3Vector res, const R3Vector &normal) {
	R3Vector tmp(res);
	for (;;) {
		double x = rand_unit();
		if (x<1.f/3.f) tmp.SetX(0);
		else if (x<2.f/3.f) tmp.SetY(0);
		else tmp.SetZ(0);
		if (tmp != R3zero_vector) break;
		else tmp = res;
	}
	R3Vector u(tmp);
	u.Cross(res);
	R3Vector v(u);
	v.Cross(res);
	R3Matrix trans(u,v,normal);
	res.Transform(trans);
	res.Normalize();
	return res;
}

R3Vector sample_cone(const R3Vector &normal, const double &theta_max, double &pdf) {
	double u1 = rand_unit();
	double u2 = rand_unit();

	pdf = 1.f/(2.f*M_PI*(1.f - cosf(theta_max)));
	double costheta = rand_interp(cosf(theta_max),1);
	double sintheta = sqrtf(1.f-costheta*costheta);
	float phi = rand_interp(0,2*M_PI);
	R3Vector res(cosf(phi)*sintheta, sinf(phi)*sintheta, costheta);

	return trans_frame(res,normal);

}
//sampling direction on hemisphere at pos, with "north pole" aligned with normal
//FIX TRANSFORMATION OF FRAME--BUGGY!
R3Vector sample_hemisphere(const R3Vector &normal, double &pdf) {
	pdf = 1.f/(2.f*M_PI);
	double u1 = rand_unit();
	double u2 = rand_unit();
	double phi = 2*M_PI*u2;
	double r  = sqrt(max(0.f,1-u1*u1));
	double x = cos(phi)*r;
	double y = sin(phi)*r;
	double z = u1;

	R3Vector delta = R3Vector(x,y,z) - R3zaxis_vector;
	return normal+delta;
}

R3Vector sample_hemisphere_cos(const R3Vector &normal, double &pdf) {

	double u1 = rand_unit();
	double u2 = rand_unit();

	double r = sqrt(u1);
	double theta = 2.f*M_PI*u2;
	double x = r*cos(theta);
	double y = r*sin(theta);
	double z = sqrt(max(0.f,1.f-x*x-y*y));
	R3Vector res(x,y,z);
/*
	R3Vector tmp(res);
	for (;;) {
		double x = rand_unit();
		if (x<1.f/3.f) tmp.SetX(0);
		else if (x<2.f/3.f) tmp.SetY(0);
		else tmp.SetZ(0);
		if (tmp != R3zero_vector) break;
		else tmp = res;
	}

	R3Vector u(tmp);
	u.Cross(res);
	R3Vector v(u);
	v.Cross(res);
	R3Matrix trans(u,v,normal);
	res.Transform(trans);
	res.Normalize();
	//assert(res.IsNormalized());
//	R3Vector delta = res - R3zaxis_vector;
	*/
	res = trans_frame(res,normal);

	pdf = abs(normal.Dot(res));
	//if (pdf>1) printf("\n%g\n",pdf);
	//assert (pdf > 0 && pdf <= 1.f);

//	return normal+delta;
	return res;
}

bool sample(const R3Sphere *sphere, const R3Point &pos, R3Point &pos2, R3Vector &normal2, double &pdf) {


	double r1 = rand_unit(), r2 = rand_unit();
	double x0 = cos(2*M_PI*r2)*2*sqrt(r1*(1-r1));
	double y0 = sin(2*M_PI*r2)*2*sqrt(r1*(1-r1));
	double z0 = 1-2*r1;
	pos2 = R3Point(x0,y0,z0)*sphere->Radius();

	normal2 = pos2.Vector();
	normal2.Normalize();
	pos2 = sphere->Center() + pos2;

	pdf = 4*M_PI*pow(sphere->Radius(),2);
	pdf = 1./pdf;

	return true;
}

bool sample(const R3Cylinder *cylinder, const R3Point &pos, R3Point &pos2, R3Vector &normal2, double &pdf) {
	//h is a random  point on the cyl axis
	R3Point h = cylinder->Axis().Start() + (rand_unit()*cylinder->Height()*cylinder->Axis().Vector());

	//find random vector u orthogonal to axis a
	R3Vector a(cylinder->Axis().Vector()),tmp,e1,e2,u1,u2,u;
	tmp = a!=R3zaxis_vector ? R3zaxis_vector : R3yaxis_vector;
	e1 = tmp; e1.Cross(a); e1.Normalize();
	e2 = e1; e2.Cross(a); e2.Normalize();
	double theta = 2*M_PI*rand_unit();
	u1 = cos(theta)*e1;
	u2 = sin(theta)*e2;
	u = u1+u2;
	//printf("%g\n",u.Length());
	assert(u.IsNormalized());
	
/*
	double costheta = cosf(theta);
	double sintheta = sqrt(1-costheta*costheta);
	double u0 = theta;//acosf(costheta);
	double u1 = asinf(sintheta);
	if (a[2] != 0) u = R3Vector(u0,u1,(-u0*a[0]-u1*a[1])/a[2]);
	else if (a[0] != 0) u = R3Vector((-u0*a[1]-u1*a[2])/a[0],u0,u1);
	else u = R3Vector(u0,(-u0*a[0]-u1*a[2])/a[1],u1);
	u.Normalize();
*/
	normal2 = u;
	pos2 = h + cylinder->Radius()*u;
	pdf = 1/(2*M_PI*cylinder->Radius()*cylinder->Height());

	return true;
}

/*
double area(R3Shape *shape) {
			R3Box *box = shape->box;
	switch(shape->type) {
		case R3_SPHERE_SHAPE:
		case R3_BOX_SHAPE:
			return abs(box->XMax()-box->XMin())*abs(box->YMax()-box->YMin()) +
				abs(box->XMax()-box->XMin())*abs(box->ZMax()-box->ZMin()) +
				abs(box->ZMax()-box->ZMin())*abs(box->YMax()-box->YMin());
		default:
			return -1;
	}
}
*/

R2Pixel brdf(R3Node *node) {
	return node->material->kd/(M_PI);
}

R2Pixel brdf_fresnel_cond(double cosi,double eta, double k) {
	if (cosi==0) return R2black_pixel;
	cosi = fabsf(cosi);
	double tmp = (eta*eta+k*k);
	double rpar = (tmp*cosi*cosi-2*eta*cosi+1)/(tmp*cosi*cosi+2*eta*cosi+1);
	double rperp = (tmp-2*eta*cosi+cosi*cosi)/(tmp+2*eta*cosi+cosi*cosi);
	double res = (rpar + rperp)/2.;
	return R2Pixel(res,res,res,1)/cosi;
}

R2Pixel brdf_fresnel_diel(double cosi,double cost,double etai,double etat){
	cosi=fabsf(cosi); cost=fabsf(cost);
	double rpar = ((etat*cosi)-(etai*cost))/((etat*cosi)+(etai*cost));
	double rperp = ((etai*cosi)-(etat*cost))/((etai*cosi)+(etat*cost));
	double res = (rpar*rpar+rperp*rperp)/2.;
	return R2Pixel(res,res,res,1)/cosi;	
}

//#define MAX_DEPTH 3
bool within_cone(const R3Point &pos, const R3Point &head, const R3Vector &axis, const double &theta) {
	assert(axis.IsNormalized());
	R3Vector v(pos-head);
	double adj_len = axis.Dot(v);
	double hyp_len = v.Length();
	return acosf(adj_len/hyp_len) < theta;
}

R2Pixel GetColor(const R3Scene *scene, const R3Ray &ray, int depth, double cutoff) {

	double t,t_tmp,pdf;
	R3Point pos,pos2,pos_tmp;
	R3Vector normal,normal2,normal_tmp;
	R3Node *node,*node_tmp;

	double epsilon = .01;
	R3Vector reflect, refract;
	R3Ray rand_ray;


	if (!R3Intersects(scene,ray,node,pos,t,normal) || !node->material)
		return R2black_pixel;
		//return R2Pixel(.7,.7,1,1);
	bool inside = normal.Dot(ray.Vector()) > 0;

	R2Pixel emission = node->material->emission;
	double max_emission = max(emission.Red(),max(emission.Green(),emission.Blue()));
	reflect = reflect_vector(ray.Vector(),normal);
	/*
	if (depth <= 0)
		if (rand_unit() < max_emission)
			emission /= max_emission;
		else return R2black_pixel;
*/
	if (depth<=0) return R2black_pixel;

/*
	R2Pixel L_d(0,0,0,0);
	R3Vector wi_L;
	double t_L,L_i;
	R2Pixel K_D = node->material->kd;
	//calculate direct lighting, using MC estimation
	//for each light obtain sample
	for (int i=0; i<scene->NLights(); i++) {
	//test sample for visibility
		//if (sample(scene,scene->Light(i),pos,wi_L,t_L,L_i))
	//if visible, add contribution to estimate
			//L_d += (K_D/M_PI)*
	}
*/

	
	
	//direct lighting
	R2Pixel L_d(R2black_pixel);
	R3Vector wi;
	for (int i=0; i<scene->NLights(); i++) {
		R3Light *light = scene->Light(i);
		//skip directional/point lights, or if this light is the same as shape being modeled
		/*if (light->type == R3_SPOT_LIGHT) {
			if (within_cone(pos,light->position,light->direction,light->angle_cutoff)) {
				wi = sample_cone(light->direction,light->angle_cutoff,pdf);
				double attenuation = 
			}
			
		} else*/ {
		if (scene->Light(i)->type != R3_AREA_LIGHT || scene->Light(i)->node==node) continue;
		switch (light->node->shape->type) {
			case R3_SPHERE_SHAPE:
				sample(light->node->shape->sphere,pos,pos2,normal2,pdf);
				break;
			case R3_BOX_SHAPE:
				sample(light->node->shape->box,pos,normal,pos2,normal2,pdf);
				//while (normal2.Dot(pos-pos2)<0);
				break;
			case R3_CYLINDER_SHAPE:
				sample(light->node->shape->cylinder,pos,pos2,normal2,pdf);
				break;
			default:
				assert(false);
		}
		wi = pos2-pos;
		//check visibility
		R3Ray wi_ray(pos+epsilon*wi,wi);
		R3Intersects(scene,wi_ray,node_tmp,pos_tmp,t_tmp,normal_tmp);
		if (node_tmp != light->node) 
			continue;
		//MC estimate
		double costheta = inside ? -normal.Dot(wi) : normal.Dot(wi);
		double costheta2 = normal2.Dot(-wi);
			R2Pixel f = brdf(node);
		L_d += f*light->node->material->emission*abs(costheta)/pow(wi.Length(),2)/pdf;
		}
	}//DBG!!!!!!!!!!!!!!!!!
	L_d.Clamp();

	//indirect lighting
	//generate a ray in a new random diretion; notation is shirley ch. 16
	/*if (inside) normal *= -1;
		double r1 = rand_unit(), r2 = rand_unit();
		R3Vector a(cos(2*M_PI*r1)*sqrt(r2),sin(2*M_PI*r1)*sqrt(r2),sqrt(1-r2));
		R3Vector w = normal;
		R3Vector s = w;
		if (s.X() <= s.Y() && s.X() <= s.Z()) s.SetX(s.X()-1); 
		else if (s.Y() <= s.Z()) s.SetY(s.Y()-1); else s.SetZ(s.Z()-1);
		R3Vector u = s;
		u.Cross(w); u.Normalize();
		R3Vector v = w;
		v.Cross(u);
		R3Vector row1(u.X(),v.X(),w.X());
		R3Vector row2(u.Y(),v.Y(),w.Y());
		R3Vector row3(u.Z(),v.Z(),w.Z());
		a = R3Vector(a.Dot(row1),a.Dot(row2),a.Dot(row3));
		rand_ray = R3Ray(pos+epsilon*a,a);
		pdf = normal.Dot(rand_ray.Vector())/M_PI;
	*/
	 rand_ray = R3Ray(pos,sample_hemisphere_cos(inside ? -normal : normal,pdf));
		R2Pixel L_e = emission;
		R2Pixel R = node->material->kd;
		R2Pixel f = brdf(node);//DBG!!!!!!!!!!!!!!
//return L_e + node->material->kd*GetColor(scene,rand_ray,depth-1,cutoff);
		R2Pixel L_i = f*GetColor(scene,rand_ray,depth-1,cutoff)*
			(inside ? -normal : normal).Dot(rand_ray.Vector())/pdf;
				//if (inside) normal *= -1;

		//	node->material->kt*GetColor(scene,R3Ray(pos,reflect),depth-1,cutoff);

	//specular reflection
		R2Pixel L_s;
		if (node->material->ks != R2Pixel(0,0,0,0)) {
		R3Vector wi = reflect_vector(-ray.Vector(),normal);
		wi.Normalize();
		double cosi = normal.Dot(wi);
		R3Ray wi_ray(pos+wi,wi);
		L_s = /*brdf_fresnel_cond(cosi,2.485,3.433)**/GetColor(scene,R3Ray(pos+epsilon*wi,wi),depth-1,cutoff)*abs(cosi);
		/*if (!(L_s.Red()==0 && L_s.Green()==0 && L_s.Blue()==0)) {
			wi.Print(); printf("\n"); assert(false);
		}*/
		}

		//transmitted light
		R2Pixel L_t(R2black_pixel);
		if (node->material->kt != R2Pixel(0,0,0,0)) {
			R3Vector wt;
			double ir = node->material->indexofrefraction;
			//if entering use given idx of refraction, otherwise reciprocal
			if (ray.Vector().Dot(normal) < 0)
				wt = refract_vector(ray.Vector(),normal,ir);
			else
				wt = refract_vector(ray.Vector(),-normal,1./ir);
			wt.Normalize();
			double cost = fabsf(wt.Dot(normal));
			double cosi = fabsf(ray.Vector().Dot(normal));
			R3Ray wt_ray(pos+wt*epsilon,wt);
			L_t = /*(R2white_pixel-brdf_fresnel_diel(cosi,cost,ir,1./ir))**/GetColor(scene,wt_ray,depth-1,cutoff)*abs(cosi);
		}

		R2Pixel L_o = L_e;//MAX_DEPTH == depth ? L_e : R2black_pixel;
		L_o += /*L_d +*/ L_i + node->material->ks*L_s + node->material->kt*L_t;
		L_o.Clamp();
		return L_o;



	//} else
	//	return R2black_pixel;

}	

//add emissive objects to list of lights in scene
void init_lights(vector<R3Light *> &lights, R3Node *node) {
	//check if this node's shape is emissive
	if (node->shape)
		if (node->material->emission != R2Pixel(0,0,0,0)) {
			R3Light *new_light = new R3Light();
			new_light->node = node;
			new_light->type = R3_AREA_LIGHT;
			lights.push_back(new_light);
		}
	//recurse on children
		for (int i=0; i<node->children.size(); i++)
			init_lights(lights,node->children[i]);
}


// Create image from scene
// This is the main ray tracing function 

R2Image *RenderImage(R3Scene *scene, int width, int height, int max_reflections, double min_luminance, int n_samples)
{


n_samples=100;																	;

















	// Allocate  image
	R2Image *image = new R2Image(width, height);
	if (!image) {
		fprintf(stderr, "Unable to allocate image\n");
		return NULL;
	}

	R3Point eye = scene->camera.eye;
	double d = scene->camera.fardist;
	double xfov = scene->camera.xfov;
	double yfov = scene->camera.yfov;
	R3Vector towards = scene->camera.towards;
	R3Vector up = scene->camera.up; 
	R3Vector right = scene->camera.right; 


	double t;
	R3Point pos;
	R3Vector normal;
	R3Node *node;

	init_lights(scene->lights,scene->root);
	GenEndpoints(eye,d,up,towards,right,xfov,yfov);










	R3Point pos2,pos_tmp; R3Vector normal2,normal_tmp; double pdf;
//sample_hemisphere_cos(R3Vector(1,0,0),pdf);
	R3Node *node_tmp;R2Pixel dbgpix;
	R3Box *box = new R3Box(R3Point(-3,-.5,3.5),R3Point(-6,.5,4.5));
	bool db = R3Intersects(scene,R3Ray(R3Point(eye),R3Point(0,0,.5)-eye),node_tmp,pos_tmp,t,normal_tmp);
	//sample(box,R3Point(-2,0,0),pos2,normal2,pdf);
	dbgpix = GetColor(scene,R3Ray(R3Point(eye),R3Point(0,0,3)-eye),7,.4);
	dbgpix = GetColor(scene,R3Ray(R3Point(-2,-.5,1),R3Vector(0,0,-1)),1,.4);



/*
	R3Ray dbg_ray = R3xaxis_ray;
	for (int i=1; i<20; i++) {
	R2Pixel dbg_pix = GetColor(scene,dbg_ray,i,0);
	printf("%d path segments: %g,%g,%g\n",i,dbg_pix.Red(),dbg_pix.Green(),dbg_pix.Blue());
	} //shriley's dbg recommendation--w/ test4.scn. refl converges to 1.
*/
//	R3Ray dbg_ray = R3xaxis_ray;
	for (int i=0; i<image->Width(); i++) {
		for (int j=0; j<image->Height(); j++) {
			for (int k=0; k<n_samples; k++) {
				//printf("%d, %d: ",i,j);
				R3Ray ray = GenRay(i,j,k,sqrt((double)n_samples),eye,up,towards,xfov,yfov,width,height);
				image->Pixel(i,j) += GetColor(scene,ray,5,0);
			}
			image->Pixel(i,j) /= n_samples;
			image->Pixel(i,j).Clamp();
		}
		fprintf(stderr,"\r%5.2f%% completed",(i*image->Width())*100.0/(image->Width()*image->Height()));
	}


	// Return image
	return image;

 }
