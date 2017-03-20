// Include file for mesh class




////////////////////////////////////////////////////////////
// DEPENDENCY INCLUDE FILES
////////////////////////////////////////////////////////////

#include <vector>
using namespace std;



////////////////////////////////////////////////////////////
// MESH VERTEX DECLARATION
////////////////////////////////////////////////////////////

struct R3MeshVertex {
  // Constructors
  R3MeshVertex(void);
  R3MeshVertex(const R3MeshVertex& vertex);
  R3MeshVertex(const R3Point& position, const R3Vector& normal, const R2Point& texcoords);

  // Data
  R3Point position;
  R3Vector normal;
  R2Point texcoords;
  int id; 
};



////////////////////////////////////////////////////////////
// MESH FACE DECLARATION
////////////////////////////////////////////////////////////

struct R3MeshFace {
  // Constructors
  R3MeshFace(void);
  R3MeshFace(const R3MeshFace& face);
  R3MeshFace(R3MeshVertex *v1, R3MeshVertex *v2, R3MeshVertex *v3);

  // Data
  R3MeshVertex *vertices[3];
  R3Plane plane;
  int id;
};



////////////////////////////////////////////////////////////
// MESH CLASS DECLARATION
////////////////////////////////////////////////////////////

struct R3Mesh {
  // Constructors
  R3Mesh(void);
  R3Mesh(const R3Mesh& mesh);
  ~R3Mesh(void);

  // Vertex and face access functions
  int NVertices(void) const;
  R3MeshVertex *Vertex(int k) const;
  int NFaces(void) const;
  R3MeshFace *Face(int k) const;

  // Warps
  void RandomNoise(double noise);
  void Translate(double dx, double dy, double dz);
  void Scale(double sx, double sy, double sz);
  void Rotate(double angle, const R3Line& axis);
  void Inflate(double offset);
  void Deform(R3Point *source_points, R3Point *target_points, int npoints, double t);
  void Fun(void);

  // Filters
  void Smooth(double sigma);
  void SmoothBilateral(double sigma);
  void Sharpen(double sigma);
  void Truncate(double t);
  void Bevel(double t);

  // Remeshing
  void Subdivide(void);
  void SubdivideLoop(void);
  void SplitLongEdges(double max_edge_length);
  void ClusterVertices(double grid_cell_size);
  void RemeshUniform(void);
  void RemeshRegular(double grid_cell_size);

  // Analysis
  void ComputeVertexNormals(void);
  void ComputeVertexCurvatures(void);

  // Topological fixup
  void FixHoles(void);
  void FixCracks(double epsilon);
  void FixIntersections(void);

  // Geometry construction
  void Fractalize(int nlevels);
  void SurfaceSweep(const R3Mesh& crosssection_polygon, const R3Mesh& centerline_curve);
  void SurfaceOfRevolution(const R3Mesh& profile_curve, 
    const R3Line& axis_of_revolution, double rotation_angle_step);

  // Boolean operations
  void Intersect(const R3Mesh& mesh);
  void Subtract(const R3Mesh& mesh);
  void Union(const R3Mesh& mesh);
  void Crop(const R3Plane& plane);

  // File input/output 
  int Read(const char *filename);
  int ReadRay(const char *filename);
  int ReadOff(const char *filename);
  int ReadObj(const char *filename);
  int ReadImage(const char *filename);
  int Write(const char *filename);
  int WriteRay(const char *filename);
  int WriteOff(const char *filename);
  int WriteObj(const char *filename);

  // Utility functions
  R3MeshVertex *CreateVertex(const R3Point& position, 
    const R3Vector& normal, const R2Point& texcoords);
  R3MeshFace *CreateFace(R3MeshVertex *v1, R3MeshVertex *v2, R3MeshVertex *v3);
  void DeleteVertex(R3MeshVertex *vertex);
  void DeleteFace(R3MeshFace *face);

  // Output functions
  void Draw(void) const;
  void Outline(void) const;
  void Print(FILE *fp = stdout) const;

  // Data
  vector<R3MeshVertex *> vertices;
  vector<R3MeshFace *> faces;
  R3Box bbox;
};



////////////////////////////////////////////////////////////
// MESH INLINE FUNCTIONS
////////////////////////////////////////////////////////////

inline int R3Mesh::
NVertices(void) const
{
  // Return number of vertices in mesh
  return vertices.size();
}



inline R3MeshVertex *R3Mesh::
Vertex(int k) const
{
  // Return kth vertex of mesh
  return vertices[k];
}



inline int R3Mesh::
NFaces(void) const
{
  // Return number of faces in mesh
  return faces.size();
}



inline R3MeshFace *R3Mesh::
Face(int k) const
{
  // Return kth face of mesh
  return faces[k];
}



