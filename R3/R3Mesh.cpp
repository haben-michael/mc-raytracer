// Source file for mesh class



// Include files

#include "R3.h"



////////////////////////////////////////////////////////////
// MESH CONSTRUCTORS/DESTRUCTORS
////////////////////////////////////////////////////////////

R3Mesh::
R3Mesh(void)
  : bbox(R3null_box)
{
}



R3Mesh::
R3Mesh(const R3Mesh& mesh)
  : bbox(R3null_box)
{
  // Create vertices
  for (int i = 0; i < mesh.NVertices(); i++) {
    R3MeshVertex *ov = mesh.Vertex(i);
    CreateVertex(ov->position, ov->normal, ov->texcoords);
    ov->id = i;
  }

  // Create faces
  for (int i = 0; i < mesh.NFaces(); i++) {
    R3MeshFace *of = mesh.Face(i);
    R3MeshVertex *nv0 = Vertex(of->vertices[0]->id);
    R3MeshVertex *nv1 = Vertex(of->vertices[1]->id);
    R3MeshVertex *nv2 = Vertex(of->vertices[2]->id);
    CreateFace(nv0, nv1, nv2);
  }
}



R3Mesh::
~R3Mesh(void)
{
  // Delete faces
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *f = Face(i);
    delete f;
  }

  // Delete vertices
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *v = Vertex(i);
    delete v;
  }
}



////////////////////////////////////////////////////////////
// MESH PROCESSING FUNCTIONS
////////////////////////////////////////////////////////////

void R3Mesh::
RandomNoise(double noise)
{
  // Add noise of a random amount and direction 
  // to the position of every vertex, where the 
  // input parameter "noise" specifies the maximum 
  // displacement of any vertex

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}



void R3Mesh::
Translate(double dx, double dy, double dz)
{
  // Translate the mesh by adding a 
  // vector (dx,dy,dz) to every vertex

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
Scale(double sx, double sy, double sz)
{
  // Scale the mesh by increasing the distance 
  // from every vertex to the origin by a factor 
  // given for each dimension (sx, sy, sz)

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
Rotate(double angle, const R3Line& axis)
{
  // Rotate the mesh counter-clockwise by an angle 
  // (in radians) around a line axis
  // Hint: use the Rotate function in the R3Point class

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
Inflate(double offset)
{
  // Move every vertex by a given offset along its 
  // normal direction (offsets can be negative)

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
Deform(R3Point *source_points, R3Point *target_points, int npoints, double t)
{
  // Warp the mesh by a deformation defined by point correspondences.  
  // The position of every vertex in the output mesh should be 
  // determined by a Gaussian weighted average of the offsets implied 
  // by the point correspondences (x,y,z) -> (x',y',z').  Use a reasonable
  // heuristic to choose sigma for the Gaussian.

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
Fun(void)
{
  // Warp a mesh using a non-linear mapping of your choice 
  // (examples are sine, bulge, swirl)

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
Smooth(double sigma)
{
  // Smooth the mesh by moving every vertex to a position 
  // determined by a weighted average of its immediate neighbors 
  // (with weights determined by a Gaussian function of distance)

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
Sharpen(double sigma)
{
  // Accentuate details in the mesh by moving every vertex away 
  // from the position determined by a weighted average of its neighbors 
  // (with weights determined by a Gaussian function of distance).  
  // This filter moves vertices in the direction opposite from smoothing.

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
Truncate(double t)
{
  // For every vertex, create a new vertex a parameter t [0-1] 
  // of the way along each of its attached edges, and then 
  // "chop off" the pyramid whose base is formed by the new vertices 
  // and whose apex is the original vertex, creating a set of faces 
  // to triangulate the hole.

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
Bevel(double t)
{
  // For every edge, create a new face whose vertices are t [0-1] 
  // of the way along each of its attached edges.  This requires 
  // first truncating all vertices by t, creating new vertices t [0-1] 
  // of the way along each of new edges, and then "chopping off" a 
  // prism for each of the original edges, creating a pair of faces 
  // to triangulate the hole.

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
SmoothBilateral(double sigma)
{
  // Smooth the mesh using a bilateral filter as in 
  // [Jones et al, Siggraph 2003] or 
  // [Fleishman et al., Siggraph 2003]

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
Subdivide(void)
{
  // Replace every triangle by four new ones according to the 
  // Loop subdivision scheme.   This requires creating a new vertex 
  // at the midpoint of every edge, removing the original triangles, 
  // and creating four new triangles to replace each of the original ones.  

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
SubdivideLoop(void)
{
  // Replace every triangle by four new ones according to the 
  // Loop subdivision scheme.   This requires creating a new vertex 
  // at the midpoint of every edge, removing the original triangles, 
  // and creating four new triangles to replace each of the original ones.  
  // The positions of all vertices should be updated according to 
  // the Loop subdivision weights

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}



void R3Mesh::
SplitLongEdges(double max_edge_length)
{
  // Iteratively split edges longer than max_edge_length.  
  // Note that every edge split produces a new vertex at the 
  // edge midpoint and replaces the two adjacent faces with four.  
  // Edges should be split repeatedly until there is none longer 
  // than the given threshold.  Note: an extra point will be given if 
  // longer edges are split first (with produces better shaped triangles).

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
ClusterVertices(double grid_cell_size)
{
  // Simplify the mesh by clustering vertices residing in the same 
  // cell of a grid defined by x, y, and z spacing parameters.  
  // All vertices within the same grid cell should be merged 
  // into a single vertex, that vertex should be placed at the 
  // centroid of the cluster vertices, and all edges and faces 
  // that collapse as a result of the vertex merging should be removed. 

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
RemeshUniform(void)
{
  // Generate a new mesh whose vertices are evenly spaced 
  // (i.e., all edges are the same length).  One way to accomplish 
  // this is to iteratively adjust vertices so that they are at 
  // the position on the surface closest to the centroid of its neighbors. 

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
RemeshRegular(double grid_cell_size)
{
  // Generate a new mesh whose vertices are on regularly 
  // spaced intervals in x, y, and z.  One way to accomplish 
  // this is to rasterize the interior of the mesh into a 
  // voxel grid with negative values inside and positive values outside 
  // and then compute the zero-set isosurface from the grid

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
ComputeVertexNormals(void)
{
  // Compute the surface normal at every vertex by taking a 
  // weighted average of the normals for the attached faces, 
  // where the weights are determined by the areas of the faces.

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
ComputeVertexCurvatures(void)
{
  // Compute the minimum and maximum curvatures of the surface at 
  // every vertex using the method described in [Rusinkiewicz, 3DPVT 2004].  
  // Store the resulting two curvatures values in the texture coordinates 
  // for every vertex

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
Fractalize(int nlevels)
{
  // Create fractal detail by recursively subdividing every triangle of this
  // mesh into four (as in the Subdivide feature above) and moving 
  // every new vertex by a random offset vector (as in the Random Noise feature above) 
  // with magnitude proportional to the original edge length.  The subdivision
  // and random displacement should be applied recursively for nlevels.  Note that the 
  // output mesh may resemble a mountain if the input is planar.

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}



void R3Mesh::
SurfaceSweep(const R3Mesh& crosssection_polygon, const R3Mesh& centerline_curve)
{
  // Create new vertices and faces by sweeping a polygon along a curve.  
  // The vertices representing a cross-section polygon are provided in 
  // the first input mesh file, and the vertices representing the sweep 
  // centerline curve are provided in the second mesh file (for both, take 
  // the vertices of the meshes in order and ignore the faces).  New vertices 
  // should be created by successively translating and rotating the vertices 
  // of the cross-section polygon to match the position and orientation of 
  // vertices/edges in the centerline, and new faces should be constructed 
  // by connecting adjacent vertices created during the sweep.  
  // Note: an extra 3 points will be awarded if your implementation avoids 
  // self-intersecting polygons in all cases.

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
SurfaceOfRevolution(const R3Mesh& profile_curve, 
  const R3Line& axis_of_revolution, double rotation_angle_step)
{
  // Add new vertices and faces to the mesh by sweeping a profile curve 
  // around an axis of revolution.  The vertices representing the profile 
  // curve are provided in the passed mesh file (take the vertices of the 
  // mesh in order and ignore the faces).  The axis of revolution and 
  // rotation angle step size are provided in the arguments.  New vertices 
  // should be created by successively rotating the original vertices around 
  // the axis by the step size and new faces should be constructed by 
  // connecting adjacent vertices to create a surface of revolution.

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}



void R3Mesh::
FixHoles(void)
{
  // Create triangles covering the holes of a mesh by connecting vertices 
  // on the boundary of every hole.  You should completely cover the hole, 
  // while doing your best to produce well-shaped triangles 
  // (e.g., by connecting closer vertices first).  
  // However, do not worry about triangles intersecting other parts of 
  // the mesh in your implementation.

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
FixCracks(double epsilon)
{
  // Merge boundary vertices and edges within a specified 
  // distance (epsilon) of one another.

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
FixIntersections(void)
{
  // Insert edges at face-face intersections and discard 
  // the smaller part of the mesh "pinched" off by new edge loops.  
  // Note: this is hard.

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
Intersect(const R3Mesh& mesh)
{
  // Intersect the solid implied by this mesh with another, 
  // keeping only the faces enclosing the intersection of the two solids.
  // This feature requires introducing edges at every face intersection 
  // and removing parts of the mesh that lie in the exterior of the 
  // solid object implied by either of the two meshes. 

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
Subtract(const R3Mesh& mesh)
{
  // Subtract the solid implied by this mesh with another, 
  // keeping only the faces enclosing the difference of the two solids.
  // This feature requires introducing edges at every face intersection 
  // and removing parts of the mesh that lie in the interior of the 
  // solid object implied by the passed mesh.

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
Union(const R3Mesh& mesh)
{
  // Union the solid implied by this mesh with another, 
  // keeping only the faces enclosing the union of the two solids.
  // This feature requires introducing edges at every face intersection 
  // and removing parts of the mesh that lie in the interior of the 
  // solid object implied by both of the two meshes. 

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




void R3Mesh::
Crop(const R3Plane& plane)
{
  // Crop the input mesh to the positive side of the plane.  
  // This feature requires clipping each polygon crossing the plane, 
  // and discarding any part of any face on the negative side of the plane.

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");
}




////////////////////////////////////////////////////////////
// MESH ELEMENT CREATION/DELETION FUNCTIONS
////////////////////////////////////////////////////////////

R3MeshVertex *R3Mesh::
CreateVertex(const R3Point& position, const R3Vector& normal, const R2Point& texcoords)
{
  // Create vertex
  R3MeshVertex *vertex = new R3MeshVertex(position, normal, texcoords);

  // Update bounding box
  bbox.Union(position);

  // Set vertex ID
  vertex->id = vertices.size();

  // Add to list
  vertices.push_back(vertex);

  // Return vertex
  return vertex;
}



R3MeshFace *R3Mesh::
CreateFace(R3MeshVertex *v1, R3MeshVertex *v2, R3MeshVertex *v3)
{
  // Create face
  R3MeshFace *face = new R3MeshFace(v1, v2, v3);

  // Set face ID
  face->id = faces.size();

  // Add to list
  faces.push_back(face);

  // Return face
  return face;
}



void R3Mesh::
DeleteVertex(R3MeshVertex *vertex)
{
  // Remove vertex from list
  for (unsigned int i = 0; i < vertices.size(); i++) {
    if (vertices[i] == vertex) {
      vertices[i] = vertices.back();
      vertices[i]->id = i;
      vertices.pop_back();
      break;
    }
  }

  // Delete vertex
  delete vertex;
}



void R3Mesh::
DeleteFace(R3MeshFace *face)
{
  // Remove face from list
  for (unsigned int i = 0; i < faces.size(); i++) {
    if (faces[i] == face) {
      faces[i] = faces.back();
      faces[i]->id = i;
      faces.pop_back();
      break;
    }
  }

  // Delete face
  delete face;
}



////////////////////////////////////////////////////////////////////////
// I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3Mesh::
Read(const char *filename)
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .ply)\n", filename);
    return 0;
  }

  // Read file of appropriate type
  if (!strncmp(extension, ".ray", 4)) 
    return ReadRay(filename);
  else if (!strncmp(extension, ".obj", 4)) 
    return ReadObj(filename);
  else if (!strncmp(extension, ".off", 4)) 
    return ReadOff(filename);
  else if (!strncmp(extension, ".jpg", 4)) 
    return ReadImage(filename);
  else if (!strncmp(extension, ".jpeg", 4)) 
    return ReadImage(filename);
  else if (!strncmp(extension, ".bmp", 4)) 
    return ReadImage(filename);
  else if (!strncmp(extension, ".ppm", 4)) 
    return ReadImage(filename);
  else {
    fprintf(stderr, "Unable to read file %s (unrecognized extension: %s)\n", filename, extension);
    return 0;
  }
}



int R3Mesh::
Write(const char *filename)
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .ply)", filename);
    return 0;
  }

  // Write file of appropriate type
  if (!strncmp(extension, ".ray", 4)) 
    return WriteRay(filename);
  else if (!strncmp(extension, ".obj", 4)) 
    return WriteObj(filename);
  else if (!strncmp(extension, ".off", 4)) 
    return WriteOff(filename);
  else {
    fprintf(stderr, "Unable to write file %s (unrecognized extension: %s)", filename, extension);
    return 0;
  }
}



////////////////////////////////////////////////////////////
// IMAGE FILE INPUT/OUTPUT
////////////////////////////////////////////////////////////

int R3Mesh::
ReadImage(const char *filename)
{
  // Create a mesh by reading an image file, 
  // constructing vertices at (x,y,luminance), 
  // and connecting adjacent pixels into triangles. 
  // That is, the image is interpretted as a height field, 
  // where the luminance of each pixel provides its z-coordinate.

  // Read image
  R2Image *image = new R2Image();
  if (!image->Read(filename)) return 0;

  // FILL IN IMPLEMENTATION HERE
  fprintf(stderr, "Function not implemented\n");

  // Delete image
  delete image;

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////
// OFF FILE INPUT/OUTPUT
////////////////////////////////////////////////////////////

int R3Mesh::
ReadOff(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s\n", filename);
    return 0;
  }

  // Read file
  int nverts = 0;
  int nfaces = 0;
  int nedges = 0;
  int line_count = 0;
  int vertex_count = 0;
  int face_count = 0;
  char buffer[1024];
  char header[64];
  while (fgets(buffer, 1023, fp)) {
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) bufferp++;

    // Skip blank lines and comments
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;

    // Check section
    if (nverts == 0) {
      // Read header keyword
      if (strstr(bufferp, "OFF")) {
        // Check if counts are on first line
        int tmp;
        if (sscanf(bufferp, "%s%d%d%d", header, &tmp, &nfaces, &nedges) == 4) {
          nverts = tmp;
        }
      }
      else {
        // Read counts from second line
        if ((sscanf(bufferp, "%d%d%d", &nverts, &nfaces, &nedges) != 3) || (nverts == 0)) {
          fprintf(stderr, "Syntax error reading header on line %d in file %s\n", line_count, filename);
          fclose(fp);
          return 0;
        }
      }
    }
    else if (vertex_count < nverts) {
      // Read vertex coordinates
      double x, y, z;
      if (sscanf(bufferp, "%lf%lf%lf", &x, &y, &z) != 3) {
        fprintf(stderr, "Syntax error with vertex coordinates on line %d in file %s\n", line_count, filename);
        fclose(fp);
        return 0;
      }

      // Create vertex
      CreateVertex(R3Point(x, y, z), R3zero_vector, R2zero_point);

      // Increment counter
      vertex_count++;
    }
    else if (face_count < nfaces) {
      // Read number of vertices in face 
      int face_nverts = 0;
      bufferp = strtok(bufferp, " \t");
      if (bufferp) face_nverts = atoi(bufferp);
      else {
        fprintf(stderr, "Syntax error with face on line %d in file %s\n", line_count, filename);
        fclose(fp);
        return 0;
      }

      // Read vertex indices for face
      R3MeshVertex *v1 = NULL;
      R3MeshVertex *v2 = NULL;
      R3MeshVertex *v3 = NULL;
      for (int i = 0; i < face_nverts; i++) {
        bufferp = strtok(NULL, " \t");
        if (bufferp) v3 = Vertex(atoi(bufferp));
        else {
          fprintf(stderr, "Syntax error with face on line %d in file %s\n", line_count, filename);
          fclose(fp);
          return 0;
        }

        // Create triangle
        if (v1) {
          if (!CreateFace(v1, v2, v3)) {
            // Must have been degeneracy (e.g., three faces sharing an edge), create new vertices
            R3MeshVertex *v1a = CreateVertex(v1->position, v1->normal, v1->texcoords);
            R3MeshVertex *v2a = CreateVertex(v2->position, v2->normal, v2->texcoords);
            R3MeshVertex *v3a = CreateVertex(v3->position, v3->normal, v3->texcoords);
            CreateFace(v1a, v2a, v3a);
          }
        }

        // Move to next triangle
        if (!v1) v1 = v2;
        v2 = v3;
      }

      // Increment counter
      face_count++;
    }
    else {
      // Should never get here
      fprintf(stderr, "Found extra text starting at line %d in file %s\n", line_count, filename);
      break;
    }
  }

  // Check whether read all vertices
  if ((vertex_count != nverts) || (NVertices() < nverts)) {
    fprintf(stderr, "Expected %d vertices, but read %d vertex lines and created %d vertices in file %s\n", 
      nverts, vertex_count, NVertices(), filename);
  }

  // Check whether read all faces
  if ((face_count != nfaces) || (NFaces() < nfaces)) {
    fprintf(stderr, "Expected %d faces, but read %d face lines and created %d faces in file %s\n", 
      nfaces, face_count, NFaces(), filename);
  }

  // Close file
  fclose(fp);

  // Return number of faces read
  return NFaces();
}



int R3Mesh::
WriteOff(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open file %s\n", filename);
    return 0;
  }

  // Write header
  fprintf(fp, "OFF\n");
  fprintf(fp, "%d %d %d\n", NVertices(), NFaces(), 0);

  // Write vertices
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    const R3Point& p = vertex->position;
    fprintf(fp, "%g %g %g\n", p.X(), p.Y(), p.Z());
    vertex->id = i;
  }

  // Write Faces
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *face = Face(i);
    fprintf(fp, "3");
    for (int j = 0; j < 3; j++) {
      fprintf(fp, " %d", face->vertices[j]->id);
    }
    fprintf(fp, "\n");
  }

  // Close file
  fclose(fp);

  // Return number of faces
  return NFaces();
}



////////////////////////////////////////////////////////////
// OBJ FILE INPUT/OUTPUT
////////////////////////////////////////////////////////////

int R3Mesh::
ReadObj(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Read body
  char buffer[1024];
  int line_count = 0;
  int triangle_count = 0;
  while (fgets(buffer, 1023, fp)) {
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) bufferp++;

    // Skip blank lines and comments
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;

    // Get keyword
    char keyword[80];
    if (sscanf(bufferp, "%s", keyword) != 1) {
      fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
      return 0;
    }

    // Check keyword
    if (!strcmp(keyword, "v")) {
      // Read vertex coordinates
      double x, y, z;
      if (sscanf(bufferp, "%s%lf%lf%lf", keyword, &x, &y, &z) != 4) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Create vertex
      CreateVertex(R3Point(x, y, z), R3zero_vector, R2zero_point);
    }
    else if (!strcmp(keyword, "f")) {
      // Read vertex indices
      int i1, i2, i3;
      if (sscanf(bufferp, "%s%d%d%d", keyword, &i1, &i2, &i3) != 4) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Get vertices
      R3MeshVertex *v1 = Vertex(i1-1);
      R3MeshVertex *v2 = Vertex(i2-1);
      R3MeshVertex *v3 = Vertex(i3-1);

      // Create face
      if (!CreateFace(v1, v2, v3)) {
        // Must have been degeneracy (e.g., three faces sharing an edge), create new vertices
        R3MeshVertex *v1a = CreateVertex(v1->position, v1->normal, v1->texcoords);
        R3MeshVertex *v2a = CreateVertex(v2->position, v2->normal, v2->texcoords);
        R3MeshVertex *v3a = CreateVertex(v3->position, v3->normal, v3->texcoords);
        CreateFace(v1a, v2a, v3a);
      }

      // Increment triangle counter
      triangle_count++;
    }
  }

  // Close file
  fclose(fp);

  // Return number of faces created
  return triangle_count;
}



int R3Mesh::
WriteObj(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "w"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Write vertices
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    const R3Point& p = vertex->position;
    fprintf(fp, "v %g %g %g\n", p.X(), p.Y(), p.Z());
    vertex->id = i;
  }

  // Write Faces
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *face = Face(i);
    R3MeshVertex *v0 = face->vertices[0];
    R3MeshVertex *v1 = face->vertices[1];
    R3MeshVertex *v2 = face->vertices[2];
    fprintf(fp, "f %d %d %d\n", v0->id + 1, v1->id + 1, v2->id + 1);
  }

  // Return number of faces written
  return NFaces();
}



////////////////////////////////////////////////////////////
// RAY FILE INPUT/OUTPUT
////////////////////////////////////////////////////////////

int R3Mesh::
ReadRay(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Read body
  char cmd[128];
  int triangle_count = 0;
  int command_number = 1;
  while (fscanf(fp, "%s", cmd) == 1) {
    if (!strcmp(cmd, "#vertex")) {
      // Read data
      double px, py, pz;
      double nx, ny, nz;
      double ts, tt;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf", &px, &py, &pz, &nx, &ny, &nz, &ts, &tt) != 8) {
        fprintf(stderr, "Unable to read vertex at command %d in file %s", command_number, filename);
        return 0;
      }

      // Create vertex
      R3Point point(px, py, pz);
      R3Vector normal(nx, ny, nz);
      R2Point texcoords(ts, tt);
      CreateVertex(point, normal, texcoords);
    }
    else if (!strcmp(cmd, "#shape_triangle")) {
      // Read data
      int m;
      int i1, i2, i3;
      if (fscanf(fp, "%d%d%d%d", &m, &i1, &i2, &i3) != 4) {
        fprintf(stderr, "Unable to read triangle at command %d in file %s", command_number, filename);
        return 0;
      }

      // Get vertices
      R3MeshVertex *v1 = Vertex(i1);
      R3MeshVertex *v2 = Vertex(i2);
      R3MeshVertex *v3 = Vertex(i3);

      // Create face
      if (!CreateFace(v1, v2, v3)) {
        // Must have been degeneracy (e.g., three faces sharing an edge), create new vertices
        R3MeshVertex *v1a = CreateVertex(v1->position, v1->normal, v1->texcoords);
        R3MeshVertex *v2a = CreateVertex(v2->position, v2->normal, v2->texcoords);
        R3MeshVertex *v3a = CreateVertex(v3->position, v3->normal, v3->texcoords);
        CreateFace(v1a, v2a, v3a);
      }

      // Increment triangle counter
      triangle_count++;
    }
	
    // Increment command number
    command_number++;
  }

  // Close file
  fclose(fp);

  // Return number of faces created
  return triangle_count;
}



int R3Mesh::
WriteRay(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "w"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Write vertices
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    const R3Point& p = vertex->position;
    const R3Vector& n = vertex->normal;
    const R2Point& t = vertex->texcoords;
    fprintf(fp, "#vertex %g %g %g  %g %g %g  %g %g\n", p.X(), p.Y(), p.Z(), 
      n.X(), n.Y(), n.Z(), t.X(), t.Y());
    vertex->id = i;
  }

  // Write Faces
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *face = Face(i);
    R3MeshVertex *v0 = face->vertices[0];
    R3MeshVertex *v1 = face->vertices[1];
    R3MeshVertex *v2 = face->vertices[2];
    fprintf(fp, "#shape_triangle 0 %d %d %d\n", v0->id, v1->id, v2->id);
  }

  // Close file
  fclose(fp);

  // Return number of faces written
  return NFaces();
}



////////////////////////////////////////////////////////////
// MESH VERTEX MEMBER FUNCTIONS
////////////////////////////////////////////////////////////

R3MeshVertex::
R3MeshVertex(void)
  : position(0, 0, 0),
    normal(0, 0, 0),
    texcoords(0, 0),
    id(0)
{
}



R3MeshVertex::
R3MeshVertex(const R3MeshVertex& vertex)
  : position(vertex.position),
    normal(vertex.normal),
    texcoords(vertex.texcoords),
    id(0)
{
}



R3MeshVertex::
R3MeshVertex(const R3Point& position, const R3Vector& normal, const R2Point& texcoords)
  : position(position),
    normal(normal),
    texcoords(texcoords),
    id(0)
{
}



////////////////////////////////////////////////////////////
// MESH FACE MEMBER FUNCTIONS
////////////////////////////////////////////////////////////

R3MeshFace::
R3MeshFace(void)
  : plane(0, 0, 0, 0),
    id(0)
{
  // Initialize vertices
  vertices[0] = NULL;
  vertices[1] = NULL;
  vertices[2] = NULL;
}



R3MeshFace::
R3MeshFace(const R3MeshFace& face)
  : plane(face.plane),
    id(0)
{
  // Initialize vertices
  vertices[0] = face.vertices[0];
  vertices[1] = face.vertices[1];
  vertices[2] = face.vertices[2];
}



R3MeshFace::
R3MeshFace(R3MeshVertex *v1, R3MeshVertex *v2, R3MeshVertex *v3)
  : plane(0, 0, 0, 0),
    id(0)
{
  // Initialize vertices
  vertices[0] = v1;
  vertices[1] = v2;
  vertices[2] = v3;

  // Initialize plane
  plane = R3Plane(v1->position, v2->position, v3->position);
}



void R3Mesh::
Draw(void) const
{
  // Draw mesh
  glBegin(GL_TRIANGLES);
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *face = Face(i);
    const R3Vector& normal = face->plane.Normal();
    glNormal3d(normal[0], normal[1], normal[2]);
    for (int j = 0; j < 3; j++) {
      R3MeshVertex *vertex = face->vertices[j];
      const R3Point& p = vertex->position;
      // const R3Vector& n = vertex->normal;
      const R2Point& t = vertex->texcoords;
      glTexCoord2d(t[0], t[1]);
      // glNormal3d(n[0], n[1], n[2]);
      glVertex3d(p[0], p[1], p[2]);
    }
  }
  glEnd();
}



void R3Mesh::
Outline(void) const
{
  // Draw mesh in wireframe
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  Draw();
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}



void R3Mesh::
Print(FILE *fp) const
{
  // Print 
  fprintf(fp, "Mesh");
}


