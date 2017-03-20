

// Include files
#ifdef _WIN32
#include <windows.h>
#endif

#include "R2/R2.h"
#include "R3/R3.h"
#include "R3Scene.h"
#include "raytrace.h"



// Program arguments

static char options[] =
"  -help\n"
"  -width <int:width>\n"
"  -height <int:height>\n"
"  -max_reflections <int:nbounces>\n"
"  -min_luminance <real:value>\n";



static void 
ShowUsage(void)
{
  // Print usage message and exit
  fprintf(stderr,  "Usage: raypro input_scene output_image [  -option [arg ...] ...]\n");
  fprintf(stderr, options);
  exit(EXIT_FAILURE);
}



static void 
CheckOption(char *option, int argc, int minargc)
{
  // Check if there are enough remaining arguments for option
  if (argc < minargc)  {
    fprintf(stderr, "Too few arguments for %s\n", option);
    ShowUsage();
    exit(-1);
  }
}


static R3Scene *
ReadScene(const char *filename, int width, int height)
{
  // Allocate scene
  R3Scene *scene = new R3Scene();
  if (!scene) {
    fprintf(stderr, "Unable to allocate scene\n");
    return NULL;
  }

  // Read scene
  if (!scene->Read(filename)) {
    fprintf(stderr, "Unable to read scene from %s\n", filename);
    return NULL;
  }

  // Provide default camera
  if (scene->camera.xfov == 0) {
    double scene_radius = scene->BBox().DiagonalRadius();
    R3Point scene_center = scene->BBox().Centroid();
    scene->camera.towards = R3Vector(0, 0, -1);
    scene->camera.up = R3Vector(0, 1, 0);
    scene->camera.right = R3Vector(1, 0, 0);
    scene->camera.eye = scene_center - 3 * scene_radius * scene->camera.towards;
    scene->camera.xfov = 0.5;
    scene->camera.yfov = 0.5;
    scene->camera.neardist = 0.01 * scene_radius;
    scene->camera.fardist = 100 * scene_radius;
  }

  // Provide default light
  if (scene->NLights() == 0) {
    // Create first directional light
    R3Light *light = new R3Light();
    light->type = R3_DIRECTIONAL_LIGHT;
    light->color = R3Rgb(1,1,1,1);
    light->position = R3Point(0, 0, 0);
    light->direction = R3Vector(-3,-4,-5);
    light->constant_attenuation = 0;
    light->linear_attenuation = 0;
    light->quadratic_attenuation = 0;
    light->angle_attenuation = 0;
    light->angle_cutoff = M_PI;
    scene->lights.push_back(light);

    // Create second directional light
    light = new R3Light();
    light->type = R3_DIRECTIONAL_LIGHT;
    light->color = R3Rgb(0.5, 0.5, 0.5, 1);
    light->position = R3Point(0, 0, 0);
    light->direction = R3Vector(3,2,3);
    light->constant_attenuation = 0;
    light->linear_attenuation = 0;
    light->quadratic_attenuation = 0;
    light->angle_attenuation = 0;
    light->angle_cutoff = M_PI;
    scene->lights.push_back(light);
  }

  // Adjust camera vertical field of view to match aspect ratio of image
  scene->camera.yfov = atan(tan(scene->camera.xfov) * (double) height / (double) width); 

  // Return scene
  return scene;
}



int 
main(int argc, char **argv)
{
  // Look for help
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "-help")) {
      ShowUsage();
    }
  }

  // Read input and output filenames
  if (argc < 3)  ShowUsage();
  argv++, argc--; // First argument is program name
  char *input_scene_name = *argv; argv++, argc--; 
  char *output_image_name = *argv; argv++, argc--; 

  // Initialize arguments to default values
  int width = 256;
  int height = 256;
  int max_reflections = 3;
  double min_luminance = 0.01;
  int n_samples = 1;
  R3Vector motion_vec(0,0,0);
  int sampling_filter = -1;
  int accel_method = 0;

  // Parse arguments 
  while (argc > 0) {
    if (!strcmp(*argv, "-width")) {
      CheckOption(*argv, argc, 2);
      width = atoi(argv[1]);
      argv += 2, argc -= 2;
    }
    else if (!strcmp(*argv, "-height")) {
      CheckOption(*argv, argc, 2);
      height = atoi(argv[1]);
      argv += 2, argc -= 2;
    }
    else if (!strcmp(*argv, "-maxreflections")) {
      CheckOption(*argv, argc, 2);
      max_reflections = atoi(argv[1]);
      argv += 2, argc -= 2;
    }
    else if (!strcmp(*argv, "-minluminance")) {
      CheckOption(*argv, argc, 2);
      min_luminance = atof(argv[1]);
      argv += 2, argc -= 2;
    }
    else if (!strcmp(*argv, "-samples")) {
      CheckOption(*argv, argc, 2);
      n_samples = atof(argv[1]);
      argv += 2, argc -= 2;
    }
    else if (!strcmp(*argv, "-motion")) {
      CheckOption(*argv, argc, 4);
	motion_vec = R3Vector(atof(argv[1]),atof(argv[2]),atof(argv[3]));
      argv += 4, argc -= 4;
    }
    else if (!strcmp(*argv, "-filter")) {
      CheckOption(*argv, argc, 2);
      sampling_filter = atof(argv[1]);
      argv += 2, argc -= 2;
    }
    else if (!strcmp(*argv, "-accel")) {
      CheckOption(*argv, argc, 2);
      accel_method = atof(argv[1]);
      argv += 2, argc -= 2;
    }

    else {
      // Unrecognized program argument
      fprintf(stderr,  "meshpro: invalid option: %s\n", *argv);
      ShowUsage();
    }
  }

  // Read scene
  R3Scene *scene = ReadScene(input_scene_name, width, height);
  if (!scene) {
    fprintf(stderr, "Unable to read scene from %s\n", input_scene_name);
    exit(-1);
  }

  // Render image
  R2Image *image = RenderImage(scene, width, height, max_reflections, min_luminance, n_samples);//, motion_vec, sampling_filter, accel_method);
  if (!image) {
    fprintf(stderr, "Unable to render image from scene\n");
    exit(-1);
  }

  // Write output image
  if (!image->Write(output_image_name)) {
    fprintf(stderr, "Unable to write image to %s\n", output_image_name);
    exit(-1);
  }

  // Delete everything
  delete scene;
  delete image;

  // Return success
  return EXIT_SUCCESS;
}



