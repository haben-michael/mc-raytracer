// Include file for ray tracing code

R2Image *RenderImage(R3Scene *scene, 
  int width, int height, 
  int max_reflections, double min_luminance, int n_samples);
