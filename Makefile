#
# Linux makefile for assignment #2
#



#
# List of source files for each of two applications
#

RAYTRACE_SRCS=raypro.cpp raytrace.cpp R3Scene.cpp
RAYVIEW_SRCS=rayview.cpp raytrace.cpp R3Scene.cpp



#
# List of source files for each of two applications
#

RAYTRACE_OBJS=$(RAYTRACE_SRCS:.cpp=.o)
RAYVIEW_OBJS=$(RAYVIEW_SRCS:.cpp=.o)



#
# Compile and link options
#

CC=g++
CPPFLAGS=-Wall -I. -g -DUSE_JPEG
LDFLAGS=-g



#
# Libraries
#

LIBS=R3/libR3.a R2/libR2.a jpeg/libjpeg.a 
#GLLIBS=-framework opengl -framework GLUT 
GLLIBS=-lglut -lGLU -lGL 



#
# Compile command
#

%.o: %.cpp
	    $(CC) $(CPPFLAGS) -c $< -o $@



#
# Make targets
#

all: $(LIBS) raytrace rayview

raytrace: $(LIBS) $(RAYTRACE_OBJS)
	    $(CC) -o raytrace $(CPPFLAGS) $(LDFLAGS) $(RAYTRACE_OBJS) $(LIBS) $(GLLIBS) -lm

rayview: $(LIBS) $(RAYVIEW_OBJS)
	    $(CC) -o rayview $(CPPFLAGS) $(LDFLAGS) $(RAYVIEW_OBJS) $(LIBS) $(GLLIBS) -lm

R3/libR3.a: 
	    cd R3; make

R2/libR2.a: 
	    cd R2; make

jpeg/libjpeg.a: 
	    cd jpeg; make

clean:
	    -  rm -f *~ *.o */*.o */*/*.o raytrace rayview $(LIBS)

