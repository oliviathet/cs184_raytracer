CC = g++
ifeq ($(shell sw_vers 2>/dev/null | grep Mac | awk '{ print $$2}'),Mac)
	CFLAGS = -g -DGL_GLEXT_PROTOTYPES -I./include/ -I/usr/X11/include -DOSX
	LDFLAGS = -framework GLUT -framework OpenGL \
    	-L"/System/Library/Frameworks/OpenGL.framework/Libraries" \
    	-lGL -lGLU -lm -lstdc++
else
	CFLAGS = -g -DGL_GLEXT_PROTOTYPES -Iglut-3.7.6-bin
	LDFLAGS = -lglut -lGLU
	FLAGS += -O3
	FLAGS += -std=c++11
	FLAGS += -D_DEBUG -g Wall
endif

all: scene 
scene: scene.o 
	$(CC) $(CFLAGS) -o scene scene.o lodepng.o $(LDFLAGS)	

scene.o: Scene.cpp
	$(CC) $(CFLAGS) -c Scene.cpp -o scene.o
	
lodepng.o: lodepng.cpp lodepng.h
	$(CC) $(CFLAGS) -c lodepng.cpp -o lodepng.o

	g++ -I Eigen/ Scene.cpp -o scene 

clean: 
	$(RM) *.o scene
 


