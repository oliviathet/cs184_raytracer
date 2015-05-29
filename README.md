# cs184_raytracer
UC Berkeley CS 184 (Computer Graphics) Raytracer partner project

Rendered scenes using Phong Illumination Model, linear transformations, reflections, multiple lights, transparency with refraction, anti-aliasing, lens effects, texture mapping.

To compile:
make

To run (example):
./as1 -ka 0.05 0.05 0.05 -kd 1 1 1 -ks 1 1 1 -sp 64 -pl 2 2 2 0.3 0.3 0.8 -dl 0 -1 -0.6 0.4 0.1 0.5 -ms 2 -ts 5

Flags:
ambient term: -ka r g b
diffuse term: -kd r g b
specular term: -ks r g b
specular exponent: -sp n
point light: -pl x y z r g b
directional light: -dl x y z r g b
multiple spheres (default setting is 1 sphere): -ms n
toon shading (default setting is all 8 bits): -ts n
