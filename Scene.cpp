#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string>
#include <limits>
#include <cfloat>
#include <cstdlib>
#include "Eigen/Geometry"

#include "Vector3.h"
#include "Vector4.h"
#include "Sample.h"
#include "Point.h"
#include "Normal.h"
#include "Matrix4.h"
#include "Transformation.h"
#include "Ray.h"
#include "DifferentialGeometry.h"
#include "Color.h"
#include "BRDFCoefficients.h"
#include "Material.h"
#include "Shape.h"
#include "Primitive.h"
#include "Intersection.h"
#include "AggregatePrimitive.h"
#include "Film.h"
#include "Camera.h"
#include "Light.h"

#include "lodepng.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <bitset>
#include <algorithm>

using namespace std;


/*
 * Scene file has:
 *   transformation1
 *   transformation2
 *   sph 0 0 0 2
 *
 *   i.e. transformation 1 and 2 transform the SPHERE INTO THE ELLIPSOID
 *
 * (1) Store an ELLIPSOID object that has scaledX, scaledY, scaledZ, translateX, translateY, translateZ, rotationObject;
 *     populate these attributes based on the transformations given.
 *
 *     Also store a vector<transformations> = [transformation1, transformation2] = M   to indicate how the sphere was transformed
 *     into the ellipsoid
 *
 * (2) In the intersect() method, transform the original input ray in world coordinates (i.e. ellipsoid coordinates)
 *     with the INVERSE transformations (M^ -1) to send the ray into OBJECT coordinates (i.e. sphere coordinates)
 *
 * (3) In the intersect() method, call the ellipsoid's PARENT's (i.e. it's corresponding sphere's) intersection method
 *     with the TRANSFORMED ray to find the intersection point of the TRANSFORMED ray and the ORIGINAL SPHERE.
 *     This intersection point is in OBJECT coordinates.
 *
 *     (NOTE: do this with super.
 *
 * (4) Perform the original stored transformation (i.e. M) on the intersection point to transform the intersection point
 *     back into WORLD (i.e. ellipsoid) coordinates.
 *
 * (5) Done.
 *
 */

/*

GeometricPrimitive
	Members:
		Transformation objToWorld, worldToObj;
		Shape* shape;
		Material* mat;

	Methods:
		bool intersect(Ray& ray, float* thit, Intersection* in)  {
			// Set worldToObj and objToWorld as identity by default
			Ray oray = worldToObj*ray;
			LocalGeo olocal;

			// If 'shape' is of type ellipsoid, then simply call its parent routine
			// with Sphere::intersect()

			// If 'shape' was of type sphere, then it will just operate as normal,
			// and oray will be the original ray anyway
			if (!shape->intersect(oray, thit, &olocal))  return false;
			in->primitive = this;
			in->local = objToWorld*olocal;
			return true;
		}

		bool intersectP(Ray& ray) {
			Ray oray = worldToObj*ray;
			return shape->intersectP(oray);
		}

		void getBRDF(LocalGeo& local, BRDF* brdf) {
			material->getBRDF(local, brdf);
		}

 */



//****************************************************
// Forward Declarations
//****************************************************
static void printRay(Ray ray);
static void printSample(Sample sample);
static void printColor(Color color);
static void printPoint(Point point);
void printCommandLineOptionVariables();
void printSamples();
Color applyShadingModel(DifferentialGeometry differentialGeometry, BRDFCoefficients brdf, Ray lightRay, Color lightColor, Point raySource);
class RayTracer;


//****************************************************
// Global Variables
//****************************************************
bool debug;
std::vector<Sample> samples;
const char * filename;
Film film;
Camera camera;
int recursionDepth;
AggregatePrimitive aggregatePrimitive;
std::vector<PointLight> point_lights;
std::vector<DirectionalLight> directional_lights;
std::vector<AmbientLight> ambient_lights;
std::vector<Point> objFileVertices;
Transformation currentlySeenTransformation;
// Most recently seen material
float currentKar, currentKag, currentKab, currentKdr, currentKdg, currentKdb, currentKsr, currentKsg, currentKsb, currentKsp, currentKrr, currentKrg, currentKrb, currentRi;
int currentSp;

// Extra Credit global variables
bool softShadows;
bool distributedRayTracing;
bool checkerboard;



class RayTracer {
public:
	void trace(Ray& ray, int depth, Color* color) {

		if (depth > recursionDepth) {
			color->r = 0;
			color->g = 0;
			color->b = 0;
			return;
		}

		//			// For testing only: shade coordinate red if ray intersects it, else shade black
		//			if (!aggregatePrimitive.intersectP(ray)) {
		//				color->r = 0;
		//				color->g = 0;
		//				color->b = 0;
		//				return;
		//			} else {
		//				//For testing purposes, simply shade red if the ray intersects the point
		//				color->r = 255;
		//				color->g = 0;
		//				color->b = 0;
		//			}

		// TODO: this may be wrong
		float tHit;
		Intersection intersection;

		// This method will populate tHit and intersection if there is an intersection with this ray and any primitive.
		if (!aggregatePrimitive.intersect(ray, &tHit, &intersection) ) {
			// If no intersection, then make the color black and return
			color->r = 0;
			color->g = 0;
			color->b = 0;
			return;
		}

		BRDFCoefficients brdf;
		// This method will populate the brdf variable with the brdf values of the intersection primitive.
		brdf = intersection.primitive->getBRDF(intersection.differentialGeometry, &brdf);

		// Initialize a new Color with R = G = B = 0.0
		// This color will be appended to with our shading model
		Color colorOfPixel;
		Ray lightRay;
		Color lightColor;

		float intersectionXCoor = intersection.differentialGeometry.position.x;
		float intersectionYCoor = intersection.differentialGeometry.position.y;
		float intersectionZCoor = intersection.differentialGeometry.position.z;

		if (ambient_lights.size() > 0) {
			// We want to add the ambient term even if our shape is blocked by a light
			for (std::vector<AmbientLight>::size_type i = 0; i < ambient_lights.size(); i++) {

				// If our checkerboard global boolean is on, then change ambient term into a checkerboard
				if (checkerboard) {
					bool x = (int)((intersectionXCoor * 1000 + 3893343)/50) % 2 == 1;
					bool y = (int)((intersectionYCoor * 1000 + 3893343)/50) % 2 == 1;
					bool z = (int)((intersectionZCoor * 1000 + 3893343)/50) % 2 == 1;

					if (x xor y xor z) {
						color->r += 2 * brdf.ka.r * ambient_lights[i].r;
						color->g += 2 * brdf.ka.g * ambient_lights[i].g;
						color->b += 2 * brdf.ka.b * ambient_lights[i].b;
					}
					// No checkerboard
				} else {
					color->r += (brdf.ka.r * ambient_lights[i].r);
					color->g += (brdf.ka.g * ambient_lights[i].g);
					color->b += (brdf.ka.b * ambient_lights[i].b);
				}
			}

		} else {
			// If our checkerboard global boolean is on, then change ambient term into a checkerboard
			if (checkerboard) {
				bool x = (int)((intersectionXCoor * 1000 + 3893343)/50) % 2 == 1;
				bool y = (int)((intersectionYCoor * 1000 + 3893343)/50) % 2 == 1;
				bool z = (int)((intersectionZCoor * 1000 + 3893343)/50) % 2 == 1;

				if (x xor y xor z) {
					color->r += 2 * brdf.ka.r;
					color->g += 2 * brdf.ka.g;
					color->b += 2 * brdf.ka.b;
				}
				// No checkerboard
			} else {
				color->r += brdf.ka.r;
				color->g += brdf.ka.g;
				color->b += brdf.ka.b;
			}
		}

		// There is an intersection, so we have to loop through all the light source
		// and consider their contributions to the intersection pixel
		for (std::vector<DirectionalLight>::size_type i = 0; i < directional_lights.size(); i++) {

			// Generate light ray from the intersection point to the light position
			// For directional lights, this is generated by subtracting the position of the intersection point
			// from the position of the light (IN WORLD COORDINATES)
			directional_lights[i].generateLightRay(intersection.differentialGeometry, &lightRay, &lightColor);

			Ray reversedLightRay = Ray(lightRay.position, lightRay.direction.scaleVector(-1.0), lightRay.t_min, lightRay.t_max);

			// ***** Regular shading + shadows *****
			// If the light ray is not blocked, we apply our shading model
			if (!aggregatePrimitive.intersectP(reversedLightRay)) {
				Color colorToAdd = applyShadingModel(
						intersection.differentialGeometry,
						brdf,
						lightRay,
						Color(directional_lights[i].r, directional_lights[i].g, directional_lights[i].b),
						ray.position);

				color->r += colorToAdd.r;
				color->g += colorToAdd.g;
				color->b += colorToAdd.b;
			}
		}

		for (std::vector<PointLight>::size_type i = 0; i < point_lights.size(); i++) {
			point_lights[i].generateLightRay(intersection.differentialGeometry, &lightRay, &lightColor);
			if (!aggregatePrimitive.intersectP(lightRay)) {

				// For now, we just ignore shadows and reflections and just apply our shading model
				// i.e. just call:
				Color colorToAdd = applyShadingModel(
						intersection.differentialGeometry,
						brdf, lightRay,
						Color(point_lights[i].r, point_lights[i].g, point_lights[i].b),
						ray.position);

				color->r += colorToAdd.r;
				color->g += colorToAdd.g;
				color->b += colorToAdd.b;
			}
		}

		// ***** Reflection but no refraction *****
		if (brdf.kr.greaterThanZero() && !brdf.isRefractive()) {
			// First, create our reflection ray...
			// r = d - 2(d . n)n

			// n
			Vector3 directional_normal_vector =
					Vector3::normalizeVector(Vector3(intersection.differentialGeometry.normal.x, intersection.differentialGeometry.normal.y, intersection.differentialGeometry.normal.z));

			Point intersectionMinusRaySource = intersection.differentialGeometry.position.subtractPoint(ray.position);

			// d (view vector) = intersection coordinate - source of ray
			// (we can't just use the intersection point, because our view is relative from the LAST HIT primitive)
			Vector3 view_vector = Vector3::normalizeVector(
					Vector3(intersectionMinusRaySource.x, intersectionMinusRaySource.y, intersectionMinusRaySource.z));


			float dDotN = view_vector.dotProduct(directional_normal_vector);
			Vector3 twoTimesdDotNTimesN = directional_normal_vector.scaleVector(dDotN * 2.0);
			Vector3 directional_reflective_vector_prenormalize = view_vector.subtractVector(twoTimesdDotNTimesN);
			// r
			Vector3 directional_reflective_vector = Vector3::normalizeVector(directional_reflective_vector_prenormalize);

			// NOTE: Not sure if this is right...
			// Offset our reflection ray's position a tiny bit in the direction of the normal
			Vector3 normalDirectionError = Vector3(intersection.differentialGeometry.normal.x, intersection.differentialGeometry.normal.y, intersection.differentialGeometry.normal.z).scaleVector(0.000001);
			Point offsetRayPosition = Point(intersection.differentialGeometry.position.x + normalDirectionError.x,
					intersection.differentialGeometry.position.y + normalDirectionError.y,
					intersection.differentialGeometry.position.z + normalDirectionError.z);

			Ray reflectionRay = Ray(offsetRayPosition, directional_reflective_vector, 0.001, FLT_MAX);

			Color tempColor;

			// Recursively call trace()
			trace(reflectionRay, depth + 1, &tempColor);

			color->r += brdf.kr.r * tempColor.r;
			color->g += brdf.kr.g * tempColor.g;
			color->b += brdf.kr.b * tempColor.b;

		// Refractive and reflective
		} else if (brdf.isRefractive() && brdf.kr.greaterThanZero()) {

			// n
			Vector3 directional_normal_vector =
					Vector3::normalizeVector(Vector3(intersection.differentialGeometry.normal.x, intersection.differentialGeometry.normal.y, intersection.differentialGeometry.normal.z));

			Point intersectionMinusRaySource = intersection.differentialGeometry.position.subtractPoint(ray.position);

			// d (view vector) = intersection coordinate - source of ray
			// (we can't just use the intersection point, because our view is relative from the LAST HIT primitive)
			Vector3 view_vector = Vector3::normalizeVector(
					Vector3(intersectionMinusRaySource.x, intersectionMinusRaySource.y, intersectionMinusRaySource.z));


			float dDotN = view_vector.dotProduct(directional_normal_vector);
			Vector3 twoTimesdDotNTimesN = directional_normal_vector.scaleVector(dDotN * 2);
			Vector3 directional_reflective_vector_prenormalize = view_vector.subtractVector(twoTimesdDotNTimesN);
			// r
			Vector3 directional_reflective_vector = Vector3::normalizeVector(directional_reflective_vector_prenormalize);

			// NOTE: Not sure if this is right...
			// Offset our reflection ray's position a tiny bit in the direction of the normal
			Vector3 normalDirectionError = Vector3(intersection.differentialGeometry.normal.x, intersection.differentialGeometry.normal.y, intersection.differentialGeometry.normal.z).scaleVector(0.00000001);
			Point offsetRayPosition = Point(intersection.differentialGeometry.position.x + normalDirectionError.x,
					intersection.differentialGeometry.position.y + normalDirectionError.y,
					intersection.differentialGeometry.position.z + normalDirectionError.z);

			Ray reflectionRay = Ray(offsetRayPosition, directional_reflective_vector, 0.0001, FLT_MAX);

			Color tempColor;

			// Recursively call trace()
			trace(reflectionRay, depth + 1, &tempColor);

			color->r += brdf.kr.r * tempColor.r;
			color->g += brdf.kr.g * tempColor.g;
			color->b += brdf.kr.b * tempColor.b;

			Vector3 directional_refractive_vector;
			float n1 = 1.0;
			float n2 = brdf.ri;
			float c = -1.0 * dDotN;
			float c1 = c;
			float nc, c2;
			// Entering sphere
			if (-c < 0) {
				nc = n1 / n2;
				c2 = sqrt(1 - nc * nc * (1- c1 * c1));
				Vector3 refractive_vector_term1 = directional_normal_vector.scaleVector(nc * c1 - c2);
				Vector3 refractive_vector_term2 = view_vector.scaleVector(nc);
				directional_refractive_vector = refractive_vector_term1.addVector(refractive_vector_term2);
			} else {
				// Leaving sphere
				nc = n2 / n1;
				c2 = sqrt(1 - nc * nc * (1- c1 * c1));
				Vector3 refractive_vector_term1 = directional_normal_vector.scaleVector(nc * c1 - c2);
				Vector3 refractive_vector_term2 = view_vector.scaleVector(nc);
				directional_refractive_vector = refractive_vector_term1.addVector(refractive_vector_term2);
				c = directional_refractive_vector.dotProduct(directional_normal_vector);
			}
			float R0 = pow((n2 - 1)/(n2 + 1), 2);
			float R = R0 + (1.0 - R0) * pow(1.0 - c, 5.0);

			if ((1.0 - R) > 0.0) {
				Ray refractionRay = Ray(intersection.differentialGeometry.position, directional_refractive_vector, 0.0001, FLT_MAX);

				Color tempColor2;

				// Recursively call trace()
				trace(refractionRay, depth + 1, &tempColor2);

				color->r += R * brdf.kr.r + (1.0 - R) * tempColor2.r;
				color->g += R * brdf.kr.g + (1.0 - R) * tempColor2.g;
				color->b += R * brdf.kr.b + (1.0 - R) * tempColor2.b;

			// Not reflective or refractive
			} else {
				return;
			}
		}
	}
};

// (Declaration of RayTracer global variable)
RayTracer rayTracer;


//****************************************************
// applies Phong shading model to differentialGeometry.position
//****************************************************
Color applyShadingModel(DifferentialGeometry differentialGeometry, BRDFCoefficients brdf, Ray lightRay, Color lightColor, Point raySource) {

	// ***** BEGIN COMPUTATION OF PHONG SHADING MODEL ***** //
	// NOTE: Ambient term (ka) is appended outisde of this model

	float resultant_rgb_sum_of_pixel_r = 0;
	float resultant_rgb_sum_of_pixel_g = 0;
	float resultant_rgb_sum_of_pixel_b = 0;

	// Set viewer vector to -1 * incoming_ray's_vector
	Vector3 viewer_vector = Vector3::normalizeVector(Vector3(-1.0 * differentialGeometry.position.x, -1.0 * differentialGeometry.position.y, -1.0 * differentialGeometry.position.z));

	// **************************************
	// For directional light
	// **************************************
	if (lightRay.t_max == FLT_MAX) {

		// Direction of light ray computed in 'generateLightRay()'
		// Pointing from "the light" (it exists at infinity though) TOWARDS the surface
		Vector3 prenormalized_directional_light_vector = lightRay.direction;

		// Change orientation of light vector to point outwards to "the light"
		Vector3 directional_light_vector = Vector3::normalizeVector(prenormalized_directional_light_vector.scaleVector(-1));

		// NOTE: (x, y, z) is in world coordinates now, not relative to center of sphere
		// NOTE: we should defer this logic to Sphere/Triangle, as follows:
		Vector3 directional_normal_vector = Vector3(differentialGeometry.normal.x, differentialGeometry.normal.y, differentialGeometry.normal.z);

		float directional_diffuse_dot_product = fmax(directional_light_vector.dotProduct(directional_normal_vector), 0);

		float directional_diffuse_r = brdf.kd.r * lightColor.r * directional_diffuse_dot_product;
		float directional_diffuse_g = brdf.kd.g * lightColor.g * directional_diffuse_dot_product;
		float directional_diffuse_b = brdf.kd.b * lightColor.b * directional_diffuse_dot_product;

		// Calculate specular term
		Vector3 directional_reflective_vector = directional_normal_vector.scaleVector(directional_light_vector.dotProduct(directional_normal_vector) * 2).subtractVector(directional_light_vector);

		float directional_specular_dot_product_term = pow(fmax(directional_reflective_vector.dotProduct(viewer_vector), 0), brdf.sp);

		float directional_specular_r = brdf.ks.r * lightColor.r * directional_specular_dot_product_term;
		float directional_specular_g = brdf.ks.g * lightColor.g * directional_specular_dot_product_term;
		float directional_specular_b = brdf.ks.b * lightColor.b * directional_specular_dot_product_term;

		// Combine three contributions together
		resultant_rgb_sum_of_pixel_r += (directional_diffuse_r + directional_specular_r);
		resultant_rgb_sum_of_pixel_g += (directional_diffuse_g + directional_specular_g);
		resultant_rgb_sum_of_pixel_b += (directional_diffuse_b + directional_specular_b);
	} else {
		// **************************************
		// For point light
		// **************************************

		// Calculate diffuse term

		// Location of point light given by command line options (i.e. x, y, z)

		Vector3 point_normal_vector = Vector3::normalizeVector(Vector3(differentialGeometry.normal.x, differentialGeometry.normal.y, differentialGeometry.normal.z));

		Vector3 prenormalized_point_light_vector = lightRay.direction;
		Vector3 point_light_vector = Vector3::normalizeVector(prenormalized_point_light_vector);

		float point_diffuse_dot_product = fmax(point_light_vector.dotProduct(point_normal_vector), 0);
		float point_diffuse_r = brdf.kd.r * lightColor.r * point_diffuse_dot_product;
		float point_diffuse_g = brdf.kd.g * lightColor.g * point_diffuse_dot_product;
		float point_diffuse_b = brdf.kd.b * lightColor.b * point_diffuse_dot_product;

		// Calculate specular term
		Vector3 point_reflective_vector = point_normal_vector.scaleVector(point_light_vector.dotProduct(point_normal_vector) * 2).subtractVector(point_light_vector);
		float point_specular_dot_product_term = pow(fmax(point_reflective_vector.dotProduct(viewer_vector), 0), brdf.sp);
		float point_specular_r = brdf.ks.r * lightColor.r * point_specular_dot_product_term;
		float point_specular_g = brdf.ks.g * lightColor.g * point_specular_dot_product_term;
		float point_specular_b = brdf.ks.b * lightColor.b * point_specular_dot_product_term;

		// Combine three contributions together
		resultant_rgb_sum_of_pixel_r += (point_diffuse_r + point_specular_r);
		resultant_rgb_sum_of_pixel_g += (point_diffuse_g + point_specular_g);
		resultant_rgb_sum_of_pixel_b += (point_diffuse_b + point_specular_b);
	}

	return Color(resultant_rgb_sum_of_pixel_r, resultant_rgb_sum_of_pixel_g, resultant_rgb_sum_of_pixel_b);
}


//****************************************************
// Debug printing functions
//****************************************************

static void printRay(Ray ray) {
	if (debug) {
		printf("Ray: (%f, %f, %f) + t(%f, %f, %f)\n", ray.position.x, ray.position.y, ray.position.z, ray.direction.x, ray.direction.y, ray.direction.z);
	}
}

static void printSample(Sample sample) {
	if (debug) {
		std::cout << "Sample: x = " << sample.x << "; y = " << sample.y << "\n";
	}
}

static void printVector(Vector3 vector) {
	if (debug) {
		std::cout << "Vector: x = " << vector.x << "; y = " << vector.y << "; z = " << vector.z << "\n";
	}
}

static void printColor(Color color) {
	if (debug) {
		printf("Color: (r, g, b) = (%f, %f, %f)\n", color.r, color.g, color.b);
	}
}

static void printPoint(Point point) {
	if (debug) {
		printf("Point: (x, y, z) = (%f, %f, %f)\n", point.x, point.y, point.z);
	}
}


void printGlobalVariables()
{
	if (debug)
	{
		std::cout << "\n***** BEGIN PRINTING GLOBAL VARIABLES *****\n";
		std::cout << "  " << "Film Width: " << film.width << "\n";
		std::cout << "  " << "Film Height: " << film.height << "\n\n";
		std::cout << "  Recursion Depth: " << recursionDepth << "\n\n";

		std::cout << "  Camera:\n";
		std::cout << "    Eye: ";
		printPoint(camera.eye);
		std::cout << "    Bottom Left Corner: ";
		printPoint(camera.imagePlaneBottomLeft);
		std::cout << "    Top Left Corner: ";
		printPoint(camera.imagePlaneTopLeft);
		std::cout << "    Bottom Right: ";
		printPoint(camera.imagePlaneBottomRight);
		std::cout << "    Top Right: ";
		printPoint(camera.imagePlaneTopRight);
		std::cout << "\n";

		std::cout << "Directional Lights:\n";
		if (directional_lights.size() == 0)
		{
			std::cout << " (none)\n";
		}
		for (std::vector<DirectionalLight>::size_type i = 0; i < directional_lights.size(); i++)
		{
			std::cout << "  " << "Light " << (i + 1) << "\n";
			std::cout << "     " << "x: " << directional_lights[i].x << " y: " << directional_lights[i].x << " z: " << directional_lights[i].x << "\n";
			std::cout << "     " << "r: " << directional_lights[i].r << " g: " << directional_lights[i].g << " b: " << directional_lights[i].b << "\n";
		}

		std::cout << "\nPoint Lights:\n";
		if (point_lights.size() == 0)
		{
			std::cout << " (none)\n";
		}
		for (std::vector<PointLight>::size_type i = 0; i < point_lights.size(); i++)
		{
			std::cout << "  " << "Light " << (i + 1) << "\n";
			std::cout << "     " << "x: " << point_lights[i].x << " y: " << point_lights[i].x << " z: " << point_lights[i].x << "\n";
			std::cout << "     " << "r: " << point_lights[i].r << " g: " << point_lights[i].g << " b: " << point_lights[i].b << "\n";
		}

		std::cout << "\n";

		std::cout << "Number of primitives: " << aggregatePrimitive.listOfPrimitives.size() << "\n\n";
		for (std::vector<GeometricPrimitive*>::size_type i = 0; i < aggregatePrimitive.listOfPrimitives.size(); i++) {
			GeometricPrimitive* currentPrimitive = aggregatePrimitive.listOfPrimitives[i];
			std::string shapeType = currentPrimitive->shape->shapeType();
			cout << shapeType << " (Shape number " << (i + 1) << " out of " << aggregatePrimitive.listOfPrimitives.size() << " total shapes)\n";

			cout << currentPrimitive->shape->printShapeInformation();

			cout << "  KA: ";
			printColor(currentPrimitive->material->constantBRDF.ka);
			cout << "  KD: ";
			printColor(currentPrimitive->material->constantBRDF.kd);
			cout << "  KR: ";
			printColor(currentPrimitive->material->constantBRDF.kr);
			cout << "  KS: ";
			printColor(currentPrimitive->material->constantBRDF.ks);
			cout << "  SP: " << currentPrimitive->material->constantBRDF.sp << "\n";
			cout << " This primitive's transformation is: \n";
			currentPrimitive->printTransformation();
		}

		std::cout << "***** FINISH PRINTING GLOBAL VARIABLES *****\n\n";
	}
}

// Prints contents of samples and buckets for debug purposes
void printSamples() {
	if (debug) {
		for (vector<Sample>::size_type i = 0; i < samples.size(); i++) {
			printSample(samples[i]);
		}
	}
}


//****************************************************
// Parsing .OBJ file specified in scene file
//****************************************************
void parseObjFile(string filename) {

	string str;
	ifstream file(filename);

	// The identifier that we're currently parsing
	string currentlyParsing;

	// current word that we're parsing on a line
	string currentWord;

	bool validLine = true;

	while (getline(file, str)) {
		// str represents the current line of the file

		validLine = true;
		int i = 0;
		istringstream iss(str);
		while (iss >> currentWord) {

			// ********** Figure out what the first word of each line is ********** //
			// We currently support:
			// (1) v ... (vertex definitions)
			// (2) f ... (face definitions)

			if ((i == 0) && (currentWord == "v")) {
				currentlyParsing = currentWord;

			} else if ((i == 0) && (currentWord == "f")) {
				currentlyParsing = currentWord;

			} else if (i == 0) {
				currentlyParsing = currentWord;
				validLine = false;
			}

			// If the current line is not valid, then just keep skipping every word in the line
			if (!validLine) {
				i++;
				continue;
			}

			// ********** After we've figured out the first word in each line, parse the rest of the line ********** //
			// If we've hit here, then we're NOT on the first word of the line anymore

			if (currentlyParsing == "v") {
				float xCoor, yCoor, zCoor;
				if (i == 0) { }
				else if (i == 1) { xCoor = stof(currentWord); }
				else if (i == 2) { yCoor = stof(currentWord); }
				else if (i == 3) { zCoor = stof(currentWord); }
				else if (i > 3) {
					cerr << "Extra parameters for " << currentlyParsing << ". Ignoring them.\n";
				}
				if (i == 3) {
					objFileVertices.push_back(Point(xCoor, yCoor, zCoor));
				}

			} else if (currentlyParsing == "f") {
				int vertexIndex1, vertexIndex2, vertexIndex3;
				if (i == 0) {}
				else if (i == 1) { vertexIndex1 = stoi(currentWord) - 1; }
				else if (i == 2) { vertexIndex2 = stoi(currentWord) - 1; }
				else if (i == 3) { vertexIndex3 = stoi(currentWord) - 1; }
				else if (i > 3) {
					cerr << "Extra parameters for " << currentlyParsing << ". Ignoring them.\n";
				}

				// Push our triangle onto the list of aggregate primitives
				if (i == 3 && vertexIndex1 < objFileVertices.size() && vertexIndex2 < objFileVertices.size() && vertexIndex3 < objFileVertices.size()) {
					GeometricPrimitive* primitiveToAdd = new GeometricPrimitive();
					BRDFCoefficients brdfToAdd = BRDFCoefficients();
					Color kaToAdd = Color(currentKar, currentKag, currentKab);
					Color kdToAdd = Color(currentKdr, currentKdg, currentKdb);
					Color ksToAdd = Color(currentKsr, currentKsg, currentKsb);
					Color krToAdd = Color(currentKrr, currentKrg, currentKrb);
					brdfToAdd.ka = kaToAdd;
					brdfToAdd.kd = kdToAdd;
					brdfToAdd.ks = ksToAdd;
					brdfToAdd.kr = krToAdd;
					brdfToAdd.sp = currentSp;
					brdfToAdd.ri = currentRi;
					Material* materialToAdd = new Material();
					materialToAdd->constantBRDF = brdfToAdd;
					primitiveToAdd->material = materialToAdd;
					primitiveToAdd->transformation = currentlySeenTransformation;
					if (currentlySeenTransformation.m.isIdentity()) {
						Triangle* triangleToAdd = new Triangle(
								objFileVertices[vertexIndex1].x, objFileVertices[vertexIndex1].y, objFileVertices[vertexIndex1].z,
								objFileVertices[vertexIndex2].x, objFileVertices[vertexIndex2].y, objFileVertices[vertexIndex2].z,
								objFileVertices[vertexIndex3].x, objFileVertices[vertexIndex3].y, objFileVertices[vertexIndex3].z);
						primitiveToAdd->shape = triangleToAdd;
					} else {
						WarpedTriangle* warpedTriangleToAdd = new WarpedTriangle(
								objFileVertices[vertexIndex1].x, objFileVertices[vertexIndex1].y, objFileVertices[vertexIndex1].z,
								objFileVertices[vertexIndex2].x, objFileVertices[vertexIndex2].y, objFileVertices[vertexIndex2].z,
								objFileVertices[vertexIndex3].x, objFileVertices[vertexIndex3].y, objFileVertices[vertexIndex3].z);
						primitiveToAdd->shape = warpedTriangleToAdd;
					}

					aggregatePrimitive.addPrimitive(primitiveToAdd);
				}
			}

			i++;
		}
	}
}



//****************************************************
// Parsing scene file from command line
//****************************************************

void parseSceneFile(string filename) {

	ifstream file(filename);
	string str;

	// The identifier that we're currently parsing
	string currentlyParsing;

	// current word that we're parsing on a line
	string currentWord;

	// Index of a word on each specific line
	int i = 0;

	// Indicates whether the current line being read is valid
	bool validLine = true;

	// ***** Variables for object creation *****//
	// Camera
	float ex, ey, ez, llx, lly, llz, lrx, lry, lrz, ulx, uly, ulz, urx, ury, urz;
	// Most recently seen material
	float kar, kag, kab, kdr, kdg, kdb, ksr, ksg, ksb, ksp, krr, krg, krb;
	float ri = 0.0;

	// .obj filename that we've found
	string objFilename;
	bool foundObjFile = false;

	while (getline(file, str)) {
		// str represents the current line of the file

		validLine = true;
		i = 0;
		istringstream iss(str);
		while (iss >> currentWord) {

			// ********** Figure out what the first word of each line is ********** //
			if ((i == 0) && (currentWord == "cam")) {
				currentlyParsing = currentWord;
			} else if ((i == 0) && (currentWord == "ltd")) {
				currentlyParsing = currentWord;
			} else if ((i == 0) && (currentWord == "ltp")) {
				currentlyParsing = currentWord;
			} else if ((i == 0) && (currentWord == "lta")) {
				currentlyParsing = currentWord;
			} else if ((i == 0) && (currentWord == "mat")) {
				currentlyParsing = currentWord;
			} else if ((i == 0) && (currentWord == "sph")) {
				currentlyParsing = currentWord;
			} else if ((i == 0) && (currentWord == "tri")) {
				currentlyParsing = currentWord;
			} else if ((i == 0) && (currentWord == "obj")) {
				currentlyParsing = currentWord;
			} else if ((i == 0) && (currentWord == "xft")) {
				currentlyParsing = currentWord;
			} else if ((i == 0) && (currentWord == "xfs")) {
				currentlyParsing = currentWord;
			} else if ((i == 0) && (currentWord == "xfz")) {
				currentlyParsing = currentWord;
			} else if ((i == 0) && (currentWord == "xfr")) {
				currentlyParsing = currentWord;

				// We've found an unspecified identifier as the first word on our line
			} else if (i == 0) {
				currentlyParsing = currentWord;
				validLine = false;
			}

			// If the current line is not valid, then just keep skipping every word in the line
			if (!validLine) {
				i++;
				continue;
			}


			// ********** After we've figured out the first word in each line, parse the rest of the line ********** //
			// If we've hit here, then we're NOT on the first word of the line anymore

			// Camera must be the first line in the scene file
			if (currentlyParsing == "cam") {
				if (i == 0) { }
				else if (i == 1) { ex = stof(currentWord); }
				else if (i == 2) { ey = stof(currentWord); }
				else if (i == 3) { ez = stof(currentWord); }
				else if (i == 4) { llx = stof(currentWord); }
				else if (i == 5) { lly = stof(currentWord); }
				else if (i == 6) { llz = stof(currentWord); }
				else if (i == 7) { lrx = stof(currentWord); }
				else if (i == 8) { lry = stof(currentWord); }
				else if (i == 9) { lrz = stof(currentWord); }
				else if (i == 10) { ulx = stof(currentWord); }
				else if (i == 11) { uly = stof(currentWord); }
				else if (i == 12) { ulz = stof(currentWord); }
				else if (i == 13) { urx = stof(currentWord); }
				else if (i == 14) { ury = stof(currentWord); }
				else if (i == 15) { urz = stof(currentWord); }
				else { cerr << "Extra parameters for " << currentlyParsing << ". Ignoring them.\n"; }

			} else if (currentlyParsing == "ltd") {
				float dx, dy, dz, r, g, b;
				if (i == 0) { }
				else if (i == 1) { dx = stof(currentWord); }
				else if (i == 2) { dy = stof(currentWord); }
				else if (i == 3) { dz = stof(currentWord); }
				else if (i == 4) { r = stof(currentWord); }
				else if (i == 5) { g = stof(currentWord); }
				else if (i == 6) { b = stof(currentWord); }
				else if (i > 6) {
					cerr << "Extra parameters for " << currentlyParsing << ". Ignoring them.\n";
				}
				if (i == 6) {
					// Add directional light to global list
					directional_lights.push_back(DirectionalLight(dx, dy, dz, r, g, b));
				}

			} else if (currentlyParsing == "ltp") {
				float px, py, pz, r, g, b, falloff;
				if (i == 0) { }
				else if (i == 1) { px = stof(currentWord); }
				else if (i == 2) { py = stof(currentWord); }
				else if (i == 3) { pz = stof(currentWord); }
				else if (i == 4) { r = stof(currentWord); }
				else if (i == 5) { g = stof(currentWord); }
				else if (i == 6) { b = stof(currentWord); }
				else if (i == 7) { falloff = stof(currentWord); }
				else if (i > 7) {
					cerr << "Extra parameters for " << currentlyParsing << ". Ignoring them.\n";
				}
				if (i == 6) {
					if (!softShadows) {
						// Add point light to global list
						point_lights.push_back(PointLight(px, py, pz, r, g, b, falloff));
					} else {
						// We want soft shadows, so generate lots of point lights

						float horizontalLightOffset = (urx - ulx) / film.width * 8;
						float verticalLightOffset = (ury - lry) / film.height * 8;

						point_lights.push_back(PointLight(px - 3 * horizontalLightOffset, py + 3.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px - 2 * horizontalLightOffset, py + 3.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px - horizontalLightOffset, py + 3.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px, py + 3.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + horizontalLightOffset, py + 3.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + 2 * horizontalLightOffset, py + 3.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + 3 * horizontalLightOffset, py + 3.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));

						point_lights.push_back(PointLight(px - 3 * horizontalLightOffset, py + 2.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px - 2 * horizontalLightOffset, py + 2.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px - horizontalLightOffset, py + 2.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px, py + 2.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + horizontalLightOffset, py + 2.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + 2 * horizontalLightOffset, py + 2.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + 3 * horizontalLightOffset, py + 2.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));

						point_lights.push_back(PointLight(px - 3 * horizontalLightOffset, py + verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px - 2 * horizontalLightOffset, py + verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px - horizontalLightOffset, py + verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px, py + verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + horizontalLightOffset, py + verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + 2 * horizontalLightOffset, py + verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + 3 * horizontalLightOffset, py + verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));

						point_lights.push_back(PointLight(px - 3 * horizontalLightOffset, py, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px - 2 * horizontalLightOffset, py, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px - horizontalLightOffset, py, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px, py, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + horizontalLightOffset, py, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + 2 * horizontalLightOffset, py, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + 3 * horizontalLightOffset, py, pz, r / 49.0, g / 49.0, b / 49.0, falloff));

						point_lights.push_back(PointLight(px - 3 * horizontalLightOffset, py - verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px - 2 * horizontalLightOffset, py - verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px - horizontalLightOffset, py - verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px, py - verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + horizontalLightOffset, py - verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + 2 * horizontalLightOffset, py - verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + 3 * horizontalLightOffset, py - verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));

						point_lights.push_back(PointLight(px - 3 * horizontalLightOffset, py - 2.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px - 2 * horizontalLightOffset, py - 2.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px - horizontalLightOffset, py - 2.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px, py - 2.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + horizontalLightOffset, py - 2.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + 2 * horizontalLightOffset, py - 2.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + 3 * horizontalLightOffset, py - 2.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));

						point_lights.push_back(PointLight(px - 3 * horizontalLightOffset, py - 3.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px - 2 * horizontalLightOffset, py - 3.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px - horizontalLightOffset, py - 3.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px, py - 3.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + horizontalLightOffset, py - 3.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + 2 * horizontalLightOffset, py - 3.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
						point_lights.push_back(PointLight(px + 3 * horizontalLightOffset, py - 3.0 * verticalLightOffset, pz, r / 49.0, g / 49.0, b / 49.0, falloff));
					}
				}
			} else if (currentlyParsing == "lta") {
				float r, g, b;
				if (i == 0) { }
				else if (i == 1) { r = stof(currentWord); }
				else if (i == 2) { g = stof(currentWord); }
				else if (i == 3) { b = stof(currentWord); }
				else if (i > 3) {
					cerr << "Extra parameters for " << currentlyParsing << ". Ignoring them.\n";
				}
				if (i == 3) {
					// Add ambient light to ambient light list
					ambient_lights.push_back(AmbientLight(r, g, b));
				}

			} else if (currentlyParsing == "mat") {
				if (i == 0) { }
				else if (i == 1) { kar = stof(currentWord); }
				else if (i == 2) { kag = stof(currentWord); }
				else if (i == 3) { kab = stof(currentWord); }
				else if (i == 4) { kdr = stof(currentWord); }
				else if (i == 5) { kdg = stof(currentWord); }
				else if (i == 6) { kdb = stof(currentWord); }
				else if (i == 7) { ksr = stof(currentWord); }
				else if (i == 8) { ksg = stof(currentWord); }
				else if (i == 9) { ksb = stof(currentWord); }
				else if (i == 10) { ksp = stof(currentWord); }
				else if (i == 11) { krr = stof(currentWord); }
				else if (i == 12) { krg = stof(currentWord); }
				else if (i == 13) { krb = stof(currentWord); }
				else if (i == 14) { ri = stof(currentWord); }
				else { cerr << "Extra parameters for " << currentlyParsing << ". Ignoring them. i is : " << i << "\n"; }
				if (i == 13 || i == 14) {
					currentKar = kar;
					currentKag = kag;
					currentKab = kab;
					currentKdr = kdr;
					currentKdg = kdg;
					currentKdb = kdb;
					currentKsr = ksr;
					currentKsg = ksg;
					currentKsb = ksb;
					currentKrr = krr;
					currentKrg = krg;
					currentKrb = krb;
					currentSp = ksp;
					currentRi = ri;
				}

			} else if (currentlyParsing == "sph") {
				float cx, cy, cz, r;
				if (i == 0) { }
				else if (i == 1) { cx = stof(currentWord); }
				else if (i == 2) { cy = stof(currentWord); }
				else if (i == 3) { cz = stof(currentWord); }
				else if (i == 4) { r = stof(currentWord); }
				else if (i > 4) {
					cerr << "Extra parameters for " << currentlyParsing << ". Ignoring them. i is : " << i << "\n";
				}
				if (i == 4) {
					// Initialize our object, with the last seen material
					GeometricPrimitive* primitiveToAdd = new GeometricPrimitive();
					BRDFCoefficients brdfToAdd = BRDFCoefficients();
					Color kaToAdd = Color(currentKar, currentKag, currentKab);
					Color kdToAdd = Color(currentKdr, currentKdg, currentKdb);
					Color ksToAdd = Color(currentKsr, currentKsg, currentKsb);
					Color krToAdd = Color(currentKrr, currentKrg, currentKrb);
					brdfToAdd.ka = kaToAdd;
					brdfToAdd.kd = kdToAdd;
					brdfToAdd.ks = ksToAdd;
					brdfToAdd.kr = krToAdd;
					brdfToAdd.ri = currentRi;
					brdfToAdd.sp = currentSp;
					Material* materialToAdd = new Material();
					materialToAdd->constantBRDF = brdfToAdd;
					primitiveToAdd->material = materialToAdd;
					primitiveToAdd->transformation = currentlySeenTransformation;

					if (!currentlySeenTransformation.m.isIdentity()) {
						Ellipsoid* ellipsoidToAdd = new Ellipsoid(cx, cy, cz, r);
						primitiveToAdd->shape = ellipsoidToAdd;
					} else {
						Sphere* sphereToAdd = new Sphere(cx, cy, cz, r);
						primitiveToAdd->shape = sphereToAdd;
					}
					aggregatePrimitive.addPrimitive(primitiveToAdd);
				}

			} else if (currentlyParsing == "tri") {
				float ax, ay, az, bx, by, bz, cx, cy, cz;
				if (i == 0) { }
				else if (i == 1) { ax = stof(currentWord); }
				else if (i == 2) { ay = stof(currentWord); }
				else if (i == 3) { az = stof(currentWord); }
				else if (i == 4) { bx = stof(currentWord); }
				else if (i == 5) { by = stof(currentWord); }
				else if (i == 6) { bz = stof(currentWord); }
				else if (i == 7) { cx = stof(currentWord); }
				else if (i == 8) { cy = stof(currentWord); }
				else if (i == 9) { cz = stof(currentWord); }
				else if (i > 9) {
					cerr << "Extra parameters for " << currentlyParsing << ". Ignoring them. i is : " << i << "\n";
				}
				if (i == 9) {
					GeometricPrimitive* primitiveToAdd = new GeometricPrimitive();
					BRDFCoefficients brdfToAdd = BRDFCoefficients();
					Color kaToAdd = Color(currentKar, currentKag, currentKab);
					Color kdToAdd = Color(currentKdr, currentKdg, currentKdb);
					Color ksToAdd = Color(currentKsr, currentKsg, currentKsb);
					Color krToAdd = Color(currentKrr, currentKrg, currentKrb);
					brdfToAdd.ka = kaToAdd;
					brdfToAdd.kd = kdToAdd;
					brdfToAdd.ks = ksToAdd;
					brdfToAdd.kr = krToAdd;
					brdfToAdd.sp = currentSp;
					brdfToAdd.ri = currentRi;
					Material* materialToAdd = new Material();
					materialToAdd->constantBRDF = brdfToAdd;
					primitiveToAdd->material = materialToAdd;
					primitiveToAdd->transformation = currentlySeenTransformation;

					if (!currentlySeenTransformation.m.isIdentity()) {
						WarpedTriangle* warpedTriangleToAdd = new WarpedTriangle(ax, ay, az, bx, by, bz, cx, cy, cz);
						primitiveToAdd->shape = warpedTriangleToAdd;
					} else {
						Triangle* triangleToAdd = new Triangle(ax, ay, az, bx, by, bz, cx, cy, cz);
						primitiveToAdd->shape = triangleToAdd;
					}
					aggregatePrimitive.addPrimitive(primitiveToAdd);
				}
			} else if (currentlyParsing == "obj") {
				if (i == 0) { }
				else if (i == 1) {
					objFilename = currentWord;
					foundObjFile = true;
				} else if (i > 1) {
					cerr << "Extra parameters for " << currentlyParsing << ". Ignoring them. i is : " << i << "\n";
				}
			} else if (currentlyParsing == "xft") {
				float x, y, z;
				if (i == 0) { }
				else if (i == 1) {x = stof(currentWord); }
				else if (i == 2) {y = stof(currentWord); }
				else if (i == 3) {z = stof(currentWord); }
				else if (i > 3) {
					cerr << "Extra parameters for " << currentlyParsing << ". Ignoring them. i is : " << i << "\n";
				}
				if (i == 3) {
					currentlySeenTransformation.appendTransformation(Matrix4::createTranslationMatrix(x, y, z));
				}
			} else if (currentlyParsing == "xfs") {
				float x, y, z;
				if (i == 0) { }
				else if (i == 1) {x = stof(currentWord); }
				else if (i == 2) {y = stof(currentWord); }
				else if (i == 3) {z = stof(currentWord); }
				else if (i > 3) {
					cerr << "Extra parameters for " << currentlyParsing << ". Ignoring them. i is : " << i << "\n";
				}
				if (i == 3) {
					currentlySeenTransformation.appendTransformation(Matrix4::createScalingMatrix(x, y, z));
				}
			} else if (currentlyParsing == "xfz") {
				if (i == 0) { currentlySeenTransformation.resetTransformation(); }
				else if (i > 0) {
					cerr << "Extra parameters for " << currentlyParsing << ". Ignoring them. i is : " << i << "\n";
				}
			} else if (currentlyParsing == "xfr") {
				float x, y, z;
				if (i == 0) { }
				else if (i == 1) {x = stof(currentWord); }
				else if (i == 2) {y = stof(currentWord); }
				else if (i == 3) {z = stof(currentWord); }
				else if (i > 3) {
					cerr << "Extra parameters for " << currentlyParsing << ". Ignoring them. i is : " << i << "\n";
				}
				if (i == 3) {
					currentlySeenTransformation.appendTransformation(Matrix4::createRotationMatrix(x, y, z));
				}
			}


			// TODO: Add more parsing here

			i++;
		}

		if (!validLine) {
			cerr << "Unsupported feature: " << currentlyParsing << ". Ignoring line.\n";
		}

	}

	if (foundObjFile) {
		parseObjFile(objFilename);
	}

	// ***** Initialize Camera global variable ***** //
	// NOTE: Film width and height MUST already be initialized
	camera = Camera(film.width, film.height, ex, ey, ez, llx, lly, llz, lrx, lry, lrz, urx, ury, urz, ulx, uly, ulz);


}



//****************************************************
// Parsing of command line options, with options:
// (1) -dimensions width height
//     adds viewport width and height attributes to Viewport global variable
// (2) -depth n
//     sets recursion depth
//
// NOTE: first command line option must be .scene file
//
// NOTE: also performs proper initialization of Film global variable
//****************************************************
void parseCommandLineOptions(int argc, char *argv[]) {

	string flag;

	int i = 2;
	while (i <= argc - 1) {
		flag = argv[i];

		if (flag == "-dimensions") {
			// Check that -dimensions has enough option parameters
			if ((i + 2) > (argc - 1)) {
				std::cout << "Invalid number of parameters for -dimensions.";
				exit(1);
			}

			int widthOfFilm = stoi(argv[i+1]);
			int heightOfFilm = stoi(argv[i+2]);

			//	  if (widthOfFilm < 1000 || heightOfFilm < 500 || widthOfFilm > 3000 || heightOfFilm > 3000) {
			//		  std::cout << "Dimensions of output file must be at least 1000x500 and no more than 3000x3000.";
			//		  exit(1);
			//	  }

			film = Film(widthOfFilm, heightOfFilm);
			i += 2;

		} else if (flag == "-depth") {
			// Check that -depth has enough option parameters
			if ((i + 1) > (argc - 1))
			{
				std::cout << "Invalid number of parameters for -depth.";
				exit(1);
			}
			recursionDepth = stoi(argv[i+1]);
			i += 1;
		} else if (flag == "-ss") {
			// Check that -depth has enough option parameters
			if (i > (argc - 1))
			{
				std::cout << "Invalid number of parameters for -ss.";
				exit(1);
			}
			softShadows = true;
		} else if (flag == "-distributed") {
			// Check that -depth has enough option parameters
			if (i > (argc - 1))
			{
				std::cout << "Invalid number of parameters for -distributed.";
				exit(1);
			}
			distributedRayTracing = true;
		} else if (flag == "-filename") {
			// Check that -depth has enough option parameters
			if ((i + 1) > (argc - 1))
			{
				std::cout << "Invalid number of parameters for -depth.";
				exit(1);
			}
			filename = argv[i+1];
			i += 1;
		} else if (flag == "-checkerboard") {
			// Check that -depth has enough option parameters
			if (i > (argc - 1))
			{
				std::cout << "Invalid number of parameters for -checkerboard.";
				exit(1);
			}
			checkerboard = true;
		}
		else {
			std::cout << "Extra parameters in command line options; terminating program.";
			exit(1);

		}

		// Advance to next flag, if one exists
		i++;
	}

	// Parse scene file
	parseSceneFile(argv[1]);
}



//****************************************************
// Main rendering loop
//
// Loops through samples, and does the following per sample:
// (1) generates a ray from the eye through the sample
// (2) traces this ray with the ray tracer
// (3) commits the color returns by the ray tracer to the film
//****************************************************
void render() {
	// Loop through all of the samples...
	for (vector<Sample>::size_type i = 0; i < samples.size(); i++) {

		if (debug && i % 1000 == 0) {
			std::cout << "Currently processing sample " << i << " out of " << samples.size() << ".\n";
		}

		// For each sample, generate a ray from the eye to the sample location
		Ray currentRay;
		Color currentSampleColor;

		if (distributedRayTracing) {
			for (int j = 0; j < 16; j++) {
				camera.generateRay(samples[i], &currentRay, distributedRayTracing);
				rayTracer.trace(currentRay, 0, &currentSampleColor);
				film.commitColor(samples[i], currentSampleColor);
			}
		} else {
			// Given the sample in Film-coordinates, tell the camera to generate a viewing ray in IMAGE PLANE [-1, 1] coordinates
			camera.generateRay(samples[i], &currentRay, distributedRayTracing);

			// Call the trace method to try to populate currentSampleColor for the currentSample
			rayTracer.trace(currentRay, 0, &currentSampleColor);

			// Commit the currentSampleColor for the currentSample onto our Film
			film.commitColor(samples[i], currentSampleColor);
		}
	}

}


//****************************************************
// Populates a list of samples based on the film height and width
//****************************************************
void initializeSampler() {
	float x, y;
	// Generates samples in the FILM's coordinates
	for (y = 0; y < film.height; y++) {
		for (x = 0; x < film.width; x++) {
			samples.push_back(Sample(x, y));
		}
	}

}



//****************************************************
// Main function
//****************************************************
int main(int argc, char *argv[]) {

	// Default the filename to 'output_image.png' if none is specified
	if (filename == NULL) {
		filename = "output_image.png";
	}

	// Turns debug mode ON or OFF
	debug = true;

	// Parse command line options
	parseCommandLineOptions(argc, argv);
	printGlobalVariables();

	// Initializes list of buckets; Buckets have a list of samples
	initializeSampler();

	render();

	// Deallocate memory
	for (std::vector<GeometricPrimitive*>::size_type i = 0; i < aggregatePrimitive.listOfPrimitives.size(); i++) {
		delete aggregatePrimitive.listOfPrimitives[i]->material;
		delete aggregatePrimitive.listOfPrimitives[i]->shape;
		delete aggregatePrimitive.listOfPrimitives[i];
	}

	film.writeImage(filename);

	return 0;
};

