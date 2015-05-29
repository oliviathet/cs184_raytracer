#ifndef SHAPE_H_
#define SHAPE_H_

#include "Matrix4.h"

static void printRay(Ray ray);
static void printVector(Vector3 vector);

// Converts float to string
std::string convertFloatToString(float number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

// NOTE: A Primitive object has both a Material and a Shape
class Shape {
	public:
		// Abstract methods that Sphere and Triangle inherit from.
		virtual bool intersect(Ray &ray, float* tHit, DifferentialGeometry* differentialGeometry) {
			return false;
		}
		virtual bool intersectP(Ray &ray) {
			return false;
		}
		virtual std::string shapeType() {
			return "";
		}
		virtual std::string printShapeInformation() {
			return "";
		}


};

class Sphere : public Shape {
	public:
		float x, y, z, r;

	Sphere() {

	}

	Sphere(float x, float y, float z, float r) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->r = r;
	}

	// Test if ray intersects with the shape or not (in object space).
	// If so, return intersection point and normal.
	// This method is passed up to the geometric primitive class, which owns a method
	// of similar signature. Only difference is that it takes in an Intersection object
	// instead of a DifferentialGeometry object.
	// NOTE: Intersection object comprises of a DifferentialGeometry and a Primitive.
	// NOTE: This method must populate tHit and differentialGeometry.
	bool intersect(Ray &ray, float* tHit, DifferentialGeometry* differentialGeometry) {

		// Sphere equation: (x - x0)^2 + (y - y0)^2 + (z - z0)^2 = r^2

		// Ray equations
		// x(t) = ft + g
		float f = ray.direction.x;
		float g = ray.position.x;

		// y(t) = ht + i
		float h = ray.direction.y;
		float i = ray.position.y;

		// z(t) = jt + k
		float j = ray.direction.z;
		float k = ray.position.z;

		// (ft + m))^2 + (ht + n))^2 + (jt + p))^2 - r^2 = 0
		float m = g - this->x;
		float n = i - this->y;
		float p = k - this->z;

		// at^2 + bt + c = 0
		float a = pow(f, 2) + pow(h, 2) + pow(j, 2);
		float b = (2 * f * m) + (2 * h * n) + (2 * j * p);
		float c = pow(m, 2) + pow(n, 2) + pow(p, 2) - pow(this->r, 2);

		float discriminant = pow(b, 2) - (4 * a * c);

		// No intersection, so return false
		if (discriminant < 0) {
			return false;
		}

		// Quadratic equation: (-b +/- sqrt(b^2 - 4ac)) / 2a
		float sqrtTerm = pow(b, 2) - (4 * a * c);

		float intersectionPoint1 = FLT_MAX;
		float intersectionPoint2 = FLT_MAX;
		float tIntersection = FLT_MAX;

		if (sqrtTerm == 0) {
			tIntersection = (b * -1.0) / (2.0 * a);
			// if our ONE intersection point is invalid, then return false
			if (!(tIntersection > ray.t_min && tIntersection < ray.t_max)) {
				return false;
			}
		// We have TWO intersection points
		} else {
			// Calculate both intersection points
			intersectionPoint1 = ((b * -1.0) + pow(sqrtTerm, 0.5)) / (2.0 * a);
			intersectionPoint2 = ((b * -1.0) - pow(sqrtTerm, 0.5)) / (2.0 * a);

			bool intersectionPoint1Valid = false;
			bool intersectionPoint2Valid = false;

			// If intersectionPoint1 is valid, according to our t_min and t_max...
			if (intersectionPoint1 > ray.t_min && intersectionPoint1 < ray.t_max) {
				intersectionPoint1Valid = true;
			}

			// If intersectionPoint1 is valid, according to our t_min and t_max...
			if (intersectionPoint2 > ray.t_min && intersectionPoint2 < ray.t_max) {
				intersectionPoint2Valid = true;
			}

			// Both points are valid
			if (intersectionPoint1Valid && intersectionPoint2Valid) {
				tIntersection = fmin(intersectionPoint1, intersectionPoint2);
			// Only intersection point 1 is valid
			} else if (intersectionPoint1Valid && !intersectionPoint2Valid) {
				tIntersection = intersectionPoint1;
			// Only intersection point 2 is valid
			} else if (!intersectionPoint1Valid && intersectionPoint2Valid) {
				tIntersection = intersectionPoint2;
			// Neither of our points are valid
			} else {
				return false;
			}
		}

		// NOTE: At this point, tIntersection is our closest intersection point

		// Now, find the (xCoor, yCoor, zCoor) that tIntersection corresponds to
		float xCoor = (f * tIntersection) + g;
		float yCoor = (h * tIntersection) + i;
		float zCoor = (j * tIntersection) + k;

		// Compute normal of the sphere and the intersection point in WORLD coordinates
		// This is NOT equal to the position of the intersection point

		differentialGeometry->normal = Normal::normalizeNormal(Normal(xCoor - this->x, yCoor - this->y, zCoor - this->z));

		differentialGeometry->position = Point(xCoor, yCoor, zCoor);

		// IMPORTANT: this used to be:
		// tHit = &tIntersection
		*tHit = tIntersection;
		return true;

	}

	// Test if ray intersects with the shape or not.
	// This method is passed up to the geometric primitive class, which owns a method
	// of the same signature.
	bool intersectP(Ray &ray) {
		// Ray equations
		// x(t) = ft + g
		float f = ray.direction.x;
		float g = ray.position.x;

		// y(t) = ht + i
		float h = ray.direction.y;
		float i = ray.position.y;

		// z(t) = jt + k
		float j = ray.direction.z;
		float k = ray.position.z;

		// (ft + m))^2 + (ht + n))^2 + (jt + p))^2 - r^2 = 0
		float m = g - this->x;
		float n = i - this->y;
		float p = k - this->z;

		// at^2 + bt + c = 0
		float a = pow(f, 2) + pow(h, 2) + pow(j, 2);
		float b = (2 * f * m) + (2 * h * n) + (2 * j * p);
		float c = pow(m, 2) + pow(n, 2) + pow(p, 2) - pow(this->r, 2);

		// discriminant = b^2 - 4ac
		float discriminant = pow(b, 2) - (4 * a * c);

		if (discriminant < 0) {
			return false;
		}

		// ***** Still need to check whether tIntersection is greater than tMin *****
		// Quadratic equation: (-b +/- sqrt(b^2 - 4ac)) / 2a
		float sqrtTerm = pow(b, 2) - (4 * a * c);

		float intersectionPoint1, intersectionPoint2;

		if (sqrtTerm == 0) {
			intersectionPoint1 = (b * -1.0) / (2.0 * a);
		} else {
			// sqrtTerm is positive
			intersectionPoint1 = ((b * -1.0) + pow(sqrtTerm, 0.5)) / (2.0 * a);
			intersectionPoint2 = ((b * -1.0) - pow(sqrtTerm, 0.5)) / (2.0 * a);
		}

		if (intersectionPoint1 >= ray.t_min && intersectionPoint1 < ray.t_max) {
			return true;
		}
		if (intersectionPoint2 >= ray.t_min && intersectionPoint2 < ray.t_max) {
			return true;
		}

		return false;
	}

	std::string shapeType() {
		return "Sphere";
	}

	std::string printShapeInformation() {
		std::string xCoor = convertFloatToString(x);
		std::string yCoor = convertFloatToString(y);
		std::string zCoor = convertFloatToString(z);
		std::string radius = convertFloatToString(r);

		std::string toReturn = "Center is (";
		toReturn.append(xCoor).append(", ").append(yCoor).append(", ").append(zCoor).append(") with radius ").append(radius).append("\n");
		return toReturn;
	}
};

class Ellipsoid : public Sphere {
public:
	float x, y, z, r;

	Ellipsoid() : Sphere() {
	}

	Ellipsoid(float x, float y, float z, float r) : Sphere(x, y, z, r) {

	}

	// If these two methods are called, then ray represents the transformed ray that has
	// been transformed into object coordinates; i.e. original sphere coordinates.
	bool intersect(Ray &ray, float* tHit, DifferentialGeometry* differentialGeometry) {
		return Sphere::intersect(ray, tHit, differentialGeometry);
	}

	bool intersectP(Ray &ray) {
		return Sphere::intersectP(ray);
	}

	std::string shapeType() {
		return "Ellipsoid";
	}

};

class Triangle : public Shape {
	public:
		Point v1, v2, v3;

	Triangle() {

	}

	Triangle(Point v1, Point v2, Point v3) {
		this->v1 = v1;
		this->v2 = v2;
		this->v3 = v3;
	}

	Triangle(float ax, float ay, float az,
			float bx, float by, float bz,
			float cx, float cy, float cz) {
		this->v1 = Point(ax, ay, az);
		this->v2 = Point(bx, by, bz);
		this->v3 = Point(cx, cy, cz);
	}

	Normal computeNormal() {
		// Calculate triangle edge vectors to calculate the normal

		// u = T.v2 - T.v1
		Vector3 u = Vector3(this->v2.x - this->v1.x, this->v2.y - this->v1.y, this->v2.z - this->v1.z);
		// v = T.v3 - T.v1
		Vector3 v = Vector3(this->v3.x - this->v1.x, this->v3.y - this->v1.y, this->v3.z - this->v1.z);

		// normal = u x v
		Vector3 normal = u.crossProduct(v);
		// (b - a) x (c - a)
//		return Normal::normalizeNormal(Normal(normal.x, normal.y, normal.z));

		// (a - b) x (a - c)
		return Normal::normalizeNormal(Normal(-1.0 * normal.x, -1.0 * normal.y, -1.0 * normal.z));
	}

	// Test if ray intersects with the shape or not (in object space).
	// If so, return intersection point and normal.
	// This method is passed up to the geometric primitive class, which owns a method
	// of similar signature. Only difference is that it takes in an Intersection object
	// instead of a DifferentialGeometry object.
	// NOTE: Intersection object comprises of a DifferentialGeometry and a Primitive.
	// NOTE: This method must populate tHit and differentialGeometry.
	bool intersect(Ray &ray, float* tHit, DifferentialGeometry* differentialGeometry) {

		// The triangle intersection occurs when:
		// F(t) = v1 + intersectionPoint1(v2 - v1) + intersectionPoint2(v3 - v1)
		// for some t, intersectionPoint1, intersectionPoint2. The intersection p = F(t)
		// intersectionPoint1 > 0, intersectionPoint2 > 0, intersectionPoint1 + intersectionPoint2 < 1 must be satisfied, otherwise the ray has hit ...
		//the plane outside of the triangle

		// Ray equations
		// x(t) = ft + g
		float f = ray.direction.x;
		float g = ray.position.x;

		// y(t) = ht + i
		float h = ray.direction.y;
		float i = ray.position.y;

		// z(t) = jt + k
		float j = ray.direction.z;
		float k = ray.position.z;


		// Calculate the triangle edge vectors
		// u = T.v2 - T.v1
		Vector3 u = Vector3(this->v2.x - this->v1.x, this->v2.y - this->v1.y, this->v2.z - this->v1.z);
		// v = T.v3 - T.v1
		Vector3 v = Vector3(this->v3.x - this->v1.x, this->v3.y - this->v1.y, this->v3.z - this->v1.z);

		// Calculuate plane normal
		// normal = u x v
		Normal normal = this->computeNormal();
		Vector3 normalVector = Vector3(normal.x, normal.y, normal.z);

		// If normal = (0, 0, 0), the triangle is degenerate, so return false
		if (normal.x == 0 && normal.y == 0 && normal.z == 0) {
			return false;
		}

		Vector3 rayTrianglePositionVector = Vector3(g - this->v1.x, i - this->v1.y, k - this->v1.z);

		Vector3 rayDirectionVector = Vector3(f, h, j);

		float rayInsideTriangle = -1.0 * (normalVector.dotProduct(rayTrianglePositionVector));
		float rayParallelTriangle = normalVector.dotProduct(rayDirectionVector);

		// If the ray is parallel to the triangle, return false
		if (std::fabs(rayParallelTriangle) < 0.001) {
			return false;
		}

		float tIntersection = rayInsideTriangle/rayParallelTriangle;
		// If the ray is heading away from the triangle, return false
		if (tIntersection < 0.0) {
			return false;
		}

		// Now, find the (xCoor, yCoor, zCoor) that tIntersection corresponds to
		float xCoor = (f * tIntersection) + g;
		float yCoor = (h * tIntersection) + i;
		float zCoor = (j * tIntersection) + k;

		// Is the intersection inside the triangle?
		float uu = u.dotProduct(u);
		float uv = u.dotProduct(v);
		float vv = v.dotProduct(v);

		Vector3 w = Vector3(xCoor - this->v1.x, yCoor - this->v1.y, zCoor - this->v1.z);

		float wu = w.dotProduct(u);
		float wv = w.dotProduct(v);

		float discriminant = uv * uv - uu * vv;

		// Check parametric coordinates
		if (tIntersection < ray.t_min || tIntersection > ray.t_max) {
			return false;
		}
		float intersectionPoint1 = (uv * wv - vv * wu) / discriminant;
		float intersectionPoint2 = (uv * wu - uu * wv) / discriminant;

		if (intersectionPoint1 < 0.0 || intersectionPoint1 > 1.0) {
			return false;
		}
		if (intersectionPoint2 < 0.0 || (intersectionPoint1 + intersectionPoint2) > 1.0) {
			return false;
		}

		differentialGeometry->normal = normal;
		differentialGeometry->position = Point(xCoor, yCoor, zCoor);

		*tHit = tIntersection;
		return true;
	}

	//Test if ray intersects with the shape or not.
	//This method is passed up to the geometric primitive class, which owns a method
	//of the same signature.
	bool intersectP(Ray &ray) {
		// Ray equations
		// x(t) = ft + g
		float f = ray.direction.x;
		float g = ray.position.x;

		// y(t) = ht + i
		float h = ray.direction.y;
		float i = ray.position.y;

		// z(t) = jt + k
		float j = ray.direction.z;
		float k = ray.position.z;


		// Calculate the triangle edge vectors
		// u = T.v2 - T.v1
		Vector3 u = Vector3(this->v2.x - this->v1.x, this->v2.y - this->v1.y, this->v2.z - this->v1.z);
		// v = T.v3 - T.v1
		Vector3 v = Vector3(this->v3.x - this->v1.x, this->v3.y - this->v1.y, this->v3.z - this->v1.z);

		// Calculuate plane normal
		// normal = u x v
		Normal normal = this->computeNormal();
		Vector3 normalVector = Vector3(normal.x, normal.y, normal.z);

		// If normal = (0, 0, 0), the triangle is degenerate, so return false
		if (normal.x == 0 && normal.y == 0 && normal.z == 0) {
			return false;
		}

		// rayTrianglePositionVector = StartingRayPosition - TriangleV1Position
		Vector3 rayTrianglePositionVector = Vector3(g - this->v1.x, i - this->v1.y, k - this->v1.z);

		// rayDirectionVector = RayDirection
		Vector3 rayDirectionVector = Vector3(f, h, j);

		float rayInsideTriangle = -1.0 * normalVector.dotProduct(rayTrianglePositionVector);
		float rayParallelTriangle = normalVector.dotProduct(rayDirectionVector);

		// If the ray is parallel to the triangle, return false
		if (std::fabs(rayParallelTriangle) < 0.001) {
			return false;
		}

		float tIntersection = rayInsideTriangle/rayParallelTriangle;
		// If the ray is heading away from the triangle, return false
		if (tIntersection < 0.0) {
			return false;
		}

		// Now, find the (xCoor, yCoor, zCoor) that tIntersection corresponds to
		float xCoor = (f * tIntersection) + g;
		float yCoor = (h * tIntersection) + i;
		float zCoor = (j * tIntersection) + k;

		// Is the intersection inside the triangle?
		float uu = u.dotProduct(u);
		float uv = u.dotProduct(v);
		float vv = v.dotProduct(v);

		Vector3 w = Vector3(this->v1.x - xCoor, this->v1.y - yCoor, this->v1.z - zCoor);

		float wu = w.dotProduct(u);
		float wv = w.dotProduct(v);

		float discriminant = uv * uv - uu * vv;

		// Check parametric coordinates
		if (tIntersection < ray.t_min || tIntersection > ray.t_max) {
			return false;
		}
		float intersectionPoint1 = (uv * wv - vv * wu) / discriminant;
		float intersectionPoint2 = (uv * wu - uu * wv) / discriminant;

		if (intersectionPoint1 < 0.0 || intersectionPoint1 > 1.0) {
			return false;
		}
		if (intersectionPoint2 < 0.0 || (intersectionPoint1 + intersectionPoint2) > 1.0) {
			return false;
		}

		return true;
	}

	std::string shapeType() {
		return "Triangle";
	}

	std::string printShapeInformation() {
		float ax = v1.x;
		float ay = v1.y;
		float az = v1.z;
		float bx = v2.x;
		float by = v2.y;
		float bz = v2.z;
		float cx = v3.x;
		float cy = v3.y;
		float cz = v3.z;

		std::string axString = convertFloatToString(ax);
		std::string ayString = convertFloatToString(ay);
		std::string azString = convertFloatToString(az);
		std::string bxString = convertFloatToString(bx);
		std::string byString = convertFloatToString(by);
		std::string bzString = convertFloatToString(bz);
		std::string cxString = convertFloatToString(cx);
		std::string cyString = convertFloatToString(cy);
		std::string czString = convertFloatToString(cz);

		std::string toReturn = "Vertex 1 is (";
		toReturn.append(axString).append(", ").append(ayString).append(", ").append(azString).append(").  ");
		toReturn.append("Vertex 3 is (").append(bxString).append(", ").append(byString).append(", ").append(bzString).append(").  ");
		toReturn.append("Vertex 3 is (").append(cxString).append(", ").append(cyString).append(", ").append(czString).append(").\n");

		return toReturn;
	}
};

class WarpedTriangle : public Triangle {
public:
	Point v1, v2, v3;

	WarpedTriangle() : Triangle() {
	}

	WarpedTriangle(Point v1, Point v2, Point v3) : Triangle(v1, v2, v3) {

	}

	WarpedTriangle(float ax, float ay, float az,
			float bx, float by, float bz,
			float cx, float cy, float cz) : Triangle(ax, ay, az, bx, by, bz, cx, cy, cz) {
	}

	// If these two methods are called, then ray represents the transformed ray that has
	// been transformed into object coordinates; i.e. original sphere coordinates.
	bool intersect(Ray &ray, float* tHit, DifferentialGeometry* differentialGeometry) {
		return Triangle::intersect(ray, tHit, differentialGeometry);
	}

	bool intersectP(Ray &ray) {
		return Triangle::intersectP(ray);
	}

	std::string shapeType() {
		return "Warped Triangle";
	}

};

#endif /* SHAPE_H_ */
