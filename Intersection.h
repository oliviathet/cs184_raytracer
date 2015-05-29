#ifndef INTERSECTION_H_
#define INTERSECTION_H_

#include "Primitive.h"

// Forward declaration
// http://stackoverflow.com/questions/28925071/c-confusion-with-forward-declaring-classes-and-using-include-statements
class GeometricPrimitive;

// This class represents an intersection point of a ray and a primitive.
class Intersection {
	public:
		DifferentialGeometry differentialGeometry;
		GeometricPrimitive *primitive;

		Intersection() {

		}

		Intersection(DifferentialGeometry differentialGeometry, GeometricPrimitive *primitive) {
			this->differentialGeometry = differentialGeometry;
			this->primitive = primitive;
		}

};

#endif /* INTERSECTION_H_ */
