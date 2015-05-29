#ifndef AGGREGATEPRIMITIVE_H_
#define AGGREGATEPRIMITIVE_H_

static void printColor(Color color);
static void printRay(Ray ray);

class AggregatePrimitive {
	public:
		std::vector<GeometricPrimitive*> listOfPrimitives;

	AggregatePrimitive() {

	}

	AggregatePrimitive(std::vector<GeometricPrimitive*> listOfPrimitives) {
		this->listOfPrimitives = listOfPrimitives;
	}

	void addPrimitive(GeometricPrimitive* primitive) {
		listOfPrimitives.push_back(primitive);
	}

	// Just loops through all the primitives in the list and calls the intersect routine
	// Compare tHit of all the intersections and return that of the nearest one (we want the first hit)
	// NOTE: in->primitive should be set to the pointer to that primitive
	bool intersect(Ray& ray, float* tHit, Intersection* in) {

		bool hit = false;
		float smallestTSeenSoFar = FLT_MAX;
		float* pointerOfSmallestTSeenSoFar;
		Intersection* pointerOfClosestIntersectionSoFar;
		bool foundAnIntersection = false;

		// For each primitive...
		for (std::vector<GeometricPrimitive*>::size_type i = 0; i < listOfPrimitives.size(); i++) {
			if (foundAnIntersection) {
				// Find whether the current primitive hits the ray...
				// ONLY populate 'in' if the tHit of the current primitive is smaller than the last argument passed in
				// In this situation, we want our tHit of the current primitive to be smaller than the previous tHit,
				// which by logic below, ends up being the smallest tHit seen so far
				if (listOfPrimitives[i]->intersectWithMaxT(ray, tHit, in, *tHit)) {
					hit = true;
					// ... if it does, then update points for "lowest" t-value and corresponding intersection object
					if (*tHit < smallestTSeenSoFar) {
						smallestTSeenSoFar = *tHit;
						pointerOfSmallestTSeenSoFar = tHit;
						pointerOfClosestIntersectionSoFar = in;
					} else {
						// since listOfPrimitives[i]->intersect(ray, tHit, in) actually POPULATES tHit and in,
						// we need to reset the pointers to our currently "smallest T" and "primitive of smallest T"
						// if the currently iterated primitive was NOT the smallest T seen so far

						// However, since an Intersection isn't a primitive object, you can't maintain an
						// Intersection closestIntersectionSoFar.
						// Instead, one way is to separate the logic into the Primitive class.
						// (see above comment about intersectWithMaxT)
						pointerOfSmallestTSeenSoFar = &smallestTSeenSoFar;
					}
				}
			} else {
				// Find whether the current primitive hits the ray...
				if (listOfPrimitives[i]->intersect(ray, tHit, in)) {
					hit = true;
					// ... if it does, then update points for "lowest" t-value and corresponding intersection object
					if (*tHit < smallestTSeenSoFar) {
						smallestTSeenSoFar = *tHit;
						pointerOfSmallestTSeenSoFar = tHit;

						pointerOfClosestIntersectionSoFar = in;
					} else {
						pointerOfSmallestTSeenSoFar = &smallestTSeenSoFar;
					}
					foundAnIntersection = true;
				}
			}
		}


		// After we've checked all of our primitives, update our 'tHit' and 'in' to the nearest combination
		*tHit = *pointerOfSmallestTSeenSoFar;
		in = pointerOfClosestIntersectionSoFar;

		// DEBUGGING
		Intersection testIntersection = *pointerOfClosestIntersectionSoFar;
		GeometricPrimitive* testPrimitive = testIntersection.primitive;

		// Return whether our ray hit any of our primitives
		return hit;
	}

	// Just loops through all the primitives in the list and calls the intersectP routine
	// Returns true if our input ray hits ANY primitive, and false otherwise
	bool intersectP(Ray& ray) {
		for (std::vector<GeometricPrimitive*>::size_type i = 0; i < listOfPrimitives.size(); i++) {
			if (listOfPrimitives[i]->intersectP(ray)) {
				return true;
			}
		}
		return false;
	}

	// Just loops through all the primitives in the list and calls the intersectP routine
	// Returns true if our input ray hits ANY primitive, and false otherwise
	// NOTE: Does not consider primitives in our aggregatePrimitive if they match 'sourcePrimitive'
	bool intersectP(Ray& ray, GeometricPrimitive *sourcePrimitive) {
		for (std::vector<GeometricPrimitive*>::size_type i = 0; i < listOfPrimitives.size(); i++) {
			if (listOfPrimitives[i] == sourcePrimitive) {
				continue;
			}
			if (listOfPrimitives[i]->intersectP(ray)) {
				return true;
			}
		}
		return false;
	}

	// NOTE: This will never get called
	void getBRDF(DifferentialGeometry& differentialGeometry, BRDFCoefficients* brdf) {
		exit(1);
	}
};

#endif /* AGGREGATEPRIMITIVE_H_ */
