#ifndef RAY_H_
#define RAY_H_

#include "Point.h"
#include "Vector3.h"

class Ray {
  public:
    Point position;
    Vector3 direction;
    float t_min, t_max;

  Ray() {

  }

  Ray(Point position, Vector3 direction, float t_min, float t_max) {
	  this->position = position;
	  this->direction = direction;
	  this->t_min = t_min;
	  this->t_max = t_max;
  }

};

#endif /* RAY_H_ */
