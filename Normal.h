#ifndef NORMAL_H_
#define NORMAL_H_

#include <math.h>

// Class that represents a normal vector
class Normal {
  public:
    float x, y, z;

  Normal() {

  }

  Normal(float x, float y, float z) {
    this->x = x;
    this->y = y;
    this->z = z;
  }

  Normal addNormal(Normal other) {
    return Normal(this->x + other.x, this->y + other.y, this->z + other.z);
  }

  Normal subtractNormal(Normal other) {
    return Normal(this->x - other.x, this->y - other.y, this->z - other.z);
  }


  static Normal normalizeNormal(Normal other) {
    float normalizedMagnitude = sqrt(pow(other.x, 2.0) + pow(other.y, 2.0) + pow(other.z, 2.0));
    return Normal(other.x / normalizedMagnitude, other.y / normalizedMagnitude, other.z / normalizedMagnitude);
  }
};

#endif /* NORMAL_H_ */
