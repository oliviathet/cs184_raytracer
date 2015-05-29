#ifndef VECTOR3_H_
#define VECTOR3_H_

#include <math.h>

class Vector3 {
  public:
    float x, y, z;

  Vector3() {

  }

  Vector3(float x, float y, float z) {
    this->x = x;
    this->y = y;
    this->z = z;
  }

  float dotProduct(Vector3 other) {
    return this->x * other.x + this->y * other.y + this->z * other.z;
  }

  Vector3 addVector(Vector3 other) {
    return Vector3(this->x + other.x, this->y + other.y, this->z + other.z);
  }

  Vector3 subtractVector(Vector3 other) {
    return Vector3(this->x - other.x, this->y - other.y, this->z - other.z);
  }

  Vector3 scaleVector(float c) {
    return Vector3(c*this->x, c*this->y, c*this->z);
  }

  Vector3 multiplyVector(Vector3 other) {
    return Vector3(this->x * other.x, this->y * other.y, this->z * other.z);
  }

  static Vector3 normalizeVector(Vector3 vector_to_normalize) {
    float magnitude = sqrt(pow(vector_to_normalize.x, 2) + pow(vector_to_normalize.y, 2) + pow(vector_to_normalize.z, 2));
    return Vector3(vector_to_normalize.x / magnitude, vector_to_normalize.y / magnitude, vector_to_normalize.z / magnitude);
  }

  Vector3 crossProduct(Vector3 other) {
    float i = this->y * other.z - this->z * other.y;
    float j = this->z * other.x - this->x * other.z;
    float k = this->x * other.y - this->y * other.x;
    return Vector3(i, j, k);
  }
};

#endif /* VECTOR3_H_ */
