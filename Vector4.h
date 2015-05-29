#ifndef VECTOR4_H_
#define VECTOR4_H_

#include <math.h>

class Vector4 {
  public:
    float x, y, z, w;

  Vector4() {

  }

  Vector4(float x, float y, float z, float w) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
  }

  float dotProduct(Vector4 other) {
    return this->x * other.x + this->y * other.y + this->z * other.z + this->w * other.w;
  }

  Vector4 addVector(Vector4 other) {
    return Vector4(this->x + other.x, this->y + other.y, this->z + other.z, this->w + other.w);
  }

  Vector4 subtractVector(Vector4 other) {
    return Vector4(this->x - other.x, this->y - other.y, this->z - other.z, this->w - other.w);
  }

  Vector4 scaleVector(float c) {
    return Vector4(c*this->x, c*this->y, c*this->z, c*this->w);
  }

  Vector4 multiplyVector(Vector4 other) {
    return Vector4(this->x * other.x, this->y * other.y, this->z * other.z, this->w * other.w);
  }
};

#endif /* VECTOR4_H_ */
