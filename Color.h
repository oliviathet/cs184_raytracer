#ifndef COLOR_H_
#define COLOR_H_

#include <math.h>

class Color {
  public:
    float r, g, b;

  Color() {
	  this->r = 0.0;
	  this->g = 0.0;
	  this->b = 0.0;
  }

  Color(float r, float g, float b) {
    this->r = r;
    this->g = g;
    this->b = b;
  }

  Color addColor(Color other) {
    return Color(this->r + other.r, this->g + other.g, this->b + other.b);
  }

  Color subtractColor(Color other) {
    return Color(this->r - other.r, this->g - other.g, this->b - other.b);
  }

  Color scaleColor(float c) {
    return Color(c*this->r, c*this->g, c*this->b);
  }

  Color multiplyColor(Color other) {
    return Color(this->r * other.r, this->g * other.g, this->b * other.b);
  }

  bool greaterThanZero() {
	  return (this->r > 0.0 || this->g > 0.0 || this->b > 0.0);
  }
};

#endif /* COLOR_H_ */
