#ifndef CAMERA_H_
#define CAMERA_H_

static void printRay(Ray ray);

static void printSample(Sample sample);

static void printVector(Vector3 vector);

class Camera {
	public:
		// Note: Camera needs Film width and height to correctly compute viewing ray
		int viewingPlaneWidth, viewingPlaneHeight;

		// z-coordinates must be the same
		// Left coordinates must align vertically
		// Right coordinates must align vertically
		// Top coordinates must align horizontally
		// Bottom coordinates must align horizontally
		Point eye;
		Point imagePlaneTopLeft;
		Point imagePlaneBottomLeft;
		Point imagePlaneTopRight;
		Point imagePlaneBottomRight;


		Camera() {

		}

		Camera(int widthOfFilm, int heightOfFilm) {
			this->viewingPlaneWidth = widthOfFilm;
			this->viewingPlaneHeight = heightOfFilm;
		}

		Camera(int widthOfFilm, int heightOfFilm,
				float ex, float ey, float ez,
				float llx, float lly, float llz,
				float lrx, float lry, float lrz,
				float urx, float ury, float urz,
				float ulx, float uly, float ulz) {
			this->viewingPlaneWidth = widthOfFilm;
			this->viewingPlaneHeight = heightOfFilm;
			this->eye = Point(ex, ey, ez);
			this->imagePlaneTopLeft = Point(ulx, uly, ulz);
			this->imagePlaneBottomLeft = Point(llx, lly, llz);
			this->imagePlaneTopRight = Point(urx, ury, urz);
			this->imagePlaneBottomRight = Point(lrx, lry, lrz);
		}

		// Given a sample in FILM coordinates, this method generates a ray from the eye (0, 0, 0)
		// to the sample in IMAGE PLANE coordinates (i.e. [-1, 1])
		void generateRay(Sample& sample, Ray* ray, bool distributedRayTracing) {

			// Textbook page 75
			ray->position = eye;

			float rectangleWidth = imagePlaneTopRight.x - imagePlaneTopLeft.x;
			float rectangleHeight = imagePlaneTopRight.y - imagePlaneBottomRight.y;

			float horizontalSampleDistance = rectangleWidth / viewingPlaneWidth;
			float verticalSampleDistance = rectangleHeight / viewingPlaneHeight;

			float imagePlaneX, imagePlaneY;

			if (distributedRayTracing) {
				// Random number between 0 and 99
				int randomNumberX = rand() % 100;
				int randomNumberY = rand() % 100;

				imagePlaneX = imagePlaneTopLeft.x + ((rectangleWidth * (sample.x + (randomNumberX / 100.0))) / viewingPlaneWidth);
				imagePlaneY = imagePlaneBottomLeft.y + ((rectangleHeight * (sample.y + (randomNumberY / 100.0))) / viewingPlaneHeight);
			} else {
				imagePlaneX = imagePlaneTopLeft.x + ((rectangleWidth * (sample.x + 0.5)) / viewingPlaneWidth);
				imagePlaneY = imagePlaneBottomLeft.y + ((rectangleHeight * (sample.y + 0.5)) / viewingPlaneHeight);
			}

			ray->direction = Vector3(imagePlaneX - eye.x, imagePlaneY - eye.y, imagePlaneTopLeft.z - eye.z);

			ray->t_min = 0.001;
			ray->t_max = FLT_MAX;
		}
};

#endif /* CAMERA_H_ */
