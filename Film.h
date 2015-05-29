#ifndef FILM_H_
#define FILM_H_

#include "lodepng.h"


class Bucket {
public:
	int numberOfSamples;
	Color color;
	Bucket() {
		Color initialColor(0.0, 0.0, 0.0);
		color = initialColor;
		numberOfSamples = 0;
	}

	void add(Color colorToAdd) {
		// TODO: Need to define this
		color = color.addColor(colorToAdd);
		numberOfSamples++;

	}
};


// Class that represents the output of our ray tracing
// The output is represented by a 2D array of Colors.
class Film {

	static void printSample(Sample sample) {
		std::cout << "Sample: x = " << sample.x << "; y = " << sample.y << "\n";
	}

	static void printColor(Color color) {
		std::cout << "Color: r = " << color.r << "; g = " << color.g << "; b = " << color.b << "\n";
	}

	public:
		std::vector<Bucket> buckets;
		int width, height;
		std::vector<unsigned char> rawImage;

	Film() {

	}

	Film(int width, int height) {
		this->width = width;
		this->height = height;

		// Might fail...
		buckets.resize(width * height);

		// Initialize one empty bucket per pixel in our image plane
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				buckets.push_back(Bucket());

			}
		}
	}

	void commitColor(Sample& sample, Color& color) {
		color.r = color.r * 255.0;
		color.g = color.g * 255.0;
		color.b = color.b * 255.0;

		// Find the index of our 'buckets' array that corresponds with the sample's pixel, and add the 'color' to this index of the array
		// NOTE: sample.x and sample.y are in the coordinate system of the viewing plane, i.e. [width -> height]
		buckets[width * sample.y + sample.x].add(color);

	}


	void convertToRawData() {
		rawImage.resize(4 * width * height);
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				double numSamplesInBucket = buckets[width * y + x].numberOfSamples;

				rawImage[4 * width * y + 4 * x + 0] = fmin((buckets[width * (height - y) + x].color.r / numSamplesInBucket), 255.0);
				rawImage[4 * width * y + 4 * x + 1] = fmin((buckets[width * (height - y) + x].color.g / numSamplesInBucket), 255.0);
				rawImage[4 * width * y + 4 * x + 2] = fmin((buckets[width * (height - y) + x].color.b / numSamplesInBucket), 255.0);
				rawImage[4 * width * y + 4 * x + 3] = 255.0;
			}
		}

	}

	// First converts our 2D array of color into raw data so LodePNG can use it
	// Then encodes the image to a file with LodePNG
	void writeImage(const char* filename) {

		// Generate rawImage from our vector of buckets
		convertToRawData();

		// Encode the image so LodePNG can play with it
		unsigned error = lodepng::encode(filename, rawImage, width, height);

		// If there's an error, display it
	    if (error) {
	    	std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
	    }
	    std::cout << "Done writing.";
	}
};

#endif /* FILM_H_ */
