#ifndef BRDFCOEFFICIENTS_H_
#define BRDFCOEFFICIENTS_H_

// This is a wrapper class for the BRDF coefficients.
// NOTE: a BRDFCoefficient is an instance variable of a Material
class BRDFCoefficients {
public:
	Color kd, ks, ka, kr;
	int sp;
	float ri;

	BRDFCoefficients() {
		this->sp = 1;
		this->ri = 0.0;
	}

	BRDFCoefficients(Color kd, Color ks, Color ka, Color kr, int sp, float ri) {
		this->kd = kd;
		this->ks = ks;
		this->ka = ka;
		this->kr = kr;
		this->sp = sp;
		this->ri = ri;
	}

	bool isRefractive() {
		return (this->ri > 0.0);
	}
};

#endif /* BRDFCOEFFICIENTS_H_ */
