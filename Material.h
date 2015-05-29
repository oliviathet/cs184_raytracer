#ifndef MATERIAL_H_
#define MATERIAL_H_

// Class that represents the material of a surface.
// NOTE: A Primitive object has both a Material and a Shape
class Material {
	public:
		BRDFCoefficients constantBRDF;

	Material() {

	}

	BRDFCoefficients getBRDF(DifferentialGeometry& differentialGeometry, BRDFCoefficients* brdf) {
		return constantBRDF;
	}
};



#endif /* MATERIAL_H_ */
