// -*- C++ -*-

#include "../../core/Definitions.h"
#include "BistableCompositePiecewise.h"
#include "../../elements/ElementTests.h"

/*
  this test does two things:
  1. checks for consistency between the reported energy, forces, and
  stiffness matrix
  2. checks the values reported against known correct solutions for a fixed
  set of inputs
*/

/*
  instructions for making a test from this template:
  1. copy, paste, and rename this file.
  2. put in some random (but admissible) values for the node positions,
  element properties, material model parameters, time, etc.
  3. run the code with makingSolution = true.  this prints out a set of
  lines with the solution.
  4. copy and paste the solution into the area at the top of main, replacing
  what's there now.
  5. set makingSolution = false;
*/
const bool makingSolution = false;

//const unsigned int numberOfQuadraturePoints = 4; No quad points for bar

typedef Elements::BistableCompositePiecewise::SingleNodeBistableCompositePiecewise   Element;
typedef Elements::BistableCompositePiecewise::Properties		Properties;
typedef Element::Vector                        Vector;
typedef Element::Point                         Point;
typedef Element::Node                          Node;

int main() {

  // ============================================================
  // *****************  < input parameters>  ********************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  //Quadrature rule is for other elements. Not for a bar
  //const QuadratureRule<Element::SpatialDimension,
  //                     numberOfQuadraturePoints> quadratureRule =
  //  Quadrature::buildGaussianQuadrature<Element::SpatialDimension,
  //                                      numberOfQuadraturePoints>();

  Node testElementNodes;
  Vector position;
  position(0) = 0.41;
  testElementNodes = NodeWithId<Point>(0, position);
  vector<double> firstSplineCoefficients;
	vector<double> secondSplineCoefficients;
	vector<double> thirdSplineCoefficients;
	Matrix<double, 2, 1> splinePoints;
  
  firstSplineCoefficients.push_back(0.);
  firstSplineCoefficients.push_back(0.);
	firstSplineCoefficients.push_back(0.0593034);
	firstSplineCoefficients.push_back(-0.018119);
	
	thirdSplineCoefficients.push_back(-0.815088);
	thirdSplineCoefficients.push_back(0.886911);
	thirdSplineCoefficients.push_back(-0.271859);
	thirdSplineCoefficients.push_back(0.0237451);
	
	secondSplineCoefficients.push_back(0.0216472);
	secondSplineCoefficients.push_back(-0.0664308);
	secondSplineCoefficients.push_back(0.132853);
	secondSplineCoefficients.push_back(-0.0521944);
	secondSplineCoefficients.push_back(0.00531799);
	
	splinePoints(0,0) = 1.091;
	splinePoints(1,0) = 5.27035-1.454;
	
	const Properties properties(firstSplineCoefficients, secondSplineCoefficients, thirdSplineCoefficients, splinePoints);
	
  const Element testElement(testElementNodes, properties);
  
	const double time = std::numeric_limits<double>::quiet_NaN();
  array<Vector, Element::NumberOfNodes> testDisplacements;
  testDisplacements[0](0) = 0.3;

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^d
  // *****************  </input parameters>  ********************
  // ============================================================

  // ============================================================
  // **********************  < solution>  ***********************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  const double correctEnergy = 4.848093e-03;
  typename Element::Forces correctForces;
  correctForces[0] << 3.068991e-02;
  typename Element::StiffnessMatrix correctStiffnessMatrix;
  correctStiffnessMatrix << 8.599260e-02;

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // **********************  </Solution>  ***********************
  // ============================================================

  Elements::runElementTests(testElement,
                            testDisplacements,
                            time,
                            correctEnergy,
                            correctForces,
                            correctStiffnessMatrix,
                            makingSolution);

  return 0;
}
