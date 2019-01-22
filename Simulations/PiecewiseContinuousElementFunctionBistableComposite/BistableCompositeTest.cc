// -*- C++ -*-

#include "../../core/Definitions.h"
#include "BistableComposite.h"
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

typedef Elements::BistableComposite            Element;
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
  double precompression = 0.;
  vector<double> coefficients;
  //coefficients.push_back(0.44);
  //coefficients.push_back(5);
  //coefficients.push_back(12);
  //coefficients.push_back(0.05);
  
  coefficients.push_back(0.47253);
  coefficients.push_back(1.7776);
  coefficients.push_back(-5.2441);
  coefficients.push_back(-3.9768);
  coefficients.push_back(-8.6078);
  coefficients.push_back(1.1779);
  coefficients.push_back(1.7494);
  coefficients.push_back(-1.0041);
  coefficients.push_back(-1.8938);
  coefficients.push_back(5.3921);
  coefficients.push_back(4.5656);
  coefficients.push_back(-0.1124);
  
  const Element testElement(testElementNodes, coefficients,precompression);
  const double time = std::numeric_limits<double>::quiet_NaN();
  array<Vector, Element::NumberOfNodes> testDisplacements;
  testDisplacements[0](0) = .3;

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^d
  // *****************  </input parameters>  ********************
  // ============================================================

  // ============================================================
  // **********************  < solution>  ***********************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  const double correctEnergy = 3.606709e-01;
  typename Element::Forces correctForces;
  correctForces[0] <<  -3.303389e+00;
  typename Element::StiffnessMatrix correctStiffnessMatrix;
  correctStiffnessMatrix << -2.594986e+01;


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
