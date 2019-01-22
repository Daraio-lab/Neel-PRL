// -*- C++ -*-
#include "../../core/Definitions.h"
#include "NonlinearMagnet.h"
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

typedef Elements::NonlinearMagnet::Magnet            Element;
typedef Elements::NonlinearMagnet::Properties        ElementProperties;
typedef Element::Vector             				 				 Vector;
typedef Element::Point              				 				 Point;
typedef Element::Node               				 				 Node;
typedef Element::NodalDisplacements 				 				 NodalDisplacements;


int main() {

  // ============================================================
  // *****************  < input parameters>  ********************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  //Quadrature rule is for other elements. Not for a bar
  //const QuadratureRule<Element::SpatialDimension,
  //                     numberOfQuadraturePoints> quadratureRule =
  //  Quadrature::buildGaussianQuadrature<Element::SpatialDimension,
  //                                      numberOfQuadraturePoints>();

  const ElementProperties testElementProperties(4.204e-06, -2.73);
  array<Node, Element::NumberOfNodes> testElementNodes;
  Vector position1;
  position1(0) = 0.41;
  Vector position2;
  position2(0) = 1.32;
  testElementNodes[0] = Node(0, position1);
  testElementNodes[1] = Node(1, position2);
  
  const Element testElement(testElementNodes, testElementProperties);
  const double time = std::numeric_limits<double>::quiet_NaN();
  NodalDisplacements testDisplacements;
  testDisplacements[0](0) = 0.3;
  testDisplacements[1](0) = 0.6;

  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  // *****************  </input parameters>  ********************
  // ============================================================

  // ============================================================
  // **********************  < solution>  ***********************
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

  const double correctEnergy = 9.745069e-01;
  typename Element::Forces correctForces;
  correctForces[0] << -8.479417e-01;
  correctForces[1] << 8.479417e-01;
  typename Element::StiffnessMatrix correctStiffnessMatrix;
  correctStiffnessMatrix << 1.625683e+00, -1.625683e+00, -1.625683e+00, 1.625683e+00;





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
