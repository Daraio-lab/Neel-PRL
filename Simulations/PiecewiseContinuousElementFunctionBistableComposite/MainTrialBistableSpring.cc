// -*- C++ -*-

// Wave propagation in magnet + bistable spring case.

#include "BistableSpring.h"
#include "NonlinearMagnet.h"
#include "../../elements/PointMass.h"
#include "../../core/Definitions.h"
#include "../../core/Utilities.h"
#include "../../core/Assemblers.h"
#include "../../core/ImplicitDynamics.h"

typedef Elements::BistableSpring::SingleNodeBistableSpring					Bistable;
typedef Elements::BistableSpring::Properties								BistableProperties;
typedef Elements::NonlinearMagnet::Magnet									Magnet;
typedef Elements::NonlinearMagnet::Properties								MagnetProperties;
typedef Elements::PointMass::PointMassND<1,1>											Mass;
typedef Mass::Properties													MassProperties;
typedef Bistable::Vector													Vector;
typedef Bistable::Node														Node;
typedef Bistable::NodalDisplacements										NodalDisplacements;
typedef Assemblers::AssemblerThreePhysicalElements<Bistable, Magnet, Mass> 	Assembler;
typedef Dynamics::State														State;

int main() {

  /////////////////////////////////////////////////////////////
  // *************** DEFINING PARAMETERS ******************* //
  // VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV //
  size_t numberOfNodes = 150;
  double mass = 1;
  MassProperties massProperty(mass);
  double timestep = 0.01;
  double alpha = 0.001;  // damps out low frequency oscillations
  double beta = 0.001;   // damps out high frequency oscillations
  double newmarkBeta = 0.25;
  double newmarkGamma = 0.5;
  double amplitude = 1.8025;
  double exponent = -2.73;
  double d = 7;
  double precompression = 0;
  
  BistableProperties bistableProperties(2,1,10,0);  
  
  MagnetProperties magnetProperties(amplitude, exponent);

  /////////////////////////////////////////////////////////////
  // ******* GENERATING PROBLEM AND ASSIGNING NODES ******** //
  // VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV //
  
  Assembler problem(numberOfNodes+1);
  
  for (unsigned int i=0; i<numberOfNodes; i++) {
	
	Node node1;
	Node node2;
	
	node1._id = i;
	node1._position << i*d;
  
	node2._id = i+1;
	node2._position << (i+1)*d;
	
	array<Node, 2> magnetNodes;
	magnetNodes[0] = node1;
	magnetNodes[1] = node2;	
	
	Bistable element1(node1, bistableProperties, precompression);
	Magnet element2(magnetNodes, magnetProperties);
	Mass element3(node1, massProperty);
	
	problem.addElement(element1);
	problem.addElement(element2);
	problem.addElement(element3);  
  
  }
  
  /////////////////////////////////////////////////////////////
  // ****************** INITIAL CONDITIONS ***************** //
  // VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV //
  VectorXd current_u(numberOfNodes+1);
  VectorXd current_v(numberOfNodes+1);
  VectorXd current_a(numberOfNodes+1);
  
  current_u.fill(0);
  current_v.fill(0);
  current_a.fill(0);
  
  current_v(numberOfNodes) = 0;
  
  State currentState;
  
  currentState._displacements = current_u;
  currentState._velocities = current_v;
  currentState._accelerations = current_a;
  
  /////////////////////////////////////////////////////////////
  // ******************* OUTPUT FILE *********************** //
  // VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV //
  
  std::ofstream DisplacementsFile;
  DisplacementsFile.open("../../../output/Displacements.txt");
  
  /////////////////////////////////////////////////////////////
  // ****************** ITERATIVE SOLVING ****************** //
  // VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV //
  
  const unsigned int iterations = 1000;
  MatrixXd Disp(numberOfNodes+1, iterations);
  
  vector<Vector> nodalDisplacements;
  
  Dynamics::Newmark<Assembler> newmarkSolver(problem,
								  timestep,
								  alpha,
								  beta,
								  newmarkBeta,
								  newmarkGamma,
								  current_v,
								  iterations);
    
  for (unsigned int timeStepIndex = 0; timeStepIndex<iterations; timeStepIndex++) {
	if (timeStepIndex % (iterations/10)==0) {
	  printf("on iteration %5u / %5u (%%%5.1f) \n", timeStepIndex, iterations, 100.*timeStepIndex / float(iterations));	
	}
	
	vector<EssentialBoundaryCondition> essentialBCs;
	essentialBCs.push_back(EssentialBoundaryCondition(0,0,0));
	
	newmarkSolver.computeNewmarkUpdate(essentialBCs, timeStepIndex, Quiet);
	
	nodalDisplacements = newmarkSolver.getNodalDisplacements();
		
	for (unsigned int j = 0; j<numberOfNodes+1; j++) {
	  
	  Disp(j,timeStepIndex) = nodalDisplacements[j](0);
	  
	}
  
  }
  
  DisplacementsFile << Disp << endl;
  
  return 0;
  
}
