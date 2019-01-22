// -*- C++ -*-

// Wave propagation in magnet + bistable composite case. Data used from Andres's experiments. Refer to: Dropbox/ae108repo/output
// Using splines to fit the experimental data using 4 points. Refer to BistableCurveFit.nb for the coefficients. Adding some quantity of randomness in the simulations. Eqpoints in each case are determined by the MATLAB program FindEqmPoint.m. Simulations for PRL paper

#include "BistableCompositePiecewise.h"
#include "NonlinearMagnet.h"
#include "../../elements/PointMass.h"
#include "../../core/Definitions.h"
#include "../../core/Utilities.h"
#include "../../core/Assemblers.h"
#include "ImplicitDynamics.h"
#include <ctime>

typedef typename Elements::BistableCompositePiecewise::SingleNodeBistableCompositePiecewise					Bistable;
typedef typename Elements::BistableCompositePiecewise::Properties							BistableProperties;
typedef typename Elements::NonlinearMagnet::Magnet														Magnet;
typedef typename Elements::NonlinearMagnet::Properties												MagnetProperties;
typedef typename Elements::PointMass::PointMassND<1,1>												Mass;
typedef typename Mass::Properties																							MassProperties;
typedef typename Bistable::Vector																							Vector;
typedef typename Bistable::Node																								Node;
typedef typename Bistable::NodalDisplacements																	NodalDisplacements;
typedef typename Assemblers::AssemblerThreePhysicalElements<Bistable, 
																														Magnet, 
																														Mass>						 	Assembler;
typedef typename Dynamics::State																							State;

int main(int argc, char *argv[]) {
	
	ignoreUnusedVariable(argc);
	
  /////////////////////////////////////////////////////////////
  // *************** DEFINING PARAMETERS ******************* //
  // VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV //
  
  /*
  	The unit system used is [gm cm ms]
  	Unit of length: cm
  	Unit of mass: gm
  	Unit of time: ms
 		Unit of force: 10 N
  	Unit of energy: 10^-1 J
  */
  
  size_t numberOfNodes = 25;
  double mass = 21.;
  MassProperties massProperty(mass);
  double timestep = 0.25;
  double alpha = atof(argv[3]);  // damps out low frequency oscillations
  double beta = 0.;   // damps out high frequency oscillations
  double newmarkBeta = 0.25;
  double newmarkGamma = 0.5;
  double amplitude = 17.51805;
  double exponent = -3.2744;
  double counter = atof(argv[1]);
  double d = atof(argv[2]); // lattice distance
  
  /////////////////////////////////////////////////////////////
  // ******* GENERATING PROBLEM AND ASSIGNING NODES ******** //
  // VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV //
  
  Assembler problem(numberOfNodes+1);
  
  vector<double> firstSplineCoefficients;
	vector<double> secondSplineCoefficients;
	vector<double> thirdSplineCoefficients;
	Matrix<double, 2, 1> splinePoints;
  
  if(counter==1){
		firstSplineCoefficients.push_back(0.);
		firstSplineCoefficients.push_back(0.);
		firstSplineCoefficients.push_back(0.0256743);
		firstSplineCoefficients.push_back(-0.00489034);
	
		thirdSplineCoefficients.push_back(-0.809805);
		thirdSplineCoefficients.push_back(4.92206);
		thirdSplineCoefficients.push_back(-0.979682);
		thirdSplineCoefficients.push_back(0.0624519);
	
		secondSplineCoefficients.push_back(0.10991);
		secondSplineCoefficients.push_back(-0.204049);
		secondSplineCoefficients.push_back(0.160139);
		secondSplineCoefficients.push_back(-0.040711);
		secondSplineCoefficients.push_back(0.00291668);
		
		splinePoints(0,0) = 1.75;
		splinePoints(1,0) = 5.465-0.236;
	}
  
  if(counter==2){
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
  }
  
  if(counter==3){
  	
		firstSplineCoefficients.push_back(0.);
		firstSplineCoefficients.push_back(0.);
		firstSplineCoefficients.push_back(0.0402355);
		firstSplineCoefficients.push_back(-0.0121484);
	
		thirdSplineCoefficients.push_back(-0.142865);
		thirdSplineCoefficients.push_back(0.409241);
		thirdSplineCoefficients.push_back(-0.213546);
		thirdSplineCoefficients.push_back(0.0272003);
	
		secondSplineCoefficients.push_back(0.139316);
		secondSplineCoefficients.push_back(-0.416266);
		secondSplineCoefficients.push_back(0.48557);
		secondSplineCoefficients.push_back(-0.208458);
		secondSplineCoefficients.push_back(0.0280114);
		
		splinePoints(0,0) = 1.104;
		splinePoints(1,0) = 3.97096-1.354;
  }
    
	BistableProperties	bistableProperties(firstSplineCoefficients, secondSplineCoefficients, thirdSplineCoefficients, splinePoints);
		
  for (unsigned int i=0; i<numberOfNodes; i++) {
  
  	Node node1;
		Node node2;
	
		node1._id = i;
		node1._position << i*d;
  
		node2._id = i+1;
		node2._position << (i+1.)*d;
	
		array<Node, 2> magnetNodes;
		magnetNodes[0] = node1;
		magnetNodes[1] = node2;	
				
		MagnetProperties magnetProperties(amplitude, exponent);
	
		Bistable element1(node2, bistableProperties);
		
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
  
  current_u.fill(0.);
  current_v.fill(0.);
  current_a.fill(0.);
  
  // current_v(0) = 0.6;
    
  State currentState;
  currentState._displacements = current_u;
  currentState._velocities = current_v;
  currentState._accelerations = current_a;
  
  /////////////////////////////////////////////////////////////
  // ******************* OUTPUT FILE *********************** //
  // VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV //
    
  string folder;
  string str1 = "../../../output/PRL2015aAndresExperimentsCorrected/";
  string str2 = "ExperimentalSplinedBistableMagneticChainDisplacementsLattice";
  string str3 = argv[2];
  string str4 = "Rail";
  string str5;
  if (counter==1.)
  {  str5 = "21p5";}
  if (counter==2.)
  {  str5 = "22";}
  if (counter==3.)
  {  str5 = "22p5";}
  string str6 = ".txt";
  
  folder.append(str1);  
  folder.append(str2);
  folder.append(str3);
  folder.append(str4);
  folder.append(str5);
  folder.append(str6);
    
  std::ofstream DisplacementsFile; 
  DisplacementsFile.open(folder);
  
  /////////////////////////////////////////////////////////////
  // ****************** ITERATIVE SOLVING ****************** //
  // VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV //
  
  const unsigned int iterations = 10000;
  MatrixXd Disp(numberOfNodes+1, iterations);
  Disp.fill(0.);
  
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
	
		essentialBCs.push_back(EssentialBoundaryCondition(0,0,atof(argv[4])));
		essentialBCs.push_back(EssentialBoundaryCondition(numberOfNodes,0,0.));
			
		newmarkSolver.computeNewmarkUpdate(essentialBCs, timeStepIndex, Quiet);
	
		nodalDisplacements = newmarkSolver.getNodalDisplacements();
		
		for (unsigned int j = 0; j<numberOfNodes+1; j++) {
	  
		  Disp(j,timeStepIndex) = nodalDisplacements[j](0);
	  
		}
  	cout << endl;
  }
  
  DisplacementsFile << Disp << endl;
  
  return 0;
  
}