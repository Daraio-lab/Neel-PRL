// -*- C++ -*-
#ifndef ELEMENT_NON_LINEAR_MAGNET
#define ELEMENT_NON_LINEAR_MAGNET

#include "../../core/Definitions.h"

namespace Elements {
namespace NonlinearMagnet {  
   
class Properties {
  public:
	 	double _amplitude;
	 	double _exponent;
	 	double _distance;
	 	Properties(const double amplitude, const double exponent) :
	  	  _amplitude(amplitude), _exponent(exponent) {
	 	}
   	};
   
class Magnet {
   
  public:
	  
	  static const unsigned int      NumberOfNodes = 2;
	  static const unsigned int      SpatialDimension = 1;
	  static const unsigned int      DegreesOfFreedom = 1;
	  
	  typedef NonlinearMagnet::Properties 				 			 Properties;
	  typedef Matrix<double, SpatialDimension, 1>        Point;
	  typedef NodeWithId<Point>                          Node;
	  typedef Matrix<double, DegreesOfFreedom, 1>        Vector;
	  typedef array<Vector, NumberOfNodes>               NodalDisplacements;
	  typedef Matrix<double, 
			NumberOfNodes*DegreesOfFreedom, 
			NumberOfNodes*DegreesOfFreedom> 			 					 StiffnessMatrix;
	  typedef Matrix<double, 
            NumberOfNodes*DegreesOfFreedom, 
            NumberOfNodes*DegreesOfFreedom>            MassMatrix;
	  typedef array<Vector, 2>                           Forces;
	  
	  Magnet(const array<Node, 2> & nodes, const Properties & properties) :
		_properties(properties) {
	  
		_nodeIds[0] = nodes[0]._id;
		_nodeIds[1] = nodes[1]._id;
		
		_distance = std::abs(nodes[1]._position(0)-nodes[0]._position(0));
	  
	   _massMatrix.fill(0.);
		
	  }
	  
	  double
	  computeEnergy(const NodalDisplacements & displacements, const double time) const {
	    ignoreUnusedVariables(time);
		
			double energy = 0.;
		
			double u = displacements[1](0)-displacements[0](0);
			
			energy = (-_properties._amplitude/(_properties._exponent+1.)*pow(_distance+u, _properties._exponent+1.) + _properties._amplitude*pow(_distance, _properties._exponent)*u);
			
			return energy;
	  
	  }
	  
	  Forces
	  computeForces(const NodalDisplacements & displacements, const double time) const {
			ignoreUnusedVariables(time);
		
			Forces forces;
			forces[0].fill(0.);
			forces[1].fill(0.);
			
			double force = 0.;
		
			double u = displacements[1](0) - displacements[0](0);
			
			force = (_properties._amplitude*pow(_distance+u, _properties._exponent) - _properties._amplitude*pow(_distance, _properties._exponent));
			
			forces[0](0) = force;
			forces[1](0) = -force; 
			
			return forces;	  
	  }
	  
	  StiffnessMatrix
	  computeStiffnessMatrix(const NodalDisplacements & displacements, const double time)
	  const {
	  
	    ignoreUnusedVariables(time);
	  
			StiffnessMatrix stiffnessMatrix;
			stiffnessMatrix.fill(0.);
			
			double stiffness = 0.;
		
			double u = displacements[1](0) - displacements[0](0);
			
			stiffness = (-_properties._exponent*_properties._amplitude*pow(_distance+u, _properties._exponent-1.));
						
			stiffnessMatrix(0,0) = stiffness;
			stiffnessMatrix(0,1) = -stiffness;
			stiffnessMatrix(1,0) = -stiffness;
			stiffnessMatrix(1,1) = stiffness;			
			
			return stiffnessMatrix;
	  }
	  
	  array<size_t, 2>
    getNodeIds() const {
      return _nodeIds;
    }
	  
	  MassMatrix
    computeConsistentMassMatrix() const {
      return _massMatrix;
    }
	  
    MassMatrix
    computeLumpedMassMatrix() const {
			return _massMatrix;
    }  
	  
	private:
	
	  Properties _properties;
    array<size_t, 2> _nodeIds;
	  MassMatrix _massMatrix;
	  double _distance;
	  array<double, 3> _coefficients;
	};
	
  }
}

#endif //ELEMENT_NON_LINEAR_MAGNET
