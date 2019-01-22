// -*- C++ -*-
#ifndef ELEMENT_BISTABLE_COMPOSITE
#define ELEMENT_BISTABLE_COMPOSITE

#include "../core/Definitions.h"

namespace Elements {
	
class BistableComposite {

  public:

  static const unsigned int     NumberOfNodes = 1;
  static const unsigned int     SpatialDimension = 1;
  static const unsigned int     DegreesOfFreedom = 1;

  typedef Matrix<double, SpatialDimension, 1>        		   		Point;
  typedef NodeWithId<Point>                          		   		Node;
  typedef Matrix<double, DegreesOfFreedom, 1>        		   		Vector;
  typedef array<Vector, 1>                               	   	NodalDisplacements;
  typedef Matrix<double, DegreesOfFreedom, DegreesOfFreedom>  StiffnessMatrix;
	typedef Matrix<double, DegreesOfFreedom, DegreesOfFreedom>  MassMatrix;
	typedef array<Vector, 1>                           		   		Forces;

  BistableComposite(const Node & node, const vector<double> coefficients, const double precompression) :
  _nodeId(node._id), _coefficients(coefficients), _precompression(precompression) {
  	 _massMatrix.fill(0);
  	 //for(unsigned int m = 0; m < _coefficients.size(); m++){
			//cout <<  _coefficients[m] << endl;
		//}
	}

  double
  computeEnergy(const NodalDisplacements & displacements, const double time) const {
	 	ignoreUnusedVariables(time);
		
    double energy = 0.;
		double u = displacements[0](0)+_precompression;
		
		for(unsigned int m = 0; m < _coefficients.size(); m++){
			energy += _coefficients[m]*pow(u,m);
		}
		
    return energy;
	}

  Forces
  computeForces(const NodalDisplacements & displacements, const double time) const {
    ignoreUnusedVariables(time);
		
		Forces forces;
 		forces[0].fill(0.);
		
 		double u = displacements[0](0)+_precompression;
		double force = 0.;
		
		for(unsigned int m = 1; m < _coefficients.size(); m++){
			force += m*_coefficients[m]*pow(u,m-1);
		}
		
		forces[0](0) = force;
		
    return forces;
	}

  StiffnessMatrix
  computeStiffnessMatrix(const NodalDisplacements & displacements, const double time) const {
    ignoreUnusedVariables(time);
		
		StiffnessMatrix stiffnessMatrix;
    stiffnessMatrix.fill(0.);
		
		double u = displacements[0](0)+_precompression;
		double stiffness = 0.;
		
		for(unsigned int m = 2; m < _coefficients.size(); m++){
			stiffness += m*(m-1)*_coefficients[m]*pow(u,m-2);
		}
		
		stiffnessMatrix(0,0) = stiffness;
		return stiffnessMatrix;
  }

  array<size_t, 1>
  getNodeIds() const {
    return (array<size_t, 1>) {{_nodeId}};
  }
	  
	MassMatrix
  computeLumpedMassMatrix() const{
		return _massMatrix;
	}
	  
  MassMatrix
  computeConsistentMassMatrix() const{
		return _massMatrix;
  }
	  
  private:

  	size_t _nodeId;
   	vector<double> _coefficients; 
  	MassMatrix _massMatrix;
  	double _precompression;
  };

}


#endif //ELEMENT_BISTABLE_COMPOSITE
