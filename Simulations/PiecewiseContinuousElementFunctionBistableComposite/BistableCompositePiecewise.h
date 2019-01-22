// -*- C++ -*-

#ifndef ELEMENT_BISTABLE_COMPOSITE_PIECEWISE
#define ELEMENT_BISTABLE_COMPOSITE_PIECEWISE

#include "../core/Definitions.h"

namespace Elements {
namespace BistableCompositePiecewise{

class Properties {
public:
  vector<double> _firstSplineCoefficients;
  vector<double> _secondSplineCoefficients;
  vector<double> _thirdSplineCoefficients;
  Matrix<double, 2, 1> _splinePoints;
  Properties(const vector<double> firstSplineCoefficients, const vector<double> secondSplineCoefficients, const vector<double> thirdSplineCoefficients, const Matrix<double, 2, 1> splinePoints) :
    _firstSplineCoefficients(firstSplineCoefficients), _secondSplineCoefficients(secondSplineCoefficients), _thirdSplineCoefficients(thirdSplineCoefficients),_splinePoints(splinePoints) {
  }
};

//template<template <class> class NodeType = NodeWithId>
class SingleNodeBistableCompositePiecewise {

public:

  static const unsigned int      NumberOfNodes = 1;
  static const unsigned int      SpatialDimension = 1;
  static const unsigned int      DegreesOfFreedom = 1;

  typedef BistableCompositePiecewise::Properties     Properties;
  typedef Matrix<double, SpatialDimension, 1>        Point;
  typedef NodeWithId<Point>                          Node;
  typedef Matrix<double, DegreesOfFreedom, 1>        Vector;
  typedef array<Vector, 1>                           NodalDisplacements;
  typedef Matrix<double, DegreesOfFreedom, DegreesOfFreedom> StiffnessMatrix;
  typedef Matrix<double,
                 NumberOfNodes*DegreesOfFreedom,
                 NumberOfNodes*DegreesOfFreedom>                  MassMatrix;
  typedef array<Vector, 1>                           Forces;

  SingleNodeBistableCompositePiecewise(const Node & node, const Properties & properties) :
    _nodeId(node._id), _properties(properties) {

    Matrix<double, 2, 1> pointsOfSplineChange = _properties._splinePoints;
		_u1 = pointsOfSplineChange(0,0);
		_u2 = pointsOfSplineChange(1,0);
    _massMatrix.fill(0.);

  }

  double
  computeEnergy(const NodalDisplacements & displacements, const double time) const {
    ignoreUnusedVariables(time);

    double energy = 0;
		double u = displacements[0](0);
	
		if (u <= _u1)
		{
			for (unsigned int i = 0; i < _properties._firstSplineCoefficients.size(); i++){
				energy += _properties._firstSplineCoefficients[i]*pow(u,i);
			}
		}
	
		if (u <= _u2 && u > _u1)
		{
			for (unsigned int i = 0; i < _properties._secondSplineCoefficients.size(); i++){
				energy += _properties._secondSplineCoefficients[i]*pow(u,i);
			}
		}
	
		if (u > _u2)
		{
			for (unsigned int i = 0; i < _properties._thirdSplineCoefficients.size(); i++){
				energy += _properties._thirdSplineCoefficients[i]*pow(u,i);
			}
		}
	
    return energy;

  }

  Forces
  computeForces(const NodalDisplacements & displacements, const double time) const {
    ignoreUnusedVariables(time);

    Forces forces;
    forces[0].fill(0);
		
		double force = 0;
		double u = displacements[0](0);
		
		if (u <= _u1)
		{
			for (unsigned int i = 1; i < _properties._firstSplineCoefficients.size(); i++){
				force += i*_properties._firstSplineCoefficients[i]*pow(u,i-1);
			}
		}
	
		if (u <= _u2 && u > _u1)
		{
			for (unsigned int i = 1; i < _properties._secondSplineCoefficients.size(); i++){
				force += i*_properties._secondSplineCoefficients[i]*pow(u,i-1);
			}
		}
	
		if (u > _u2)
		{
			for (unsigned int i = 1; i < _properties._thirdSplineCoefficients.size(); i++){
				force += i*_properties._thirdSplineCoefficients[i]*pow(u,i-1);
			}	
		}
		
		forces[0](0) = force;

    return forces;
  }

  StiffnessMatrix
  computeStiffnessMatrix(const NodalDisplacements & displacements, const double time)
    const {
    ignoreUnusedVariables(time);

    StiffnessMatrix stiffnessMatrix;
    stiffnessMatrix.fill(0);
		
		double stiffness = 0;
		double u = displacements[0](0);
		
		if (u <= _u1)
		{
			for (unsigned int i = 2; i < _properties._firstSplineCoefficients.size(); i++){
				stiffness += i*(i-1)*_properties._firstSplineCoefficients[i]*pow(u,i-2);
			}
		}
	
		if (u <= _u2 && u > _u1)
		{
			for (unsigned int i = 2; i < _properties._secondSplineCoefficients.size(); i++){
				stiffness += i*(i-1)*_properties._secondSplineCoefficients[i]*pow(u,i-2);
			}
		}
	
		if (u > _u2)
		{
			for (unsigned int i = 2; i < _properties._thirdSplineCoefficients.size(); i++){
				stiffness += i*(i-1)*_properties._thirdSplineCoefficients[i]*pow(u,i-2);
			}
		}
		
		stiffnessMatrix.fill(stiffness);

    return stiffnessMatrix;
  }

  array<size_t, 1>
  getNodeIds() const {
    return (array<size_t, 1>) {{_nodeId}};
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
  size_t _nodeId;
	double _u1;
	double _u2;
  Properties _properties;
  MassMatrix _massMatrix;

};

}
}

#endif //ELEMENT_BISTABLE_COMPOSITE_PIECEWISE