// -*- C++ -*-
#ifndef ELEMENT_BISTABLE_SPRING
#define ELEMENT_BISTABLE_SPRING

#include "../core/Definitions.h"

namespace Elements {
namespace BistableSpring{

class Properties {
public:
  double _delta;
  double _lmeta;
  double _stiffness;
  double _density;
  Properties(const double delta, const double lmeta, const double stiffness, const double density) :
    _delta(delta), _lmeta(lmeta), _stiffness(stiffness), _density(density) {
  }
};

template<template <class> class NodeType = NodeWithId>
class SingleNodeBistableSpring {

public:

  static const unsigned int      NumberOfNodes = 1;
  static const unsigned int      SpatialDimension = 1;
  static const unsigned int      DegreesOfFreedom = 1;

  typedef BistableSpring::Properties                 Properties;
  typedef Matrix<double, SpatialDimension, 1>        Point;
  typedef NodeType<Point>                            Node;
  typedef Matrix<double, DegreesOfFreedom, 1>        Vector;
  typedef array<Vector, 1>                           NodalDisplacements;
  typedef Matrix<double, DegreesOfFreedom, DegreesOfFreedom> StiffnessMatrix;
  typedef Matrix<double,
                 NumberOfNodes*DegreesOfFreedom,
                 NumberOfNodes*DegreesOfFreedom>                  MassMatrix;
  typedef array<Vector, 1>                           Forces;


  SingleNodeBistableSpring(const Node & node, const Properties & properties, const double precompression) :
    _nodeId(node._id), _properties(properties), _precompression(precompression) {

    _l = _properties._lmeta;
    _L = pow(pow((_properties._delta/2),2)+pow(_l,2),0.5);

    _massMatrix(0,0) = _properties._density*_L;

  }

  double
  computeEnergy(const NodalDisplacements & displacements, const double time) const {
    ignoreUnusedVariables(time);

    double energy = 0;

    double u = displacements[0](0)- _precompression;
    double Lp = pow((pow((_l+u),2) + pow((_properties._delta/2),2)),0.5);

    energy = _properties._stiffness*pow((Lp-_L),2);

    return energy;

  }

  Forces
  computeForces(const NodalDisplacements & displacements, const double time) const {
    ignoreUnusedVariables(time);

    Forces forces;
    forces[0].fill(0);

    double u = displacements[0](0)-_precompression;
    double Lp = pow((pow((_l+u),2) + pow((_properties._delta/2),2)),0.5);
    double X = Lp - _L;

    forces[0](0) = 2*_properties._stiffness*X*(_l+u)/Lp;

    return forces;
  }

  StiffnessMatrix
  computeStiffnessMatrix(const NodalDisplacements & displacements, const double time)
    const {
    ignoreUnusedVariables(time);

    StiffnessMatrix stiffnessMatrix;
    stiffnessMatrix.fill(0);

    double u = displacements[0](0)-_precompression;
    double Lp = pow((pow((_l+u),2) + pow((_properties._delta/2),2)),0.5);
    stiffnessMatrix(0,0) = 2*_properties._stiffness*(1-(_L*pow((_properties._delta/2),2))/pow(Lp,3));

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
  Properties _properties;
  double _precompression;
  double _L;
  double _l;
  MassMatrix _massMatrix;

};

}
}

#endif //ELEMENT_NON_LINEAR_SPRING
