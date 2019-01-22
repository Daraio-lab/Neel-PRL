// -*- C++ -*-
#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "Definitions.h"
#include "Utilities.h"
#include "MeshUtilities.h"

namespace Assemblers {

namespace Utilities {

template <class Element>
void
assembleEnergyUtility(const vector<typename Element::Vector> & displacements,
                      const vector<Element> & elements,
                      const double time,
                      double * energy) {
  double localEnergy = 0.0;
  for (size_t elementIndex = 0; elementIndex < elements.size();
       ++elementIndex) {
    // get element displacements
    const Element & element = elements[elementIndex];
    const array<size_t, Element::NumberOfNodes> elementNodeIds =
      element.getNodeIds();
    const array<typename Element::Vector, Element::NumberOfNodes>
      elementDisplacements =
      ::Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds,
                                                                  displacements);

    // assemble local contribution
    localEnergy += element.computeEnergy(elementDisplacements, time);
  }
  *energy += localEnergy;
}

template <class Element>
void
updateInternalVariablesUtility(const vector<typename Element::Vector> & displacements,
                               const double time,
                               vector<Element> * elements) {
  for (size_t elementIndex = 0; elementIndex < elements->size(); ++elementIndex) {
    Element & element = elements->at(elementIndex);
    const array<size_t, Element::NumberOfNodes> elementNodeIds =
      element.getNodeIds();
    const array<typename Element::Vector,
                Element::NumberOfNodes> elementDisplacements =
      ::Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds,
                                                                  displacements);
    element.updateInternalVariables(elementDisplacements, time);

  }
}

template <class Element>
void
updateOrderParametersUtility(const vector<typename Element::OrderParameter> & orderParameters,
                             vector<Element> * elements) {
  for (size_t elementIndex = 0; elementIndex < elements->size(); ++elementIndex) {
    Element & element = elements->at(elementIndex);
    const array<size_t, Element::NumberOfNodes> elementNodeIds =
      element.getNodeIds();
    const array<typename Element::OrderParameter,
                Element::NumberOfNodes> elementOrderParameters =
      ::Utilities::getElementOrderParametersFromGlobalList<Element>(elementNodeIds,
                                                                    orderParameters);
    element.updateOrderParameters(elementOrderParameters);

  }
}

template <class Element>
void
assembleForceVectorUtility(const vector<typename Element::Vector> & displacements,
                           const vector<Element> & elements,
                           const double time,
                           Eigen::VectorXd * forces) {
  const size_t numberOfDofs = displacements.size() * Element::DegreesOfFreedom;
  if (forces->size() != numberOfDofs) {
    errorStatement("cannot assembleForceVectorUtility with forces->size() = %zu "
                   "and numberOfDofs = %zu\n", forces->size(),
                   numberOfDofs);
    exit(1);
  }

  for (size_t elementIndex = 0; elementIndex < elements.size();
       ++elementIndex) {
    // get element displacements
    const Element & element = elements[elementIndex];
    const array<size_t, Element::NumberOfNodes> elementNodeIds =
      element.getNodeIds();
    const array<typename Element::Vector, Element::NumberOfNodes>
      elementDisplacements =
      ::Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds,
                                                                  displacements);

    try {
      // assemble local contribution
      typename Element::Forces elementForces =
        element.computeForces(elementDisplacements, time);
      for (size_t nodeIndex = 0; nodeIndex < elementNodeIds.size();
           ++nodeIndex) {
        size_t nodeId = elementNodeIds[nodeIndex];
        for (size_t i = 0; i < Element::DegreesOfFreedom; ++i) {
          (*forces)(nodeId * Element::DegreesOfFreedom + i) +=
            elementForces[nodeIndex](i);
        }
      }
    } catch (std::exception & e) {
      errorStatement("exception caught trying to assemble element %zu's force vector\n",
                     elementIndex);
      throw;
    }
  }
}

template <class Element>
void
assembleOrderParameterForceVectorUtility(const vector<typename Element::Vector> & displacements,
                                         const vector<Element> & elements,
                                         const double time,
                                         Eigen::VectorXd * forces) {
  const size_t numberOfDofs = displacements.size() * Element::OrderParameterDofs;
  if (forces->size() != numberOfDofs) {
    errorStatement("cannot assembleForceVectorUtility with forces->size() = %zu "
                   "and numberOfDofs = %zu\n", forces->size(),
                   numberOfDofs);
    exit(1);
  }

  for (size_t elementIndex = 0; elementIndex < elements.size();
       ++elementIndex) {
    // get element displacements
    const Element & element = elements[elementIndex];
    const array<size_t, Element::NumberOfNodes> elementNodeIds =
      element.getNodeIds();
    const array<typename Element::Vector, Element::NumberOfNodes>
      elementDisplacements =
      ::Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds,
                                                                  displacements);

    try {
      // assemble local contribution
      typename Element::OrderParameterForces elementForces =
        element.computeOrderParameterForces(elementDisplacements, time);
      for (size_t nodeIndex = 0; nodeIndex < elementNodeIds.size();
           ++nodeIndex) {
        size_t nodeId = elementNodeIds[nodeIndex];
        for (size_t i = 0; i < Element::OrderParameterDofs; ++i) {
          (*forces)(nodeId * Element::OrderParameterDofs + i) +=
            elementForces[nodeIndex](i);
        }
      }
    } catch (std::exception & e) {
      errorStatement("exception caught trying to assemble element %zu's force vector\n",
                     elementIndex);
      throw;
    }
  }
}

template <class Element>
struct StiffnessMatrixExtractor {

  typedef typename Element::StiffnessMatrix MatrixType;

  const vector<typename Element::Vector> & _displacements;
  const double _time;
  StiffnessMatrixExtractor(const vector<typename Element::Vector> & displacements,
                           const double time) :
    _displacements(displacements), _time(time) {
  }

  MatrixType
  operator()(const Element & element) const {
    const array<size_t, Element::NumberOfNodes> elementNodeIds =
      element.getNodeIds();
    const array<typename Element::Vector, Element::NumberOfNodes>
      elementDisplacements =
      ::Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds,
                                                                  _displacements);
    return element.computeStiffnessMatrix(elementDisplacements, _time);
  }
private:
  StiffnessMatrixExtractor();
};

template <class Element>
struct LumpedMassMatrixExtractor {

  typedef typename Element::MassMatrix MatrixType;

  MatrixType
  operator()(const Element & element) const {
    return element.computeLumpedMassMatrix();
  }
};

template <class Element>
struct ConsistentMassMatrixExtractor {

  typedef typename Element::MassMatrix MatrixType;

  MatrixType
  operator()(const Element & element) const {
    return element.computeConsistentMassMatrix();
  }
};

template <class Element>
struct LumpedMobilityMatrixExtractor {

  typedef typename Element::MobilityMatrix MatrixType;

  MatrixType
  operator()(const Element & element) const {
    return element.computeLumpedMobilityMatrix();
  }
};

template <class Element>
struct ConsistentMobilityMatrixExtractor {

  typedef typename Element::MobilityMatrix MatrixType;

  MatrixType
  operator()(const Element & element) const {
    return element.computeConsistentMobilityMatrix();
  }
};

template <class Element, class MatrixExtractor>
void
assembleMatrixUtility(const size_t numberOfNodes,
                      const vector<Element> & elements,
                      const MatrixExtractor & matrixExtractor,
                      Eigen::SparseMatrix<double> * matrix) {
  const size_t numberOfDofs = numberOfNodes * Element::DegreesOfFreedom;
  // doesn't matter now for sparse matricies does it?
  if (matrix->rows() != numberOfDofs || matrix->cols() != numberOfDofs) {
    throwException("cannot assembleMatrix with matrix->rows() = %d "
                   "and matrix->rows() = %d and numberOfDofs = %zu,"
                   "numberOfNodes = %lu,"
                   "DegreesOfFreedom = %u\n",
                   matrix->rows(), matrix->cols(), numberOfDofs,
                   numberOfNodes,
                   Element::DegreesOfFreedom);
  }
  std::vector<Eigen::Triplet<double> > tripletList;
  Eigen::SparseMatrix<double> tempStorageMatrix(matrix->rows(),matrix->cols());
  // esimate the size based on the number of non zero entries
  // in the first matrix
  if (elements.size() > 0){
    const Element & exampleElement = elements[0];
    const typename MatrixExtractor::MatrixType elementMatrixExample =
      matrixExtractor(exampleElement);
    // SO this could be better thought out but it's the best I have right now
    tripletList.reserve(elements.size() * 20);
  }
  for (size_t elementIndex = 0; elementIndex < elements.size();
       ++elementIndex) {
    const Element & element = elements[elementIndex];
    const array<size_t, Element::NumberOfNodes> elementNodeIds =
      element.getNodeIds();

    try {
      const typename MatrixExtractor::MatrixType elementMatrix =
        matrixExtractor(element);

      // assemble local contribution
      for (size_t nodeIndex1 = 0; nodeIndex1 < elementNodeIds.size();
           ++nodeIndex1) {
        const size_t nodeId1 = elementNodeIds[nodeIndex1];
        for (size_t nodeIndex2 = 0; nodeIndex2 < elementNodeIds.size();
             ++nodeIndex2) {
          const size_t nodeId2 = elementNodeIds[nodeIndex2];
          for (size_t i = 0; i < Element::DegreesOfFreedom; ++i) {
            for (size_t j = 0; j < Element::DegreesOfFreedom; ++j) {
              // heaven help me if this is wrong ... seriously
              tripletList.push_back(Eigen::Triplet<double>
                                    (nodeId1 * Element::DegreesOfFreedom + i,
                                     nodeId2 * Element::DegreesOfFreedom + j,
                                     elementMatrix(nodeIndex1 *
                                                   Element::DegreesOfFreedom + i,
                                                   nodeIndex2 *
                                                   Element::DegreesOfFreedom + j)));
            }
          }
        }
      }
    } catch (std::exception & e) {
      errorStatement("exception caught trying to assemble element %zu's matrix\n",
                     elementIndex);
      throw;
    }
  }
  tempStorageMatrix.setFromTriplets(tripletList.begin(),tripletList.end());
  (*matrix) += tempStorageMatrix;
}

template <class Element, class MatrixExtractor>
void
assembleMechanicalMatrixUtility(const size_t numberOfNodes,
                                const vector<Element> & elements,
                                const MatrixExtractor & matrixExtractor,
                                Eigen::SparseMatrix<double> * matrix) {
  const size_t numberOfDofs = numberOfNodes * Element::MechanicalDegreesOfFreedom;
  // doesn't matter now for sparse matricies does it?
  if (matrix->rows() != numberOfDofs || matrix->cols() != numberOfDofs) {
    errorStatement("cannot assembleMatrix with matrix->rows() = %d "
                   "and matrix->rows() = %d and numberOfDofs = %zu\n",
                   matrix->rows(), matrix->cols(), numberOfDofs);
    exit(1);
  }
  std::vector<Eigen::Triplet<double> > tripletList;
  Eigen::SparseMatrix<double> tempStorageMatrix(matrix->rows(),matrix->cols());
  // esimate the size based on the number of non zero entries
  // in the first matrix
  if (elements.size() > 0){
    const Element & exampleElement = elements[0];
    const typename MatrixExtractor::MatrixType elementMatrixExample =
      matrixExtractor(exampleElement);
    // SO this could be better thought out but it's the best I have right now
    tripletList.reserve(elements.size() * 20);
  }
  for (size_t elementIndex = 0; elementIndex < elements.size();
       ++elementIndex) {
    const Element & element = elements[elementIndex];
    const array<size_t, Element::NumberOfNodes> elementNodeIds =
      element.getNodeIds();

    try {
      const typename MatrixExtractor::MatrixType elementMatrix =
        matrixExtractor(element);

      // assemble local contribution
      for (size_t nodeIndex1 = 0; nodeIndex1 < elementNodeIds.size();
           ++nodeIndex1) {
        const size_t nodeId1 = elementNodeIds[nodeIndex1];
        for (size_t nodeIndex2 = 0; nodeIndex2 < elementNodeIds.size();
             ++nodeIndex2) {
          const size_t nodeId2 = elementNodeIds[nodeIndex2];
          for (size_t i = 0; i < Element::MechanicalDegreesOfFreedom; ++i) {
            for (size_t j = 0; j < Element::MechanicalDegreesOfFreedom; ++j) {
              // heaven help me if this is wrong ... seriously
              tripletList.push_back(Eigen::Triplet<double>
                                    (nodeId1 * Element::MechanicalDegreesOfFreedom + i,
                                     nodeId2 * Element::MechanicalDegreesOfFreedom + j,
                                     elementMatrix(nodeIndex1 *
                                                   Element::MechanicalDegreesOfFreedom + i,
                                                   nodeIndex2 *
                                                   Element::MechanicalDegreesOfFreedom + j)));
            }
          }
        }
      }
    } catch (std::exception & e) {
      errorStatement("exception caught trying to assemble element %zu's matrix\n",
                     elementIndex);
      throw;
    }
  }
  tempStorageMatrix.setFromTriplets(tripletList.begin(),tripletList.end());
  (*matrix) += tempStorageMatrix;
}

template <class Element, class MatrixExtractor>
void
assembleOrderParameterMatrixUtility(const size_t numberOfNodes,
                                    const vector<Element> & elements,
                                    const MatrixExtractor & matrixExtractor,
                                    Eigen::SparseMatrix<double> * matrix) {
  const size_t numberOfDofs = numberOfNodes * Element::OrderParameterDofs;
  // doesn't matter now for sparse matricies does it?
  if (matrix->rows() != numberOfDofs || matrix->cols() != numberOfDofs) {
    errorStatement("cannot assembleMatrix with matrix->rows() = %d "
                   "and matrix->rows() = %d and numberOfDofs = %zu\n",
                   matrix->rows(), matrix->cols(), numberOfDofs);
    exit(1);
  }
  std::vector<Eigen::Triplet<double> > tripletList;
  Eigen::SparseMatrix<double> tempStorageMatrix(matrix->rows(),matrix->cols());
  // esimate the size based on the number of non zero entries
  // in the first matrix
  if (elements.size() > 0){
    const Element & exampleElement = elements[0];
    const typename MatrixExtractor::MatrixType elementMatrixExample =
      matrixExtractor(exampleElement);
    // SO this could be better thought out but it's the best I have right now
    tripletList.reserve(elements.size() * 20);
  }
  for (size_t elementIndex = 0; elementIndex < elements.size();
       ++elementIndex) {
    const Element & element = elements[elementIndex];
    const array<size_t, Element::NumberOfNodes> elementNodeIds =
      element.getNodeIds();

    try {
      const typename MatrixExtractor::MatrixType elementMatrix =
        matrixExtractor(element);

      // assemble local contribution
      for (size_t nodeIndex1 = 0; nodeIndex1 < elementNodeIds.size();
           ++nodeIndex1) {
        const size_t nodeId1 = elementNodeIds[nodeIndex1];
        for (size_t nodeIndex2 = 0; nodeIndex2 < elementNodeIds.size();
             ++nodeIndex2) {
          const size_t nodeId2 = elementNodeIds[nodeIndex2];
          for (size_t i = 0; i < Element::OrderParameterDofs; ++i) {
            for (size_t j = 0; j < Element::OrderParameterDofs; ++j) {
              // heaven help me if this is wrong ... seriously
              tripletList.push_back(Eigen::Triplet<double>
                                    (nodeId1 * Element::OrderParameterDofs + i,
                                     nodeId2 * Element::OrderParameterDofs + j,
                                     elementMatrix(nodeIndex1 *
                                                   Element::OrderParameterDofs + i,
                                                   nodeIndex2 *
                                                   Element::OrderParameterDofs + j)));
            }
          }
        }
      }
    } catch (std::exception & e) {
      errorStatement("exception caught trying to assemble element %zu's matrix\n",
                     elementIndex);
      throw;
    }
  }
  tempStorageMatrix.setFromTriplets(tripletList.begin(),tripletList.end());
  (*matrix) += tempStorageMatrix;
}

template <class Element>
void
assembleStiffnessMatrixUtility(const vector<typename Element::Vector> & displacements,
                               const vector<Element> & elements,
                               const double time,
                               Eigen::SparseMatrix<double> * matrix) {
  assembleMatrixUtility<Element,
                        StiffnessMatrixExtractor<Element> >(displacements.size(),
                                                            elements,
                                                            StiffnessMatrixExtractor<Element>(displacements, time),
                                                            matrix);
}

template <class Element>
void
assembleLumpedMassMatrixUtility(const size_t numberOfNodes,
                                const vector<Element> & elements,
                                Eigen::SparseMatrix<double> * matrix) {
  assembleMatrixUtility<Element,
                        LumpedMassMatrixExtractor<Element> >(numberOfNodes,
                                                             elements,
                                                             LumpedMassMatrixExtractor<Element>(),
                                                             matrix);
}

template <class Element>
void
assembleConsistentMassMatrixUtility(const size_t numberOfNodes,
                                    const vector<Element> & elements,
                                    Eigen::MatrixXd * matrix) {
  assembleMatrixUtility<Element,
                        ConsistentMassMatrixExtractor<Element> >(numberOfNodes,
                                                                 elements,
                                                                 ConsistentMassMatrixExtractor<Element>(),
                                                                 matrix);
}

template <class Element>
void
assembleMechanicalLumpedMassMatrixUtility(const size_t numberOfNodes,
                                          const vector<Element> & elements,
                                          Eigen::SparseMatrix<double> * matrix) {
  assembleMechanicalMatrixUtility<Element,
                                  LumpedMassMatrixExtractor<Element> >(numberOfNodes,
                                                                       elements,
                                                                       LumpedMassMatrixExtractor<Element>(),
                                                                       matrix);
}

template <class Element>
void
assembleMechanicalConsistentMassMatrixUtility(const size_t numberOfNodes,
                                              const vector<Element> & elements,
                                              Eigen::MatrixXd * matrix) {
  assembleMechanicalMatrixUtility<Element,
                                  ConsistentMassMatrixExtractor<Element> >(numberOfNodes,
                                                                           elements,
                                                                           ConsistentMassMatrixExtractor<Element>(),
                                                                           matrix);
}

template <class Element>
void
assembleLumpedMobilityMatrixUtility(const size_t numberOfNodes,
                                    const vector<Element> & elements,
                                    Eigen::SparseMatrix<double> * matrix) {
  assembleOrderParameterMatrixUtility<Element,
                                      LumpedMobilityMatrixExtractor<Element> >(numberOfNodes,
                                                                               elements,
                                                                               LumpedMobilityMatrixExtractor<Element>(),
                                                                               matrix);
}

template <class Element>
void
assembleConsistentMobilityMatrixUtility(const size_t numberOfNodes,
                                        const vector<Element> & elements,
                                        Eigen::MatrixXd * matrix) {
  assembleOrderParameterMatrixUtility<Element,
                                      ConsistentMobilityMatrixExtractor<Element> >(numberOfNodes,
                                                                                   elements,
                                                                                   ConsistentMobilityMatrixExtractor<Element>(),
                                                                                   matrix);
}

template <class Element>
vector<typename Element::Stress>
computeElementStresses(const vector<typename Element::Vector> & displacements,
                       const vector<Element> & elements,
                       const double time) {

  vector<typename Element::Stress> allElementStresses(elements.size());
  for (size_t elementIndex = 0; elementIndex < elements.size(); ++elementIndex) {
    const Element & element = elements[elementIndex];
    array<size_t, Element::NumberOfNodes> elementNodeIds = element.getNodeIds();
    array<typename Element::Vector, Element::NumberOfNodes> elementDisplacements =
      ::Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds,
                                                                  displacements);
    array<typename Element::Stress, Element::QuadPoints> elementStresses =
      element.computeStressesAtGaussPoints(elementDisplacements, time);
    typename Element::Stress average = Element::Stress::Zero();
    for (size_t qpIndex = 0; qpIndex < Element::QuadPoints; ++qpIndex) {
      average += elementStresses[qpIndex];
    }
    average /= Element::QuadPoints;
    allElementStresses[elementIndex] = average;
  }
  return allElementStresses;
}

template <class Element>
vector<typename Element::Stress>
computeNodalStresses(const vector<typename Element::Vector> & displacements,
                     const vector<Element> & elements,
                     const double time) {

  const size_t numberOfNodes = displacements.size();

  vector<typename Element::Stress> nodalStresses(numberOfNodes,
                                                 Element::Stress::Zero());
  vector<double> volumeSums(numberOfNodes, 0);

  vector<typename Element::Stress> elementStresses =
    computeElementStresses(displacements,elements, time);

  for (size_t elementIndex = 0; elementIndex < elements.size(); ++elementIndex) {
    const Element & element = elements[elementIndex];
    array<size_t, Element::NumberOfNodes> elementNodeIds = element.getNodeIds();
    array<double, Element::NumberOfNodes> elementWeights =
      element.computeNodalWeights();
    const typename Element::Stress & elementStress = elementStresses[elementIndex];

    for (size_t nodeIndex = 0; nodeIndex < elementNodeIds.size(); nodeIndex++) {
      nodalStresses[elementNodeIds[nodeIndex]] +=
        elementStress / elementWeights[nodeIndex];
      volumeSums[elementNodeIds[nodeIndex]] += 1./elementWeights[nodeIndex];
    }
  }

  for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++) {
    nodalStresses[nodeIndex] /= volumeSums[nodeIndex];
  }
  return nodalStresses;
}

template <class Element>
vector<typename Element::Strain>
computeElementStrains(const vector<typename Element::Vector> & displacements,
                      const vector<Element> & elements) {

  vector<typename Element::Strain> allElementStrains(elements.size());
  for (size_t elementIndex = 0; elementIndex < elements.size(); ++elementIndex) {
    const Element & element = elements[elementIndex];
    array<size_t, Element::NumberOfNodes> elementNodeIds = element.getNodeIds();
    array<typename Element::Vector, Element::NumberOfNodes> elementDisplacements =
      ::Utilities::getElementDisplacementsFromGlobalList<Element>(elementNodeIds,
                                                                  displacements);
    array<typename Element::Strain, Element::QuadPoints> elementStrains =
      element.computeStrainsAtGaussPoints(elementDisplacements);
    typename Element::Strain average = Element::Strain::Zero();
    for (size_t qpIndex = 0; qpIndex < Element::QuadPoints; ++qpIndex) {
      average += elementStrains[qpIndex];
    }
    average /= Element::QuadPoints;
    allElementStrains[elementIndex] = average;
  }
  return allElementStrains;
}

template <class Element>
vector<typename Element::Strain>
computeNodalStrains(const vector<typename Element::Vector> & displacements,
                    const vector<Element> & elements) {

  const size_t numberOfNodes = displacements.size();

  vector<typename Element::Strain> nodalStrains(numberOfNodes,
                                                Element::Strain::Zero());
  vector<double> volumeSums(numberOfNodes, 0);

  vector<typename Element::Strain> elementStrains =
    computeElementStrains(displacements, elements);

  for (size_t elementIndex = 0; elementIndex < elements.size(); ++elementIndex) {
    const Element & element = elements[elementIndex];
    array<size_t, Element::NumberOfNodes> elementNodeIds = element.getNodeIds();
    array<double, Element::NumberOfNodes> elementWeights =
      element.computeNodalWeights();
    const typename Element::Strain & elementStrain = elementStrains[elementIndex];

    for (size_t nodeIndex = 0; nodeIndex < elementNodeIds.size(); nodeIndex++) {
      nodalStrains[elementNodeIds[nodeIndex]] +=
        elementStrain / elementWeights[nodeIndex];
      volumeSums[elementNodeIds[nodeIndex]] += 1./elementWeights[nodeIndex];
    }
  }

  for (size_t nodeIndex = 0; nodeIndex < numberOfNodes; nodeIndex++){
    nodalStrains[nodeIndex] /= volumeSums[nodeIndex];
  }
  return nodalStrains;
}

}


template <class ... Elements>
class ExternalForceAssembler {
};

template <class E0>
class ExternalForceAssembler <E0> {

public:

  typedef typename E0::Vector            ElementVector;

  vector<E0> _elements0;

  void
  addElement(const E0 & e) {
    _elements0.push_back(e);
  }

  ExternalForceAssembler () {
  }

  void
  assembleEnergy(const vector<ElementVector> & displacements,
                 const double time,
                 double * energy) const {
    Assemblers::Utilities::assembleEnergyUtility(displacements, _elements0, time, energy);
  }

  void
  assembleForceVector(const vector<ElementVector> & displacements,
                      const double time,
                      Eigen::VectorXd * forceVector) const {
    Assemblers::Utilities::assembleForceVectorUtility(displacements, _elements0, time,
                                                      forceVector);
  }

  void
  assembleStiffnessMatrix(const vector<ElementVector> & displacements,
                          const double time,
                          Eigen::SparseMatrix<double> * stiffnessMatrix) const {
    Assemblers::Utilities::assembleStiffnessMatrixUtility(displacements, _elements0, time,
                                                          stiffnessMatrix);
  }

};

template <class E0, class E1>
class ExternalForceAssembler <E0, E1> {

public:

  typedef typename E0::Vector    ElementVector;

  vector<E0> _elements0;
  vector<E1> _elements1;

  void
  addElement(const E0 & e) {
    _elements0.push_back(e);
  }

  void
  addElement(const E1 & e) {
    _elements1.push_back(e);
  }

  ExternalForceAssembler () {
  }

  void
  assembleEnergy(const vector<ElementVector> & displacements,
                 const double time,
                 double * energy) const {
    Assemblers::Utilities::assembleEnergyUtility(displacements, _elements0, time, energy);
    Assemblers::Utilities::assembleEnergyUtility(displacements, _elements1, time, energy);
  }

  void
  assembleForceVector(const vector<ElementVector> & displacements,
                      const double time,
                      Eigen::VectorXd * forceVector) const {
    Assemblers::Utilities::assembleForceVectorUtility(displacements, _elements0, time,
                                                      forceVector);
    Assemblers::Utilities::assembleForceVectorUtility(displacements, _elements1, time,
                                                      forceVector);
  }

  void
  assembleStiffnessMatrix(const vector<ElementVector> & displacements,
                          const double time,
                          Eigen::SparseMatrix<double> * stiffnessMatrix) const {
    Assemblers::Utilities::assembleStiffnessMatrixUtility(displacements, _elements0, time,
                                                          stiffnessMatrix);
    Assemblers::Utilities::assembleStiffnessMatrixUtility(displacements, _elements1, time,
                                                          stiffnessMatrix);
  }

};

template <class E0, class E1, class E2>
class ExternalForceAssembler <E0, E1, E2> {

public:

  typedef typename E0::Vector    ElementVector;

  vector<E0> _elements0;
  vector<E1> _elements1;
  vector<E2> _elements2;

  void
  addElement(const E0 & e) {
    _elements0.push_back(e);
  }

  void
  addElement(const E1 & e) {
    _elements1.push_back(e);
  }

  void
  addElement(const E2 & e) {
    _elements2.push_back(e);
  }

  ExternalForceAssembler () {
  }

  void
  assembleEnergy(const vector<ElementVector> & displacements,
                 const double time,
                 double * energy) const {
    Assemblers::Utilities::assembleEnergyUtility(displacements, _elements0, time, energy);
    Assemblers::Utilities::assembleEnergyUtility(displacements, _elements1, time, energy);
    Assemblers::Utilities::assembleEnergyUtility(displacements, _elements2, time, energy);
  }

  void
  assembleForceVector(const vector<ElementVector> & displacements,
                      const double time,
                      Eigen::VectorXd * forceVector) const {
    Assemblers::Utilities::assembleForceVectorUtility(displacements, _elements0, time,
                                                      forceVector);
    Assemblers::Utilities::assembleForceVectorUtility(displacements, _elements1, time,
                                                      forceVector);
    Assemblers::Utilities::assembleForceVectorUtility(displacements, _elements2, time,
                                                      forceVector);
  }

  void
  assembleStiffnessMatrix(const vector<ElementVector> & displacements,
                          const double time,
                          Eigen::SparseMatrix<double> * stiffnessMatrix) const {
    Assemblers::Utilities::assembleStiffnessMatrixUtility(displacements, _elements0, time,
                                                          stiffnessMatrix);
    Assemblers::Utilities::assembleStiffnessMatrixUtility(displacements, _elements1, time,
                                                          stiffnessMatrix);
    Assemblers::Utilities::assembleStiffnessMatrixUtility(displacements, _elements2, time,
                                                          stiffnessMatrix);
  }

};

template <class Element>
class ExternalForceAssemblerPlacebo {

public:

  void
  addElement(const Element & e) {
  }

  ExternalForceAssemblerPlacebo () {
  }

  void
  assembleEnergy(const vector<typename Element::Vector> & displacements,
                 const double time,
                 double * energy) const {
  }

  void
  assembleForceVector(const vector<typename Element::Vector> & displacements,
                      const double time,
                      Eigen::VectorXd * forceVector) const {
    ignoreUnusedVariables(displacements, time, forceVector);
  }

  void
  assembleStiffnessMatrix(const vector<typename Element::Vector> & displacements,
                          const double time,
                          Eigen::SparseMatrix<double> * stiffnessMatrix) const {
    ignoreUnusedVariables(displacements, time, stiffnessMatrix);
  }

};



template <class PhysicalElementType,
          class ExternalForceAssembler = ExternalForceAssemblerPlacebo<PhysicalElementType> >
class Assembler {

public:
  typedef PhysicalElementType                           PhysicalElement;
  typedef typename PhysicalElementType::Vector          ElementVector;
  typedef typename PhysicalElementType::Stress          ElementStress;
  typedef typename PhysicalElementType::Strain          ElementStrain;
  typedef typename PhysicalElementType::TangentMatrix   ElementTangentMatrix;
  typedef typename PhysicalElementType::StiffnessMatrix ElementStiffnessMatrix;
  typedef typename PhysicalElementType::Point           ElementPoint;
  typedef typename PhysicalElementType::Node            ElementNode;
  static const unsigned int      ElementNumberOfNodes = PhysicalElementType::NumberOfNodes;
  static const unsigned int      SpatialDimension     = PhysicalElementType::SpatialDimension;
  static const unsigned int      DegreesOfFreedom     = PhysicalElementType::DegreesOfFreedom;

  Assembler(const size_t numberOfNodes) : _numberOfNodes(numberOfNodes) {
  }

  template <class QuadratureRule, class MaterialModel>
  Assembler(const SingleElementMesh<PhysicalElementType> & mesh,
            const typename PhysicalElementType::Properties & elementProperties,
            const QuadratureRule & quadratureRule,
            const MaterialModel & materialModel) :
    _numberOfNodes(mesh._nodes.size()) {
    addAllElementsFromMesh(mesh, elementProperties,
                           quadratureRule, materialModel);
  }

  template <class QuadratureRule, class MaterialModel>
  void
  addAllElementsFromMesh(const SingleElementMesh<PhysicalElementType> & mesh,
                         const typename PhysicalElementType::Properties & elementProperties,
                         const QuadratureRule & quadratureRule,
                         const MaterialModel & materialModel) {
    const unsigned int numberOfElements = mesh._connectivity.size();
    for (unsigned int elementIndex = 0; elementIndex < numberOfElements;
         ++elementIndex) {
      const array<typename PhysicalElementType::Node, PhysicalElementType::NumberOfNodes> elementNodes =
        MeshUtilities::getNodesFromMesh<PhysicalElementType>(mesh, elementIndex);
      addElement(PhysicalElementType(elementNodes, elementProperties,
                                     &quadratureRule, &materialModel));
    }
  }

  // external force assembler has to handle the general case
  template <class E>
  void
  addElement(const E & e) {
    _externalForceAssembler.addElement(e);
  }

  // specialization for the physical element type
  void
  addElement(const PhysicalElementType & e) {
    _physicalElements.push_back(e);
  }

  void
  assembleEnergy(const vector<ElementVector> & displacements,
                 const double time,
                 double * energy) const {
    *energy = 0;
    try {
      Assemblers::Utilities::assembleEnergyUtility(displacements, _physicalElements,
                                                   time, energy);
    } catch (std::exception & e) {
      errorStatement("exception caught trying to assemble the energy for "
                     "the physical elements\n");
      throw;
    }
    try {
      _externalForceAssembler.assembleEnergy(displacements, time, energy);
    } catch (std::exception & e) {
      errorStatement("exception caught trying to assemble the energy for "
                     "the external force assembler\n");
      throw;
    }
  }

  void
  assembleForceVector(const vector<ElementVector> & displacements,
                      const double time,
                      Eigen::VectorXd * forceVector) const {
    forceVector->fill(0);
    try {
      Assemblers::Utilities::assembleForceVectorUtility(displacements,
                                                        _physicalElements, time,
                                                        forceVector);
    } catch (std::exception & e) {
      errorStatement("exception caught trying to assemble the force vector for "
                     "the physical elements\n");
      throw;
    }
    try {
      _externalForceAssembler.assembleForceVector(displacements, time,
                                                  forceVector);
    } catch (std::exception & e) {
      errorStatement("exception caught trying to assemble the force vector for "
                     "the external force assembler\n");
      throw;
    }
  }

  void
  assembleStiffnessMatrix(const vector<ElementVector> & displacements,
                          const double time,
                          Eigen::SparseMatrix<double> * stiffnessMatrix) const {
    stiffnessMatrix->setZero();
    try {
      Assemblers::Utilities::assembleStiffnessMatrixUtility(displacements,
                                                            _physicalElements,
                                                            time,
                                                            stiffnessMatrix);
    } catch (std::exception & e) {
      errorStatement("exception caught trying to assemble the stiffness matrix for "
                     "the physical elements\n");
      throw;
    }
    try {
      _externalForceAssembler.assembleStiffnessMatrix(displacements, time,
                                                      stiffnessMatrix);
    } catch (std::exception & e) {
      errorStatement("exception caught trying to assemble the stiffness matrix for "
                     "the external force assembler\n");
      throw;
    }
  }

  void
  updateInternalVariables(const vector<ElementVector> & displacements,
                          const double time)  {
    Assemblers::Utilities::updateInternalVariablesUtility(displacements, time,
                                                          &_physicalElements);
  }

  Eigen::SparseMatrix<double>
  assembleLumpedMassMatrix() const {
    size_t numberOfDofs = _numberOfNodes * DegreesOfFreedom;
    Eigen::SparseMatrix<double> massMatrix(numberOfDofs,numberOfDofs);
    Assemblers::Utilities::assembleLumpedMassMatrixUtility(_numberOfNodes,
                                                           _physicalElements,
                                                           &massMatrix);
    return massMatrix;
  }

  Eigen::SparseMatrix<double>
  assembleConsistentMassMatrix() const {
    size_t numberOfDofs = _numberOfNodes * DegreesOfFreedom;
    Eigen::SparseMatrix<double> massMatrix(numberOfDofs, numberOfDofs);
    massMatrix.setZero();
    Assemblers::Utilities::assembleConsistentMassMatrixUtility(_numberOfNodes,
                                                               _physicalElements,
                                                               &massMatrix);
    return massMatrix;
  }

  vector<ElementStress>
  computeElementStresses(const vector<ElementVector> & displacements,
                         const double time) const {
    vector<ElementStress> elementStresses =
      Assemblers::Utilities::computeElementStresses(displacements,
                                                    _physicalElements, time);
    return elementStresses;
  }

  vector<ElementStress>
  computeNodalStresses(const vector<ElementVector> & displacements,
                       const double time) const {
    vector<ElementStress> nodalStresses =
      Assemblers::Utilities::computeNodalStresses(displacements,
                                                  _physicalElements, time);
    return nodalStresses;
  }

  vector<ElementStrain>
  computeElementStrains(const vector<ElementVector> & displacements) const {
    vector<ElementStrain> elementStrains =
      Assemblers::Utilities::computeElementStrains(displacements,
                                                   _physicalElements);
    return elementStrains;
  }

  vector<ElementStrain>
  computeNodalStrains(const vector<ElementVector> & displacements) const {
    vector<ElementStrain> nodalStrains =
      Assemblers::Utilities::computeNodalStrains(displacements,
                                                 _physicalElements);
    return nodalStrains;
  }

  const vector<PhysicalElementType> &
  getPhysicalElements() const {
    return _physicalElements;
  }

  size_t
  getNumberOfNodes() const {
    return _numberOfNodes;
  }

  Eigen::VectorXd
  allocateZeroedGlobalForceVector() const {
    const size_t numberOfDOFs = _numberOfNodes * DegreesOfFreedom;
    Eigen::VectorXd temp(numberOfDOFs);
    temp.fill(0);
    return temp;
  }

  Eigen::SparseMatrix<double>
  allocateZeroedStiffnessMatrix() const {
    const size_t numberOfDOFs = _numberOfNodes * DegreesOfFreedom;
    Eigen::SparseMatrix<double> temp(numberOfDOFs, numberOfDOFs);
    return temp;
  }

protected:
  vector<PhysicalElementType> _physicalElements;
  size_t _numberOfNodes;
  ExternalForceAssembler _externalForceAssembler;
};

template <class PhysicalElementType1,
          class PhysicalElementType2,
          class ExternalForceAssembler = ExternalForceAssemblerPlacebo<PhysicalElementType1> >
class AssemblerTwoPhysicalElements {

public:
  typedef PhysicalElementType1                           PhysicalElement;
  typedef typename PhysicalElementType1::Vector          ElementVector;
  typedef typename PhysicalElementType1::Stress          ElementStress;
  typedef typename PhysicalElementType1::Strain          ElementStrain;
  typedef typename PhysicalElementType1::TangentMatrix   PrimaryElementTangentMatrix;
  typedef typename PhysicalElementType1::TangentMatrix   SecondaryElementTangentMatrix;
  typedef typename PhysicalElementType1::StiffnessMatrix PrimaryElementStiffnessMatrix;
  typedef typename PhysicalElementType1::StiffnessMatrix SecondaryElementStiffnessMatrix;
  typedef typename PhysicalElementType1::Point           ElementPoint;
  typedef typename PhysicalElementType1::Node            ElementNode;
  static const unsigned int      PrimaryElementNumberOfNodes = PhysicalElementType1::NumberOfNodes;
  static const unsigned int      SecondaryElementNumberOfNodes = PhysicalElementType2::NumberOfNodes;
  static const unsigned int      SpatialDimension     = PhysicalElementType1::SpatialDimension;
  static const unsigned int      DegreesOfFreedom     = PhysicalElementType1::DegreesOfFreedom;

  AssemblerTwoPhysicalElements(const size_t numberOfNodes) : _numberOfNodes(numberOfNodes) {
  }

  // external force assembler has to handle the general case
  template <class E>
  void
  addElement(const E & e) {
    _externalForceAssembler.addElement(e);
  }

  // specialization for the physical element type
  void
  addElement(const PhysicalElementType1 & e) {
    _primaryPhysicalElements.push_back(e);
  }

  void
  addElement(const PhysicalElementType2 & e) {
    _secondaryPhysicalElements.push_back(e);
  }

  void
  assembleEnergy(const vector<ElementVector> & displacements,
                 const double time,
                 double * energy) const {
    Assemblers::Utilities::assembleEnergyUtility(displacements, _primaryPhysicalElements,
                                                 time, energy);
    Assemblers::Utilities::assembleEnergyUtility(displacements, _secondaryPhysicalElements,
                                                 time, energy);
    _externalForceAssembler.assembleEnergy(displacements, time, energy);
  }

  void
  assembleForceVector(const vector<ElementVector> & displacements,
                      const double time,
                      Eigen::VectorXd * forceVector) const {
    Assemblers::Utilities::assembleForceVectorUtility(displacements,
                                                      _primaryPhysicalElements, time,
                                                      forceVector);
    Assemblers::Utilities::assembleForceVectorUtility(displacements,
                                                      _secondaryPhysicalElements, time,
                                                      forceVector);
    _externalForceAssembler.assembleForceVector(displacements, time,
                                                forceVector);
  }

  void
  assembleStiffnessMatrix(const vector<ElementVector> & displacements,
                          const double time,
                          Eigen::SparseMatrix<double> * stiffnessMatrix) const {
    Assemblers::Utilities::assembleStiffnessMatrixUtility(displacements,
                                                          _primaryPhysicalElements, time,
                                                          stiffnessMatrix);
    Assemblers::Utilities::assembleStiffnessMatrixUtility(displacements,
                                                          _secondaryPhysicalElements, time,
                                                          stiffnessMatrix);
    _externalForceAssembler.assembleStiffnessMatrix(displacements, time,
                                                    stiffnessMatrix);
  }

  void
  updateInternalVariables(const vector<ElementVector> & displacements,
                          const double time)  {
    Assemblers::Utilities::updateInternalVariablesUtility(displacements,time,
                                                          &_primaryPhysicalElements);
    Assemblers::Utilities::updateInternalVariablesUtility(displacements,time,
                                                          &_secondaryPhysicalElements);
  }

  Eigen::SparseMatrix<double>
  assembleLumpedMassMatrix() const {
    size_t numberOfDofs = _numberOfNodes * DegreesOfFreedom;
    Eigen::SparseMatrix<double> massMatrix(numberOfDofs,numberOfDofs);
    Assemblers::Utilities::assembleLumpedMassMatrixUtility(_numberOfNodes,
                                                           _primaryPhysicalElements,
                                                           &massMatrix);
    Assemblers::Utilities::assembleLumpedMassMatrixUtility(_numberOfNodes,
                                                           _secondaryPhysicalElements,
                                                           &massMatrix);
    return massMatrix;
  }

  Eigen::SparseMatrix<double>
  assembleConsistentMassMatrix() const {
    size_t numberOfDofs = _numberOfNodes * DegreesOfFreedom;
    Eigen::SparseMatrix<double> massMatrix(numberOfDofs, numberOfDofs);
    Assemblers::Utilities::assembleConsistentMassMatrixUtility(_numberOfNodes,
                                                               _primaryPhysicalElements,
                                                               &massMatrix);
    Assemblers::Utilities::assembleConsistentMassMatrixUtility(_numberOfNodes,
                                                               _secondaryPhysicalElements,
                                                               &massMatrix);
    return massMatrix;
  }

  vector<ElementStress>
  computeElementStresses(const vector<ElementVector> & displacements,
                         const double time) {
    // TODO THIS NEEDS TO BE UPDATED FOR TWO ELEMENTS
    errorStatement("computeElementStresses needs to be updated for the two element "
                   "assembler class\n");
    vector<ElementStress> elementStresses =
      Assemblers::Utilities::computeElementStresses(displacements,_primaryPhysicalElements, time);
    return elementStresses;
  }

  vector<ElementStress>
  computeNodalStresses(const vector<ElementVector> & displacements,
                       const double time) {
    // TODO THIS NEEDS TO BE UPDATED FOR TWO ELEMENTS
    errorStatement("computeNodalStresses needs to be updated for the two element "
                   "assembler class\n");
    vector<ElementStress> nodalStresses =
      Assemblers::Utilities::computeNodalStresses(displacements,_primaryPhysicalElements, time);
    return nodalStresses;
  }

  vector<ElementStrain>
  computeElementStrains(const vector<ElementVector> & displacements) {
    // TODO THIS NEEDS TO BE UPDATED FOR TWO ELEMENTS
    errorStatement("computeElementStrains needs to be updated for the two element "
                   "assembler class\n");
    vector<ElementStrain> elementStrains =
      Assemblers::Utilities::computeElementStrains(displacements,_primaryPhysicalElements);
    return elementStrains;
  }

  vector<ElementStrain>
  computeNodalStrains(const vector<ElementVector> & displacements) {
    // TODO THIS NEEDS TO BE UPDATED FOR TWO ELEMENTS
    errorStatement("computeNodalStrains needs to be updated for the two element "
                   "assembler class\n");
    vector<ElementStrain> nodalStrains =
      Assemblers::Utilities::computeNodalStrains(displacements,_primaryPhysicalElements);
    return nodalStrains;
  }

  const vector<PhysicalElementType1> &
  getPrimaryPhysicalElements() const {
    return _primaryPhysicalElements;
  }

  const vector<PhysicalElementType2> &
  getSecondaryPhysicalElements() const {
    return _secondaryPhysicalElements;
  }

  size_t
  getNumberOfNodes() const {
    return _numberOfNodes;
  }

  Eigen::VectorXd
  allocateZeroedGlobalForceVector() const {
    const size_t numberOfDOFs = _numberOfNodes * DegreesOfFreedom;
    Eigen::VectorXd temp(numberOfDOFs);
    temp.fill(0);
    return temp;
  }

  Eigen::SparseMatrix<double>
  allocateZeroedStiffnessMatrix() const {
    const size_t numberOfDOFs = _numberOfNodes * DegreesOfFreedom;
    Eigen::SparseMatrix<double> temp(numberOfDOFs, numberOfDOFs);
    return temp;
  }

protected:
  vector<PhysicalElementType1> _primaryPhysicalElements;
  vector<PhysicalElementType2> _secondaryPhysicalElements;
  size_t _numberOfNodes;
  ExternalForceAssembler _externalForceAssembler;
};

template <class PhysicalElementType1,
          class PhysicalElementType2,
          class PhysicalElementType3,
          class ExternalForceAssembler = ExternalForceAssemblerPlacebo<PhysicalElementType1> >
class AssemblerThreePhysicalElements {

public:
  typedef typename PhysicalElementType1::Vector          ElementVector;
  //typedef typename PhysicalElementType1::Stress          ElementStress;
  //typedef typename PhysicalElementType1::Strain          ElementStrain;
  //typedef typename PhysicalElementType1::TangentMatrix   PrimaryElementTangentMatrix;
  //typedef typename PhysicalElementType1::TangentMatrix   SecondaryElementTangentMatrix;
  typedef typename PhysicalElementType1::StiffnessMatrix PrimaryElementStiffnessMatrix;
  typedef typename PhysicalElementType1::StiffnessMatrix SecondaryElementStiffnessMatrix;
  typedef typename PhysicalElementType1::Point           ElementPoint;
  typedef typename PhysicalElementType1::Node            ElementNode;
  static const unsigned int      PrimaryElementNumberOfNodes = PhysicalElementType1::NumberOfNodes;
  static const unsigned int      SecondaryElementNumberOfNodes = PhysicalElementType2::NumberOfNodes;
  static const unsigned int      ThirdSecondaryElementNumberOfNodes = PhysicalElementType3::NumberOfNodes;
  static const unsigned int      SpatialDimension     = PhysicalElementType1::SpatialDimension;
  static const unsigned int      DegreesOfFreedom     = PhysicalElementType1::DegreesOfFreedom;

  AssemblerThreePhysicalElements(const size_t numberOfNodes) : _numberOfNodes(numberOfNodes) {
  }

  // external force assembler has to handle the general case
  template <class E>
  void
  addElement(const E & e) {
    _externalForceAssembler.addElement(e);
  }

  // specialization for the physical element type
  void
  addElement(const PhysicalElementType1 & e) {
    _primaryPhysicalElements.push_back(e);
  }

  void
  addElement(const PhysicalElementType2 & e) {
    _secondaryPhysicalElements.push_back(e);
  }

  void
  addElement(const PhysicalElementType3 & e) {
    _tertiaryPhysicalElements.push_back(e);
  }

  void
  assembleEnergy(const vector<ElementVector> & displacements,
                 const double time,
                 double * energy) const {
    Assemblers::Utilities::assembleEnergyUtility(displacements, _primaryPhysicalElements,
                                                 time, energy);
    Assemblers::Utilities::assembleEnergyUtility(displacements, _secondaryPhysicalElements,
                                                 time, energy);
    Assemblers::Utilities::assembleEnergyUtility(displacements, _tertiaryPhysicalElements,
                                                 time, energy);
    _externalForceAssembler.assembleEnergy(displacements, time, energy);
  }

  void
  assembleForceVector(const vector<ElementVector> & displacements,
                      const double time,
                      Eigen::VectorXd * forceVector) const {
    Assemblers::Utilities::assembleForceVectorUtility(displacements,
                                                      _primaryPhysicalElements, time,
                                                      forceVector);
    Assemblers::Utilities::assembleForceVectorUtility(displacements,
                                                      _secondaryPhysicalElements, time,
                                                      forceVector);
    Assemblers::Utilities::assembleForceVectorUtility(displacements,
                                                      _tertiaryPhysicalElements, time,
                                                      forceVector);
    _externalForceAssembler.assembleForceVector(displacements, time,
                                                forceVector);
  }

  void
  assembleStiffnessMatrix(const vector<ElementVector> & displacements,
                          const double time,
                          Eigen::SparseMatrix<double> * stiffnessMatrix) const {
    Assemblers::Utilities::assembleStiffnessMatrixUtility(displacements,
                                                          _primaryPhysicalElements, time,
                                                          stiffnessMatrix);
    Assemblers::Utilities::assembleStiffnessMatrixUtility(displacements,
                                                          _secondaryPhysicalElements, time,
                                                          stiffnessMatrix);
    Assemblers::Utilities::assembleStiffnessMatrixUtility(displacements,
                                                          _tertiaryPhysicalElements, time,
                                                          stiffnessMatrix);
    _externalForceAssembler.assembleStiffnessMatrix(displacements, time,
                                                    stiffnessMatrix);
  }

  void
  updateInternalVariables(const vector<ElementVector> & displacements,
                          const double time)  {
    Assemblers::Utilities::updateInternalVariablesUtility(displacements,time,
                                                          &_primaryPhysicalElements);
    Assemblers::Utilities::updateInternalVariablesUtility(displacements,time,
                                                          &_secondaryPhysicalElements);
    Assemblers::Utilities::updateInternalVariablesUtility(displacements,time,
                                                          &_tertiaryPhysicalElements);
  }

  Eigen::SparseMatrix<double>
  assembleLumpedMassMatrix() const {
    size_t numberOfDofs = _numberOfNodes * DegreesOfFreedom;
    Eigen::SparseMatrix<double> massMatrix(numberOfDofs,numberOfDofs);
    Assemblers::Utilities::assembleLumpedMassMatrixUtility(_numberOfNodes,
                                                           _primaryPhysicalElements,
                                                           &massMatrix);
    Assemblers::Utilities::assembleLumpedMassMatrixUtility(_numberOfNodes,
                                                           _secondaryPhysicalElements,
                                                           &massMatrix);
    Assemblers::Utilities::assembleLumpedMassMatrixUtility(_numberOfNodes,
                                                           _tertiaryPhysicalElements,
                                                           &massMatrix);
    return massMatrix;
  }

  Eigen::SparseMatrix<double>
  assembleConsistentMassMatrix() const {
    size_t numberOfDofs = _numberOfNodes * DegreesOfFreedom;
    Eigen::SparseMatrix<double> massMatrix(numberOfDofs, numberOfDofs);
    Assemblers::Utilities::assembleConsistentMassMatrixUtility(_numberOfNodes,
                                                               _primaryPhysicalElements,
                                                               &massMatrix);
    Assemblers::Utilities::assembleConsistentMassMatrixUtility(_numberOfNodes,
                                                               _secondaryPhysicalElements,
                                                               &massMatrix);
    Assemblers::Utilities::assembleConsistentMassMatrixUtility(_numberOfNodes,
                                                               _tertiaryPhysicalElements,
                                                               &massMatrix);
    return massMatrix;
  }

  /*vector<ElementStress>
    computeElementStresses(const vector<ElementVector> & displacements,
    const double time) {
    // TODO THIS NEEDS TO BE UPDATED FOR TWO ELEMENTS
    errorStatement("computeElementStresses needs to be updated for the three element "
    "assembler class\n");
    vector<ElementStress> elementStresses =
    Assemblers::Utilities::computeElementStresses(displacements,_primaryPhysicalElements, time);
    return elementStresses;
    }

    vector<ElementStress>
    computeNodalStresses(const vector<ElementVector> & displacements,
    const double time) {
    // TODO THIS NEEDS TO BE UPDATED FOR TWO ELEMENTS
    errorStatement("computeNodalStresses needs to be updated for the three element "
    "assembler class\n");
    vector<ElementStress> nodalStresses =
    Assemblers::Utilities::computeNodalStresses(displacements,_primaryPhysicalElements, time);
    return nodalStresses;
    }

    vector<ElementStrain>
    computeElementStrains(const vector<ElementVector> & displacements) {
    // TODO THIS NEEDS TO BE UPDATED FOR TWO ELEMENTS
    errorStatement("computeElementStrains needs to be updated for the three element "
    "assembler class\n");
    vector<ElementStrain> elementStrains =
    Assemblers::Utilities::computeElementStrains(displacements,_primaryPhysicalElements);
    return elementStrains;
    }

    vector<ElementStrain>
    computeNodalStrains(const vector<ElementVector> & displacements) {
    // TODO THIS NEEDS TO BE UPDATED FOR TWO ELEMENTS
    errorStatement("computeNodalStrains needs to be updated for the three element "
    "assembler class\n");
    vector<ElementStrain> nodalStrains =
    Assemblers::Utilities::computeNodalStrains(displacements,_primaryPhysicalElements);
    return nodalStrains;
    }*/

  const vector<PhysicalElementType1> &
  getPrimaryPhysicalElements() const {
    return _primaryPhysicalElements;
  }

  const vector<PhysicalElementType2> &
  getSecondaryPhysicalElements() const {
    return _secondaryPhysicalElements;
  }

  const vector<PhysicalElementType3> &
  getTertiaryPhysicalElements() const {
    return _tertiaryPhysicalElements;
  }

  size_t
  getNumberOfNodes() const {
    return _numberOfNodes;
  }

  Eigen::VectorXd
  allocateZeroedGlobalForceVector() const {
    const size_t numberOfDOFs = _numberOfNodes * DegreesOfFreedom;
    Eigen::VectorXd temp(numberOfDOFs);
    temp.fill(0);
    return temp;
  }

  Eigen::SparseMatrix<double>
  allocateZeroedStiffnessMatrix() const {
    const size_t numberOfDOFs = _numberOfNodes * DegreesOfFreedom;
    Eigen::SparseMatrix<double> temp(numberOfDOFs, numberOfDOFs);
    return temp;
  }

protected:
  vector<PhysicalElementType1> _primaryPhysicalElements;
  vector<PhysicalElementType2> _secondaryPhysicalElements;
  vector<PhysicalElementType3> _tertiaryPhysicalElements;
  size_t _numberOfNodes;
  ExternalForceAssembler _externalForceAssembler;
};

template <class PhysicalElementType,
          class ExternalForceAssembler = ExternalForceAssemblerPlacebo<PhysicalElementType> >
class AssemblerWithOrderParameters : public Assembler<PhysicalElementType, ExternalForceAssembler> {

  typedef Assembler<PhysicalElementType, ExternalForceAssembler> Base;
public:
  typedef typename PhysicalElementType::Vector          ElementVector;
  typedef typename PhysicalElementType::OrderParameter  OrderParameter;
  static const unsigned int      OrderParameterDofs   = PhysicalElementType::OrderParameterDofs;

  AssemblerWithOrderParameters(const size_t numberOfNodes) : Base(numberOfNodes) {
  }

  template <class QuadratureRule, class MaterialModel>
  AssemblerWithOrderParameters(const SingleElementMesh<PhysicalElementType> & mesh,
                               const typename PhysicalElementType::Properties & elementProperties,
                               const QuadratureRule & quadratureRule,
                               const MaterialModel & materialModel) :
    Base(mesh, elementProperties, quadratureRule, materialModel) {
  }

  void
  updateOrderParameters(const vector<OrderParameter> & orderParameters)  {
    Assemblers::Utilities::updateOrderParametersUtility(orderParameters,
                                                        &this->_physicalElements);
  }

  void
  assembleOrderParameterForceVector(const vector<ElementVector> & displacements,
                                    const double time,
                                    Eigen::VectorXd * forceVector) const {
    forceVector->fill(0);
    try {
      Assemblers::Utilities::assembleOrderParameterForceVectorUtility(displacements,
                                                                      this->_physicalElements, time,
                                                                      forceVector);
    } catch (std::exception & e) {
      errorStatement("exception caught trying to assemble the force vector for "
                     "the physical elements\n");
      throw;
    }
    try {
      this->_externalForceAssembler.assembleForceVector(displacements, time,
                                                        forceVector);
    } catch (std::exception & e) {
      errorStatement("exception caught trying to assemble the force vector for "
                     "the external force assembler\n");
      throw;
    }
  }

  Eigen::SparseMatrix<double>
  assembleLumpedMobilityMatrix() const {
    const size_t numberOfDofs = this->_numberOfNodes * OrderParameterDofs;
    Eigen::SparseMatrix<double> mobilityMatrix(numberOfDofs,numberOfDofs);
    Assemblers::Utilities::assembleLumpedMobilityMatrixUtility(this->_numberOfNodes,
                                                               this->_physicalElements,
                                                               &mobilityMatrix);
    return mobilityMatrix;
  }

  Eigen::SparseMatrix<double>
  assembleConsistentMobilityMatrix() const {
    const size_t numberOfDofs = this->_numberOfNodes * OrderParameterDofs;
    Eigen::SparseMatrix<double> mobilityMatrix(numberOfDofs, numberOfDofs);
    mobilityMatrix.setZero();
    Assemblers::Utilities::assembleConsistentMobilityMatrixUtility(this->_numberOfNodes,
                                                                   this->_physicalElements,
                                                                   &mobilityMatrix);
    return mobilityMatrix;
  }

};

template <class PhysicalElementType,
          class ExternalForceAssembler = ExternalForceAssemblerPlacebo<PhysicalElementType> >
class ElectroMechanicalAssembler : public Assembler<PhysicalElementType, ExternalForceAssembler> {

  typedef Assembler<PhysicalElementType, ExternalForceAssembler> Base;
public:
  typedef typename PhysicalElementType::Vector          ElementVector;
  static const unsigned int MechanicalDegreesOfFreedom = PhysicalElementType::MechanicalDegreesOfFreedom;

  ElectroMechanicalAssembler(const size_t numberOfNodes) : Base(numberOfNodes) {
  }

  template <class QuadratureRule, class MaterialModel>
  ElectroMechanicalAssembler(const SingleElementMesh<PhysicalElementType> & mesh,
                             const typename PhysicalElementType::Properties & elementProperties,
                             const QuadratureRule & quadratureRule,
                             const MaterialModel & materialModel) :
    Base(mesh, elementProperties, quadratureRule, materialModel) {
  }

  Eigen::SparseMatrix<double>
  assembleLumpedMassMatrix() const {
    size_t numberOfDofs = this->_numberOfNodes * MechanicalDegreesOfFreedom;
    Eigen::SparseMatrix<double> massMatrix(numberOfDofs,numberOfDofs);
    Assemblers::Utilities::assembleMechanicalLumpedMassMatrixUtility(this->_numberOfNodes,
                                                                     this->_physicalElements,
                                                                     &massMatrix);
    return massMatrix;
  }

  Eigen::SparseMatrix<double>
  assembleConsistentMassMatrix() const {
    size_t numberOfDofs = this->_numberOfNodes * MechanicalDegreesOfFreedom;
    Eigen::SparseMatrix<double> massMatrix(numberOfDofs, numberOfDofs);
    massMatrix.setZero();
    Assemblers::Utilities::assembleMechanicalConsistentMassMatrixUtility(this->_numberOfNodes,
                                                                         this->_physicalElements,
                                                                         &massMatrix);
    return massMatrix;
  }

};

}

#endif  // ASSEMBLER_H
