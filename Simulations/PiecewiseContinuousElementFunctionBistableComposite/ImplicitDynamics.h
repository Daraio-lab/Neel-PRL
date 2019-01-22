// -*- C++ -*-
#ifndef IMPLICIT_DYNAMICS_H
#define IMPLICIT_DYNAMICS_H
#include "../core/Definitions.h"
#include "../core/Utilities.h"
#include <ads/timer.h>
#include <Eigen/LU>
#include <Eigen/SuperLUSupport>

namespace Dynamics {

// A structure to store the displacements, velocities, and accelerations
struct State {

  VectorXd _displacements;
  VectorXd _velocities;
  VectorXd _accelerations;

  State() {
  }

  State(const unsigned int totalNumberOfDOFs) :
    _displacements(totalNumberOfDOFs, 1),
    _velocities(totalNumberOfDOFs, 1),
    _accelerations(totalNumberOfDOFs, 1) {
    _displacements.fill(0);
    _velocities.fill(0);
    _accelerations.fill(0);
  }
};

template <class Assembler>
class Newmark {

  typedef typename Assembler::ElementVector                         Vector;

public:

  Newmark(const Assembler & assembler,
          const double timestep,
          const double dampingAlpha,
          const double dampingBeta,
          const double newmarkBeta,
          const double newmarkGamma,
          const VectorXd initialVelocities,
          const unsigned int maxNumberOfIterations = 1000,
          const double iterationTolerance = 1e-8) :
    _assembler(assembler),
    _timestep(timestep),
    _dampingAlpha(dampingAlpha),
    _dampingBeta(dampingBeta),
    _newmarkBeta(newmarkBeta),
    _newmarkGamma(newmarkGamma),
    _maxNumberOfIterations(maxNumberOfIterations),
    _iterationTolerance(iterationTolerance) {

    const unsigned int numberOfNodes = _assembler.getNumberOfNodes();
    const unsigned int totalNumberOfDOFs =
      _assembler.getNumberOfNodes() * Assembler::DegreesOfFreedom;

    if (initialVelocities.size() != numberOfNodes) {
      throwException("Cannot instantiate a newmark integrator with %zu "
                     "initial velocities and %u nodes in the assembler",
                     initialVelocities.size(), numberOfNodes);
    }

    // initialize the state
    _state._displacements.resize(totalNumberOfDOFs, 1);
    _state._displacements.fill(0);
    _state._velocities = initialVelocities;
    _state._accelerations.resize(totalNumberOfDOFs, 1);
    _state._accelerations.fill(0);

  }



  Newmark(const Assembler & assembler,
          const double timestep,
          const double dampingAlpha,
          const double dampingBeta,
          const double newmarkBeta,
          const double newmarkGamma,
          const unsigned int maxNumberOfIterations = 1000,
          const double iterationTolerance = 1e-8,
          const Vector initialVelocity = Vector::Zero()) :
    _assembler(assembler),
    _timestep(timestep),
    _dampingAlpha(dampingAlpha),
    _dampingBeta(dampingBeta),
    _newmarkBeta(newmarkBeta),
    _newmarkGamma(newmarkGamma),
    _maxNumberOfIterations(maxNumberOfIterations),
    _iterationTolerance(iterationTolerance) {

    const unsigned int numberOfNodes = _assembler.getNumberOfNodes();
    const unsigned int totalNumberOfDOFs =
      _assembler.getNumberOfNodes() * Assembler::DegreesOfFreedom;

    // initialize the state
    _state._displacements.resize(totalNumberOfDOFs, 1);
    _state._displacements.fill(0);
    _state._velocities.resize(totalNumberOfDOFs, 1);
    _state._velocities.fill(0);
    for (unsigned int nodeIndex = 0; nodeIndex < numberOfNodes; ++nodeIndex) {
      for (unsigned int dofIndex = 0; dofIndex < Assembler::DegreesOfFreedom;
           ++dofIndex) {
        _state._velocities(nodeIndex * Assembler::DegreesOfFreedom + dofIndex) =
          initialVelocity(dofIndex);
      }
    }
    _state._accelerations.resize(totalNumberOfDOFs, 1);
    _state._accelerations.fill(0);

  }

  Newmark(const Assembler & assembler,
          const double timestep,
          const double dampingAlpha,
          const double dampingBeta,
          const double newmarkBeta,
          const double newmarkGamma,
          const vector<Vector> & initialVelocities,
          const unsigned int maxNumberOfIterations = 1000,
          const double iterationTolerance = 1e-8) :
    _assembler(assembler),
    _timestep(timestep),
    _dampingAlpha(dampingAlpha),
    _dampingBeta(dampingBeta),
    _newmarkBeta(newmarkBeta),
    _newmarkGamma(newmarkGamma),
    _maxNumberOfIterations(maxNumberOfIterations),
    _iterationTolerance(iterationTolerance) {

    const unsigned int numberOfNodes = _assembler.getNumberOfNodes();
    const unsigned int totalNumberOfDOFs =
      _assembler.getNumberOfNodes() * Assembler::DegreesOfFreedom;

    if (initialVelocities.size() != numberOfNodes) {
      throwException("Cannot instantiate a newmark integrator with %zu "
                     "initial velocities and %u nodes in the assembler",
                     initialVelocities.size(), numberOfNodes);
    }

    // initialize the state
    _state._displacements.resize(totalNumberOfDOFs, 1);
    _state._displacements.fill(0);
    _state._velocities = initialVelocities;
    _state._accelerations.resize(totalNumberOfDOFs, 1);
    _state._accelerations.fill(0);

  }

  virtual
  ~Newmark() {
  }

  void
  computeNewmarkUpdate(const vector<EssentialBoundaryCondition> & essentialBCs,
                       const unsigned int timestepIndex,
                       const VerboseFlag verboseFlag = Quiet){
    const vector<EssentialBoundaryCondition> velocityBCs;
    const vector<EssentialBoundaryCondition> accelerationBCs;
    computeNewmarkUpdate(essentialBCs,
                         velocityBCs,
                         accelerationBCs,
                         timestepIndex,
                         verboseFlag);
  }
  

  void
  computeNewmarkUpdate(const vector<EssentialBoundaryCondition> & essentialBCs,
                       const vector<EssentialBoundaryCondition> & velocityBCs,
                       const vector<EssentialBoundaryCondition> & accelerationBCs,
                       const unsigned int timestepIndex,
                       const VerboseFlag verboseFlag = Quiet) {
    ignoreUnusedVariables(velocityBCs,accelerationBCs);

    const double currentTime = timestepIndex * _timestep;

    // Solving for the updated displacements using Newton-Raphson iterations
    if (verboseFlag == Verbose) {
      printf("Implicit Dynamics solver trying to achieve a tolerance of %e in %u "
             "maximum iterations\n", _iterationTolerance, _maxNumberOfIterations);
    }

    typedef typename Assembler::ElementVector ElementVector;
    size_t numberOfDOFs = _state._displacements.size();
    double DegreesOfFreedom = Assembler::DegreesOfFreedom;
    // We start with an initial guess of the current displacements:
    VectorXd solution(numberOfDOFs);
    solution = _state._displacements;
    // Set essential boundary conditions inside the solution vector
    for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
      const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
      const size_t dofIndex = bc._nodeId * DegreesOfFreedom + bc._coordinate;
      solution(dofIndex) = bc._constraint;
    }

    // Creating the vector of nodal-wise displacements
    vector<ElementVector> nodalDisplacements =
      Utilities::distributeGlobalVectorToLocalVectors<Assembler>(solution);
    for (size_t vectorIndex = 0; vectorIndex < nodalDisplacements.size();
         vectorIndex++){
      if (std::isfinite(nodalDisplacements[vectorIndex].norm()) == false){
        throwException("tstep %6u, at beginning of computeNewmarkUpdate "
                       "nodalDisplacements[%lu] is not a finite number",
                       timestepIndex, vectorIndex);
      }
    }

    // Lumped mass matrix
    Eigen::SparseMatrix<double> lumpedMassMatrix =
      _assembler.assembleLumpedMassMatrix();
    if (std::isfinite(lumpedMassMatrix.norm()) == false){
      throwException("tstep %6u, the lumpedMassMatrix computed before doing "
                     "the newmark update has a non-finite number",
                     timestepIndex);
    }

    // Tangent matrix
    Eigen::SparseMatrix<double> tangentMatrix(numberOfDOFs,numberOfDOFs);
    _assembler.assembleStiffnessMatrix(nodalDisplacements, currentTime,
                                       &tangentMatrix);
    if (std::isfinite(tangentMatrix.norm()) == false){
      throwException("tstep %6u, the tangentMatrix computed before doing "
                     "the newmark update has a non-finite number",
                     timestepIndex);
    }

    // Damping matrix
    Eigen::SparseMatrix<double> dampingMatrix =
      _dampingAlpha * lumpedMassMatrix + _dampingBeta*tangentMatrix;

    // Force Vector
    Eigen::VectorXd forceVector(numberOfDOFs);
    forceVector.fill(0);
    _assembler.assembleForceVector(nodalDisplacements, currentTime, &forceVector);
    if (std::isfinite(forceVector.norm()) == false){
      throwException("tstep %6u, the forceVector computed before doing "
                     "the newmark update has a non-finite number",
                     timestepIndex);
    }

    // Computing the "effective force vector"
    Eigen::VectorXd effectiveForceVector = forceVector +
      ( (1. / (_newmarkBeta * pow(_timestep,2.) ) ) * lumpedMassMatrix +
        ( _newmarkGamma / (_newmarkBeta * _timestep) ) * dampingMatrix ) * solution
      - lumpedMassMatrix * ( (1./(_newmarkBeta * pow(_timestep,2.) ) ) *
                             _state._displacements + (1. / (_newmarkBeta * _timestep) ) *
                             _state._velocities + ( (1. / (2. * _newmarkBeta) ) -1.) *
                             _state._accelerations)
      - dampingMatrix * ( (_newmarkGamma / (_newmarkBeta * _timestep) ) *
                          _state._displacements + (_newmarkGamma / _newmarkBeta - 1.) *
                          _state._velocities + _timestep *
                          (_newmarkGamma / (2. * _newmarkBeta) - 1.) * _state._accelerations);

    // Computing the residue vector
    Eigen::VectorXd residueVector = effectiveForceVector;

    // Zero out entries of residue vector that are not considered
    for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
      const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
      const size_t dofIndex = bc._nodeId * DegreesOfFreedom + bc._coordinate;
      effectiveForceVector[dofIndex] = bc._constraint;
      residueVector[dofIndex] = 0.;
    }

    // Computing the residue
    double residue = residueVector.norm();
    if (verboseFlag == Verbose){
      printf("tstep %6u, Initial residue = %9.3e\n", timestepIndex, residue);
    }

    // While the residue > tolerance compute Newton-Raphson iterations
    unsigned int numberOfIterations = 0;
    VectorXd previousIterationDisplacements = _state._displacements;
    Eigen::SparseMatrix<double> effectiveTangentMatrix;
    Eigen::SparseMatrix<double> copyOfEffectiveTangentMatrix;
    while(residue > _iterationTolerance &&
          numberOfIterations < _maxNumberOfIterations) {

      // Computing the "effective tangent matrix"
      effectiveTangentMatrix =
        (1. / (_newmarkBeta * pow(_timestep,2.) ) ) * lumpedMassMatrix +
        (_newmarkGamma / (_newmarkBeta * _timestep) ) * dampingMatrix + tangentMatrix;

      copyOfEffectiveTangentMatrix =
        createCopyOfEffectiveTangentMatrixWithBoundaryConditionsApplied(effectiveTangentMatrix,
                                                                        essentialBCs);

      Eigen::SparseLU<Eigen::SparseMatrix<double> > sparseLUofEffectiveTangentMatrix(copyOfEffectiveTangentMatrix);
      const Eigen::ComputationInfo infoAboutLUdecomposition = sparseLUofEffectiveTangentMatrix.info();
      if (infoAboutLUdecomposition != Eigen::Success){
        errorStatement("LU decomposition failed\n\n");
        printf("The error code is: ");
        if (infoAboutLUdecomposition == Eigen::NumericalIssue){
          printf("Numerical Issue\n");
        }
        else if (infoAboutLUdecomposition == Eigen::NoConvergence){
          printf("No Convergence\n");
        }
        else if (infoAboutLUdecomposition == Eigen::InvalidInput){
          printf("Invalid Input\n");
        }
      }
      else{
        if (verboseFlag == true){
          printf("LU decomposition successfull!\n");
        }
      }
      // Update solution using the Newmark method update rule
      solution -= sparseLUofEffectiveTangentMatrix.solve(effectiveForceVector);
      // Set essential boundary conditions inside the solution vector
      for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
        const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
        const size_t dofIndex = bc._nodeId * DegreesOfFreedom + bc._coordinate;
        solution(dofIndex) = bc._constraint;
      }

      // Reassemble displacements from solution
      nodalDisplacements =
        Utilities::distributeGlobalVectorToLocalVectors<Assembler>(solution);
        
        cout << "NodalDisplacement 1:" << nodalDisplacements[0](0) << endl;
        cout << "NodalDisplacement 2:" << nodalDisplacements[1](0) << endl;
        cout << "NodalDisplacement 3:" << nodalDisplacements[2](0) << endl;

      // Have the problem compute the tangent matrix
      tangentMatrix.setZero();
      _assembler.assembleStiffnessMatrix(nodalDisplacements, currentTime,
                                         &tangentMatrix);

      // Compute the damping matrix
      dampingMatrix =
        _dampingAlpha * lumpedMassMatrix + _dampingBeta * tangentMatrix;

      // Have the problem compute the forces vector
      forceVector.fill(0);
      _assembler.assembleForceVector(nodalDisplacements, currentTime, & forceVector);
      if (std::isfinite(forceVector.norm()) == false){
        throwException("tstep %6u, newmark iteration %u, "
                       "forceVector has a non-finite number",
                       timestepIndex, numberOfIterations);
      }

      // Computing the effective force vector
      effectiveForceVector = forceVector +
        ( (1. / (_newmarkBeta * pow(_timestep,2.) ) ) * lumpedMassMatrix +
          (_newmarkGamma / (_newmarkBeta * _timestep) ) * dampingMatrix) * solution
        - lumpedMassMatrix * ( (1. / (_newmarkBeta * pow(_timestep,2.) ) ) *
                               _state._displacements + (1. / (_newmarkBeta * _timestep)) *
                               _state._velocities + ( (1. / (2. * _newmarkBeta) ) -1.) *
                               _state._accelerations)
        - dampingMatrix * ( (_newmarkGamma / (_newmarkBeta * _timestep) ) *
                            _state._displacements + (_newmarkGamma / _newmarkBeta -1.) *
                            _state._velocities + _timestep *
                            (_newmarkGamma / (2.*_newmarkBeta) -1.) * _state._accelerations);

      // Computing the residue vector
      //residueVector = (previousIterationDisplacements - solution);
      residueVector = effectiveForceVector;

      // zero out entries of effective residue vector that are not considered
      for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
        const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
        const size_t dofIndex = bc._nodeId * DegreesOfFreedom + bc._coordinate;
        effectiveForceVector[dofIndex] = 0;//bc._constraint;
        residueVector[dofIndex] = 0.;
      }

      // Checking the residue
      residue = residueVector.norm()/double(numberOfDOFs);

      if (std::isfinite(residue) == false){
        throwException("tstep %6u, newmark iteration %u, "
                       "residue is not a finite number",
                       timestepIndex, numberOfIterations);
      }

      if (verboseFlag == Verbose) {
        printf("tstep %6u, Newton Raphson iteration %4u, "
               "residue = %8.3e\n", timestepIndex, numberOfIterations, residue);
      }

      numberOfIterations++;
      previousIterationDisplacements = solution;
    }

    if (numberOfIterations == _maxNumberOfIterations) {
      throwException("tstep %6u, Newton Raphson solver could not converge "
                     "in %u iterations.", timestepIndex, _maxNumberOfIterations);
    }
    // Update the state
    State updatedState(numberOfDOFs);
    updatedState._displacements = solution;
    updatedState._accelerations = ( 1./ (_newmarkBeta * pow(_timestep,2.) ) ) *
      (solution - _state._displacements - _timestep*_state._velocities) -
      (1. - 2. * _newmarkBeta) / (2. * _newmarkBeta) * _state._accelerations;
    updatedState._velocities = _state._velocities + _timestep *
      (_newmarkGamma * updatedState._accelerations + (1. - _newmarkGamma) *
       _state._accelerations);

    _state._displacements = updatedState._displacements;
    _state._velocities = updatedState._velocities;
    _state._accelerations = updatedState._accelerations;
  }

  const State &
  getState() const {
    return _state;
  }

  vector<typename Assembler::ElementVector>
  getNodalDisplacements() const {
    return Utilities::distributeGlobalVectorToLocalVectors<Assembler>(_state._displacements);
  }

  string
  getName() const {
    return string("ImplicitDynamics");
  }

protected:
  void
  computeNewmarkUpdateWithConstantMatrices(const vector<EssentialBoundaryCondition> & essentialBCs,
                                           const Eigen::SparseMatrix<double> & lumpedMassMatrix,
                                           const Eigen::SparseMatrix<double> & dampingMatrix,
                                           const Eigen::SparseMatrix<double> & preComputedLumpedAndDampedAddition,
                                           const Eigen::SuperLU<Eigen::SparseMatrix<double> > & lUofEffectiveTangentMatrix,
                                           const unsigned int timestepIndex,
                                           const VerboseFlag verboseFlag = Quiet) {

    const size_t totalNumberOfDOFs = _state._displacements.size();
    const unsigned int DegreesOfFreedom = Assembler::DegreesOfFreedom;
    const double currentTime = timestepIndex * _timestep;

    // We start with an initial guess of the current displacements:
    VectorXd solution = _state._displacements;

    // Set essential boundary conditions inside the solution vector
    for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
      const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
      const size_t dofIndex = bc._nodeId * DegreesOfFreedom + bc._coordinate;
      solution(dofIndex) = bc._constraint;
    }

    // Creating the vector of nodal-wise displacements
    vector<Vector> nodalDisplacements =
      Utilities::distributeGlobalVectorToLocalVectors<Assembler>(solution);
    for (unsigned int nodeIndex = 0; nodeIndex < nodalDisplacements.size();
         ++nodeIndex){
      if (std::isfinite(nodalDisplacements[nodeIndex].norm()) == false){
        throwException("tstep %6u nodalDisplacements[%u] is not a finite number",
                       timestepIndex, nodeIndex);
      }
    }

    // Force Vector
    Eigen::VectorXd forceVector(totalNumberOfDOFs);
    forceVector.fill(0);
    _assembler.assembleForceVector(nodalDisplacements, currentTime, &forceVector);
    if (std::isfinite(forceVector.norm()) == false){
      throwException("tstep %6u initial forceVector has a non-finite number",
                     timestepIndex);
    }

    const MatrixXd preComputedLumpedMultiplication = lumpedMassMatrix *
      ( (1. / (_newmarkBeta * _timestep * _timestep ) ) * _state._displacements +
        (1. / (_newmarkBeta * _timestep)) * _state._velocities +
        (1. / (2. * _newmarkBeta) - 1.) * _state._accelerations) + dampingMatrix *
      ( (_newmarkGamma / (_newmarkBeta * _timestep) ) * _state._displacements +
        (_newmarkGamma / _newmarkBeta - 1.) * _state._velocities +
        _timestep * (_newmarkGamma / (2. * _newmarkBeta) - 1.) * _state._accelerations);

    // Computing the "effective force vector"
    Eigen::VectorXd effectiveForceVector = forceVector +
      preComputedLumpedAndDampedAddition * solution
      - preComputedLumpedMultiplication;

    // Computing the residue vector
    Eigen::VectorXd residueVector = effectiveForceVector;

    // Zero out entries of residue vector that are not considered
    for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
      const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
      const size_t dofIndex = bc._nodeId * DegreesOfFreedom + bc._coordinate;
      effectiveForceVector[dofIndex] = bc._constraint;
      residueVector[dofIndex] = 0.;
    }

    // Computing the residue
    double residue = residueVector.norm();
    if (verboseFlag == Verbose) {
      printf("tstep %6u Initial residue = %e\n", timestepIndex, residue);
    }

    //////////////////////////////////////////////////////////////////
    // ***************** START ITERATIVE SOLVING *******************//
    //VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV//
    // While the residue > tolerance compute Newton-Raphson iterations
    unsigned int iterationIndex = 0;
    while(residue > _iterationTolerance &&
          iterationIndex < _maxNumberOfIterations) {

      // Update solution using the Newmark method update rule
      solution -= lUofEffectiveTangentMatrix.solve(effectiveForceVector);

      // Set essential boundary conditions inside the solution vector
      for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
        const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
        const size_t dofIndex = bc._nodeId * DegreesOfFreedom + bc._coordinate;
        solution(dofIndex) = bc._constraint;
      }

      // Reassemble displacements from solution
      nodalDisplacements =
        Utilities::distributeGlobalVectorToLocalVectors<Assembler>(solution);

      // Have the assembler compute the forces vector
      forceVector.fill(0);
      _assembler.assembleForceVector(nodalDisplacements, currentTime, &forceVector);
      if (std::isfinite(forceVector.norm()) == false) {
        throwException("tstep %6u forceVector has a non-finite number on "
                       "iteration %u", timestepIndex, iterationIndex);
      }

      // Computing the effective force vector
      effectiveForceVector = forceVector +
        preComputedLumpedAndDampedAddition * solution -
        preComputedLumpedMultiplication;

      // Computing the residue vector
      residueVector = effectiveForceVector/forceVector.norm();

      // zero out entries of effective residue vector that are not considered
      for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); ++bcIndex) {
        const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
        const size_t dofIndex = bc._nodeId * DegreesOfFreedom + bc._coordinate;
        effectiveForceVector[dofIndex] = bc._constraint;
        residueVector[dofIndex] = 0.;
      }

      // Checking the residue
      residue = residueVector.norm();
      if (std::isfinite(residue) == false){
        throwException("tstep %6u residue is not a finite number on iteration %u",
                       timestepIndex, iterationIndex);
      }

      if (verboseFlag == Verbose) {
        printf("tstep %6u Newton Raphson iteration %4u, residue = %8.2e\n",
               timestepIndex, iterationIndex, residue);
      }

      iterationIndex++;
    }
    if (iterationIndex == _maxNumberOfIterations) {
      throwException("tstep %6u Newton Raphson solver could not converge in %u "
                     "iterations. Final residual: %e\n", timestepIndex,
                     _maxNumberOfIterations, residue);
    }

    // Update the state
    State updatedState(totalNumberOfDOFs);
    updatedState._displacements = solution;
    updatedState._accelerations = ( 1./ (_newmarkBeta * _timestep * _timestep ) ) *
      (solution - _state._displacements - _timestep*_state._velocities) -
      (1. - 2. * _newmarkBeta) / (2. * _newmarkBeta) * _state._accelerations;
    updatedState._velocities = _state._velocities + _timestep *
      (_newmarkGamma * updatedState._accelerations + (1. - _newmarkGamma) *
       _state._accelerations);
    _state._displacements = updatedState._displacements;
    _state._velocities = updatedState._velocities;
    _state._accelerations = updatedState._accelerations;
  }

  Eigen::SparseMatrix<double>
  createCopyOfEffectiveTangentMatrixWithBoundaryConditionsApplied(const Eigen::SparseMatrix<double> & effectiveTangentMatrix,
                                                                  const vector<EssentialBoundaryCondition> & essentialBCs) {
    Eigen::SparseMatrix<double> copy =
      effectiveTangentMatrix;
    vector<size_t> dofsToConstrain;
    dofsToConstrain.reserve(essentialBCs.size());
    for (size_t bcIndex = 0; bcIndex < essentialBCs.size(); bcIndex++){
      const EssentialBoundaryCondition & bc = essentialBCs[bcIndex];
      const size_t dofIndex = bc._nodeId * Assembler::DegreesOfFreedom + bc._coordinate;
      dofsToConstrain.push_back(dofIndex);
    }
    copy.prune(SparseEigenMatrixRowZeroer(dofsToConstrain));
    // add the ones on the diagonal
    for (size_t vectorIndex = 0; vectorIndex < essentialBCs.size(); vectorIndex++){
      // this is an okay coeffRef
      copy.coeffRef(dofsToConstrain[vectorIndex],
                    dofsToConstrain[vectorIndex]) = 1;
    }
    copy.makeCompressed();
    return copy;
  }



  const Assembler & _assembler;
  const double _timestep;
  const double _dampingAlpha;
  const double _dampingBeta;
  const double _newmarkBeta;
  const double _newmarkGamma;
  const unsigned int _maxNumberOfIterations;
  const double _iterationTolerance;

  State _state;
};

template <class Assembler>
class NewmarkConstantStiffnessAndMass : public Newmark<Assembler> {

  typedef typename Assembler::ElementVector                         Vector;
  typedef Newmark<Assembler>                                        Base;

  void
  precomputeMatrices(const double startingTime) {

    const unsigned int numberOfNodes = this->_assembler.getNumberOfNodes();
    const unsigned int totalNumberOfDOFs =
      this->_assembler.getNumberOfNodes() * Assembler::DegreesOfFreedom;

    const double dampingAlpha = this->_dampingAlpha;
    const double dampingBeta = this->_dampingBeta;
    const double newmarkBeta = this->_newmarkBeta;
    const double newmarkGamma = this->_newmarkGamma;
    const double timestep = this->_timestep;

    // assemble the lumped mass matrix
    _lumpedMassMatrix = this->_assembler.assembleLumpedMassMatrix();

    // assemble the tangent matrix
    Eigen::SparseMatrix<double> tangentMatrix(totalNumberOfDOFs,
                                              totalNumberOfDOFs);
    this->_assembler.assembleStiffnessMatrix(vector<Vector>(numberOfNodes, Vector::Zero()),
                                             startingTime,
                                             &tangentMatrix);

    // form the damping matrix used in newmark updates
    _dampingMatrix =
      dampingAlpha * _lumpedMassMatrix + dampingBeta * tangentMatrix;

    // precompute the lumped and damped sum used in newmark updates
    _preComputedLumpedAndDampedAddition =
      (1. / (newmarkBeta * timestep * timestep ) ) * _lumpedMassMatrix +
      (newmarkGamma / (newmarkBeta * timestep) ) * _dampingMatrix;

    // compute the effective tangent matrix used in newmark updates
    _effectiveTangentMatrix =
      _preComputedLumpedAndDampedAddition + tangentMatrix;

  }

public:

  NewmarkConstantStiffnessAndMass(const Assembler & assembler,
                                  const double startingTime,
                                  const double timestep,
                                  const double dampingAlpha,
                                  const double dampingBeta,
                                  const double newmarkBeta,
                                  const double newmarkGamma,
                                  const unsigned int maxNumberOfIterations = 1000,
                                  const double iterationTolerance = 1e-8,
                                  const Vector initialVelocity = Vector::Zero()) :
    Base(assembler, timestep, dampingAlpha, dampingBeta,
         newmarkBeta, newmarkGamma, maxNumberOfIterations, iterationTolerance,
         initialVelocity) {

    precomputeMatrices(startingTime);

  }

  NewmarkConstantStiffnessAndMass(const Assembler & assembler,
                                  const double startingTime,
                                  const double timestep,
                                  const double dampingAlpha,
                                  const double dampingBeta,
                                  const double newmarkBeta,
                                  const double newmarkGamma,
                                  const VectorXd & initialVelocities,
                                  const unsigned int maxNumberOfIterations = 1000,
                                  const double iterationTolerance = 1e-8) :
    Base(assembler, timestep, dampingAlpha, dampingBeta,
         newmarkBeta, newmarkGamma, initialVelocities, maxNumberOfIterations,
         iterationTolerance) {

    if (initialVelocities.size() != assembler.getNumberOfNodes()) {
      throwException("Cannot instantiate a newmark integrator with %zu "
                     "initial velocities and %u nodes in the assembler",
                     initialVelocities.size(), assembler.getNumberOfNodes());
    }

    precomputeMatrices(startingTime);

  }

  void
  computeNewmarkUpdate(const vector<EssentialBoundaryCondition> & essentialBCs,
                       const unsigned int timestepIndex,
                       const VerboseFlag verboseFlag = Quiet) {

    const Eigen::SparseMatrix<double> copyOfEffectiveTangentMatrixWithBoundaryConditionsApplied =
      this->createCopyOfEffectiveTangentMatrixWithBoundaryConditionsApplied(_effectiveTangentMatrix,
                                                                            essentialBCs);

    const Eigen::SuperLU<Eigen::SparseMatrix<double> > lUofEffectiveTangentMatrix(copyOfEffectiveTangentMatrixWithBoundaryConditionsApplied);
    Base::computeNewmarkUpdateWithConstantMatrices(essentialBCs,
                                                   _lumpedMassMatrix,
                                                   _dampingMatrix,
                                                   _preComputedLumpedAndDampedAddition,
                                                   lUofEffectiveTangentMatrix,
                                                   timestepIndex,
                                                   verboseFlag);
  }

private:
  Eigen::SparseMatrix<double> _dampingMatrix;
  Eigen::SparseMatrix<double> _lumpedMassMatrix;
  Eigen::SparseMatrix<double> _preComputedLumpedAndDampedAddition;
  Eigen::SparseMatrix<double> _effectiveTangentMatrix;

};

template <class Assembler>
class NewmarkConstantStiffnessAndMassAndBoundaryConditions : public Newmark<Assembler> {

  typedef typename Assembler::ElementVector                         Vector;
  typedef Newmark<Assembler>                                        Base;

  void
  precomputeMatrices(const double startingTime) {

    const unsigned int numberOfNodes = this->_assembler.getNumberOfNodes();
    const unsigned int totalNumberOfDOFs =
      this->_assembler.getNumberOfNodes() * Assembler::DegreesOfFreedom;

    const double dampingAlpha = this->_dampingAlpha;
    const double dampingBeta = this->_dampingBeta;
    const double newmarkBeta = this->_newmarkBeta;
    const double newmarkGamma = this->_newmarkGamma;
    const double timestep = this->_timestep;

    // assemble the lumped mass matrix
    _lumpedMassMatrix = this->_assembler.assembleLumpedMassMatrix();

    // assemble the tangent matrix
    Eigen::SparseMatrix<double> tangentMatrix(totalNumberOfDOFs,
                                              totalNumberOfDOFs);
    this->_assembler.assembleStiffnessMatrix(vector<Vector>(numberOfNodes,
                                                            Vector::Zero()),
                                             startingTime,
                                             &tangentMatrix);

    // form the damping matrix used in newmark updates
    _dampingMatrix =
      dampingAlpha * _lumpedMassMatrix + dampingBeta * tangentMatrix;

    // precompute the lumped and damped sum used in newmark updates
    _preComputedLumpedAndDampedAddition =
      (1. / (newmarkBeta * timestep * timestep ) ) * _lumpedMassMatrix +
      (newmarkGamma / (newmarkBeta * timestep) ) * _dampingMatrix;

    // compute the effective tangent matrix used in newmark updates
    const Eigen::SparseMatrix<double> effectiveTangentMatrix =
      _preComputedLumpedAndDampedAddition + tangentMatrix;

    const Eigen::SparseMatrix<double> copyOfEffectiveTangentMatrixWithBoundaryConditionsApplied =
      this->createCopyOfEffectiveTangentMatrixWithBoundaryConditionsApplied(effectiveTangentMatrix,
                                                                            _essentialBCs);

    // ugggg-leee!
    // unfortunately, the constructor of this object is the only way to set
    //  its properties, so we have to do something like this.
    _lUofEffectiveTangentMatrix =
      new Eigen::SuperLU<Eigen::SparseMatrix<double> >(copyOfEffectiveTangentMatrixWithBoundaryConditionsApplied);

  }

public:

  NewmarkConstantStiffnessAndMassAndBoundaryConditions(const Assembler & assembler,
                                                       const vector<EssentialBoundaryCondition> & essentialBCs,
                                                       const double startingTime,
                                                       const double timestep,
                                                       const double dampingAlpha,
                                                       const double dampingBeta,
                                                       const double newmarkBeta,
                                                       const double newmarkGamma,
                                                       const unsigned int maxNumberOfIterations = 1000,
                                                       const double iterationTolerance = 1e-8,
                                                       const Vector initialVelocity = Vector::Zero()) :
    Base(assembler, timestep, dampingAlpha, dampingBeta,
         newmarkBeta, newmarkGamma, maxNumberOfIterations, iterationTolerance,
         initialVelocity),
    _essentialBCs(essentialBCs) {

    precomputeMatrices(startingTime);

  }

  NewmarkConstantStiffnessAndMassAndBoundaryConditions(const Assembler & assembler,
                                                       const vector<EssentialBoundaryCondition> & essentialBCs,
                                                       const double startingTime,
                                                       const double timestep,
                                                       const double dampingAlpha,
                                                       const double dampingBeta,
                                                       const double newmarkBeta,
                                                       const double newmarkGamma,
                                                       const vector<Vector> & initialVelocities,
                                                       const unsigned int maxNumberOfIterations = 1000,
                                                       const double iterationTolerance = 1e-8) :
    Base(assembler, timestep, dampingAlpha, dampingBeta,
         newmarkBeta, newmarkGamma, initialVelocities, maxNumberOfIterations,
         iterationTolerance),
    _essentialBCs(essentialBCs) {

    if (initialVelocities.size() != assembler.getNumberOfNodes()) {
      throwException("Cannot instantiate a newmark integrator with %zu "
                     "initial velocities and %u nodes in the assembler",
                     initialVelocities.size(), assembler.getNumberOfNodes());
    }

    precomputeMatrices(startingTime);

  }

  ~NewmarkConstantStiffnessAndMassAndBoundaryConditions() {
    delete _lUofEffectiveTangentMatrix;
  }

  void
  computeNewmarkUpdate(const unsigned int timestepIndex,
                       const VerboseFlag verboseFlag = Quiet) {

    Base::computeNewmarkUpdateWithConstantMatrices(_essentialBCs,
                                                   _lumpedMassMatrix,
                                                   _dampingMatrix,
                                                   _preComputedLumpedAndDampedAddition,
                                                   *_lUofEffectiveTangentMatrix,
                                                   timestepIndex,
                                                   verboseFlag);
  }

private:
  const vector<EssentialBoundaryCondition> & _essentialBCs;
  Eigen::SparseMatrix<double> _dampingMatrix;
  Eigen::SparseMatrix<double> _lumpedMassMatrix;
  Eigen::SparseMatrix<double> _preComputedLumpedAndDampedAddition;
  const Eigen::SuperLU<Eigen::SparseMatrix<double> > * _lUofEffectiveTangentMatrix;

};

}
#endif // IMPLICIT_DYNAMICS_H
