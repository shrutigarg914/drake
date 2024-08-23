#include "drake/geometry/optimization/iris.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <optional>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "drake/common/symbolic/expression.h"
#include "drake/geometry/optimization/affine_ball.h"
#include "drake/geometry/optimization/cartesian_product.h"
#include "drake/geometry/optimization/convex_set.h"
#include "drake/geometry/optimization/iris_internal.h"
#include "drake/geometry/optimization/minkowski_sum.h"
#include "drake/geometry/optimization/vpolytope.h"
#include "drake/math/autodiff_gradient.h"
#include "drake/multibody/plant/multibody_plant.h"
#include "drake/multibody/rational/rational_forward_kinematics.h"
#include "drake/solvers/choose_best_solver.h"
#include "drake/solvers/ipopt_solver.h"
#include "drake/solvers/snopt_solver.h"

namespace drake {
namespace geometry {
namespace optimization {

using Eigen::MatrixXd;
using Eigen::Ref;
using Eigen::Vector3d;
using Eigen::VectorXd;
using math::RigidTransform;
using multibody::Frame;
using multibody::JacobianWrtVariable;
using multibody::MultibodyPlant;
using solvers::Binding;
using solvers::Constraint;
using solvers::MathematicalProgram;
using symbolic::Expression;
using systems::Context;

HPolyhedron Iris(const ConvexSets& obstacles, const Ref<const VectorXd>& sample,
                 const HPolyhedron& domain, const IrisOptions& options) {
  const int dim = sample.size();
  const int N = obstacles.size();
  DRAKE_DEMAND(domain.ambient_dimension() == dim);
  for (int i = 0; i < N; ++i) {
    DRAKE_DEMAND(obstacles[i]->ambient_dimension() == dim);
  }
  DRAKE_DEMAND(domain.IsBounded());
  const double kEpsilonEllipsoid = 1e-2;
  Hyperellipsoid E = options.starting_ellipse.value_or(
      Hyperellipsoid::MakeHypersphere(kEpsilonEllipsoid, sample));
  HPolyhedron P = domain;
  if (options.termination_func && options.termination_func(P)) {
    throw std::runtime_error(
        "Iris: The options.termination_func() returned true for the initial "
        "region (defined by the domain argument).  Please check the "
        "implementation of your termination_func.");
  }

  if (options.bounding_region) {
    DRAKE_DEMAND(options.bounding_region->ambient_dimension() == dim);
    P = P.Intersection(*options.bounding_region);
  }

  const int num_initial_constraints = P.A().rows();

  // On each iteration, we will build the collision-free polytope represented as
  // {x | A * x <= b}.  Here we pre-allocate matrices of the maximum size.
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A(
      P.A().rows() + N, dim);
  VectorXd b(P.A().rows() + N);
  A.topRows(P.A().rows()) = P.A();
  b.head(P.A().rows()) = P.b();
  // Use pairs {scale, index}, so that I can back out the indices after a sort.
  std::vector<std::pair<double, int>> scaling(N);
  MatrixXd closest_points(dim, N);

  double best_volume = E.Volume();
  int iteration = 0;
  MatrixXd tangent_matrix;

  while (true) {
    DRAKE_ASSERT(best_volume > 0);
    // Find separating hyperplanes
    for (int i = 0; i < N; ++i) {
      const auto touch = E.MinimumUniformScalingToTouch(*obstacles[i]);
      scaling[i].first = touch.first;
      scaling[i].second = i;
      closest_points.col(i) = touch.second;
    }
    std::sort(scaling.begin(), scaling.end());

    int num_constraints = num_initial_constraints;
    tangent_matrix = 2.0 * E.A().transpose() * E.A();
    for (int i = 0; i < N; ++i) {
      // Only add a constraint if this obstacle still has overlap with the set
      // that has been constructed so far on this iteration.
      if (HPolyhedron(A.topRows(num_constraints), b.head(num_constraints))
              .IntersectsWith(*obstacles[scaling[i].second])) {
        // Add the tangent to the (scaled) ellipsoid at this point as a
        // constraint.
        const VectorXd point = closest_points.col(scaling[i].second);
        A.row(num_constraints) =
            (tangent_matrix * (point - E.center())).normalized();
        b[num_constraints] = A.row(num_constraints) * point;
        num_constraints++;
      }
    }

    if (options.require_sample_point_is_contained &&
        ((A.topRows(num_constraints) * sample).array() >=
         b.head(num_constraints).array())
            .any()) {
      break;
    }

    if (options.termination_func &&
        options.termination_func(
            HPolyhedron(A.topRows(num_constraints), b.head(num_constraints)))) {
      break;
    }

    P = HPolyhedron(A.topRows(num_constraints), b.head(num_constraints));

    iteration++;
    if (iteration >= options.iteration_limit) {
      break;
    }

    E = P.MaximumVolumeInscribedEllipsoid();
    const double volume = E.Volume();
    const double delta_volume = volume - best_volume;
    if (delta_volume <= options.termination_threshold) {
      break;
    }
    if (delta_volume / best_volume <= options.relative_termination_threshold) {
      break;
    }
    best_volume = volume;
  }

  return P;
}

namespace {
// Constructs a ConvexSet for each supported Shape and adds it to the set.
class IrisConvexSetMaker final : public ShapeReifier {
 public:
  DRAKE_NO_COPY_NO_MOVE_NO_ASSIGN(IrisConvexSetMaker)

  IrisConvexSetMaker(const QueryObject<double>& query,
                     std::optional<FrameId> reference_frame)
      : query_{query}, reference_frame_{reference_frame} {};

  void set_reference_frame(const FrameId& reference_frame) {
    DRAKE_DEMAND(reference_frame.is_valid());
    *reference_frame_ = reference_frame;
  }

  void set_geometry_id(const GeometryId& geom_id) { geom_id_ = geom_id; }

  using ShapeReifier::ImplementGeometry;

  void ImplementGeometry(const Box&, void* data) {
    DRAKE_DEMAND(geom_id_.is_valid());
    auto& set = *static_cast<copyable_unique_ptr<ConvexSet>*>(data);
    // Note: We choose HPolyhedron over VPolytope here, but the IRIS paper
    // discusses a significant performance improvement using a "least-distance
    // programming" instance from CVXGEN that exploited the VPolytope
    // representation.  So we may wish to revisit this.
    set = std::make_unique<HPolyhedron>(query_, geom_id_, reference_frame_);
  }

  void ImplementGeometry(const Capsule&, void* data) {
    DRAKE_DEMAND(geom_id_.is_valid());
    auto& set = *static_cast<copyable_unique_ptr<ConvexSet>*>(data);
    set = std::make_unique<MinkowskiSum>(query_, geom_id_, reference_frame_);
  }

  void ImplementGeometry(const Cylinder&, void* data) {
    DRAKE_DEMAND(geom_id_.is_valid());
    auto& set = *static_cast<copyable_unique_ptr<ConvexSet>*>(data);
    set =
        std::make_unique<CartesianProduct>(query_, geom_id_, reference_frame_);
  }

  void ImplementGeometry(const Ellipsoid&, void* data) {
    DRAKE_DEMAND(geom_id_.is_valid());
    auto& set = *static_cast<copyable_unique_ptr<ConvexSet>*>(data);
    set = std::make_unique<Hyperellipsoid>(query_, geom_id_, reference_frame_);
  }

  void ImplementGeometry(const HalfSpace&, void* data) {
    DRAKE_DEMAND(geom_id_.is_valid());
    auto& set = *static_cast<copyable_unique_ptr<ConvexSet>*>(data);
    set = std::make_unique<HPolyhedron>(query_, geom_id_, reference_frame_);
  }

  void ImplementGeometry(const Sphere&, void* data) {
    DRAKE_DEMAND(geom_id_.is_valid());
    auto& set = *static_cast<copyable_unique_ptr<ConvexSet>*>(data);
    set = std::make_unique<Hyperellipsoid>(query_, geom_id_, reference_frame_);
  }

  void ImplementGeometry(const Convex&, void* data) {
    DRAKE_DEMAND(geom_id_.is_valid());
    auto& set = *static_cast<copyable_unique_ptr<ConvexSet>*>(data);
    set = std::make_unique<VPolytope>(query_, geom_id_, reference_frame_);
  }

  void ImplementGeometry(const Mesh&, void* data) {
    DRAKE_DEMAND(geom_id_.is_valid());
    auto& set = *static_cast<copyable_unique_ptr<ConvexSet>*>(data);
    set = std::make_unique<VPolytope>(query_, geom_id_, reference_frame_);
  }

 private:
  const QueryObject<double>& query_{};
  std::optional<FrameId> reference_frame_{};
  GeometryId geom_id_{};
};

}  // namespace

ConvexSets MakeIrisObstacles(const QueryObject<double>& query_object,
                             std::optional<FrameId> reference_frame) {
  const SceneGraphInspector<double>& inspector = query_object.inspector();
  const std::vector<GeometryId> geom_ids =
      inspector.GetAllGeometryIds(Role::kProximity);
  ConvexSets sets(geom_ids.size());

  IrisConvexSetMaker maker(query_object, reference_frame);
  int count = 0;
  for (GeometryId geom_id : geom_ids) {
    maker.set_geometry_id(geom_id);
    inspector.GetShape(geom_id).Reify(&maker, &sets[count++]);
  }
  return sets;
}

namespace {

// Takes an abstract cspace variable s, p_AA, and p_BB and enforces that
// p_WA(s) == p_WB(s). We assume that s = f(q) where f is an invertible map.
class SamePointConstraintAbstractCSpace : public Constraint {
 public:
  DRAKE_NO_COPY_NO_MOVE_NO_ASSIGN(SamePointConstraintAbstractCSpace)

  SamePointConstraintAbstractCSpace(const MultibodyPlant<double>* plant,
                                    const Context<double>& context)
      : Constraint(3, plant->num_positions() + 6, Vector3d::Zero(),
                   Vector3d::Zero()),
        plant_(plant),
        context_(plant->CreateDefaultContext()) {
    DRAKE_DEMAND(plant_ != nullptr);
    context_->SetTimeStateAndParametersFrom(context);
  }

  ~SamePointConstraintAbstractCSpace() override {}

  void set_frameA(const multibody::Frame<double>* frame) { frameA_ = frame; }

  void set_frameB(const multibody::Frame<double>* frame) { frameB_ = frame; }

  void EnableSymbolic() {
    if (symbolic_plant_ != nullptr) {
      return;
    }
    symbolic_plant_ = systems::System<double>::ToSymbolic(*plant_);
    symbolic_context_ = symbolic_plant_->CreateDefaultContext();
    symbolic_context_->SetTimeStateAndParametersFrom(*context_);
  }

  virtual inline VectorXd ExtractAbstractCSpaceVariableFromQ(
      const Ref<const VectorXd>& q) const = 0;
  virtual inline VectorXd ExtractQFromAbstractCSpaceVariable(
      const Ref<const VectorXd>& cspace_variable) const = 0;
  virtual inline AutoDiffVecXd ExtractQFromAbstractCSpaceVariable(
      const Ref<const AutoDiffVecXd>& cspace_variable) const = 0;
  virtual inline VectorX<symbolic::Expression>
  ExtractQFromAbstractCSpaceVariable(
      const Ref<const VectorX<symbolic::Expression>>& cspace_variable)
      const = 0;

  const MultibodyPlant<double>* const plant_;
  const multibody::Frame<double>* frameA_{nullptr};
  const multibody::Frame<double>* frameB_{nullptr};
  std::unique_ptr<Context<double>> context_;

  std::unique_ptr<MultibodyPlant<Expression>> symbolic_plant_{nullptr};
  std::unique_ptr<Context<Expression>> symbolic_context_{nullptr};

 protected:
  // Js_v_WF = Jq_v_WF * dqdot/dsdot
  virtual inline void Compute_Js_v_WF(const Ref<const AutoDiffVecXd>& s,
                                      const Ref<const MatrixXd>& Jq_v_WF,
                                      EigenPtr<MatrixXd> Js_v_WF) const = 0;

 private:
  void DoEval(const Ref<const VectorXd>& x, VectorXd* y) const override {
    DRAKE_DEMAND(frameA_ != nullptr);
    DRAKE_DEMAND(frameB_ != nullptr);
    VectorXd s = x.head(plant_->num_positions());
    VectorXd q = ExtractQFromAbstractCSpaceVariable(s);
    Vector3d p_AA = x.template segment<3>(plant_->num_positions()),
             p_BB = x.template tail<3>();
    Vector3d p_WA, p_WB;
    plant_->SetPositions(context_.get(), q);
    plant_->CalcPointsPositions(*context_, *frameA_, p_AA,
                                plant_->world_frame(), &p_WA);
    plant_->CalcPointsPositions(*context_, *frameB_, p_BB,
                                plant_->world_frame(), &p_WB);
    *y = p_WA - p_WB;
  }

  // p_WA = X_WA(s)*p_AA
  // dp_WA = Js_v_WA*ds + X_WA(s)*dp_AA
  void DoEval(const Ref<const AutoDiffVecXd>& x,
              AutoDiffVecXd* y) const override {
    DRAKE_DEMAND(frameA_ != nullptr);
    DRAKE_DEMAND(frameB_ != nullptr);
    AutoDiffVecXd s = x.head(plant_->num_positions());
    // discard the gradient as the gradients of s are not needed to compute
    // Js_v_WA, Js_v_WB
    VectorXd q = ExtractQFromAbstractCSpaceVariable(ExtractDoubleOrThrow(s));
    Vector3<AutoDiffXd> p_AA = x.template segment<3>(plant_->num_positions()),
                        p_BB = x.template tail<3>();

    plant_->SetPositions(context_.get(), q);

    const RigidTransform<double>& X_WA =
        plant_->EvalBodyPoseInWorld(*context_, frameA_->body());
    const RigidTransform<double>& X_WB =
        plant_->EvalBodyPoseInWorld(*context_, frameB_->body());

    Eigen::Matrix3Xd Jq_v_WA(3, plant_->num_positions()),
        Jq_v_WB(3, plant_->num_positions());
    plant_->CalcJacobianTranslationalVelocity(
        *context_, JacobianWrtVariable::kQDot, *frameA_,
        ExtractDoubleOrThrow(p_AA), plant_->world_frame(),
        plant_->world_frame(), &Jq_v_WA);
    plant_->CalcJacobianTranslationalVelocity(
        *context_, JacobianWrtVariable::kQDot, *frameB_,
        ExtractDoubleOrThrow(p_BB), plant_->world_frame(),
        plant_->world_frame(), &Jq_v_WB);

    Eigen::Matrix3Xd Js_v_WA(3, plant_->num_positions()),
        Js_v_WB(3, plant_->num_positions());
    Eigen::MatrixXd dqdot_ds =
        MatrixXd::Identity(plant_->num_positions(), plant_->num_positions());

    Compute_Js_v_WF(s, Jq_v_WA, &Js_v_WA);
    Compute_Js_v_WF(s, Jq_v_WB, &Js_v_WB);

    const Eigen::Vector3d y_val =
        X_WA * math::ExtractValue(p_AA) - X_WB * math::ExtractValue(p_BB);
    Eigen::Matrix3Xd dy(3, plant_->num_positions() + 6);
    dy << Js_v_WA - Js_v_WB, X_WA.rotation().matrix(),
        -X_WB.rotation().matrix();
    *y = math::InitializeAutoDiff(y_val, dy * math::ExtractGradient(x));
  }

  void DoEval(const Ref<const VectorX<symbolic::Variable>>& x,
              VectorX<symbolic::Expression>* y) const override {
    DRAKE_DEMAND(symbolic_plant_ != nullptr);
    DRAKE_DEMAND(frameA_ != nullptr);
    DRAKE_DEMAND(frameB_ != nullptr);
    const Frame<Expression>& frameA =
        symbolic_plant_->get_frame(frameA_->index());
    const Frame<Expression>& frameB =
        symbolic_plant_->get_frame(frameB_->index());

    VectorX<Expression> s = x.head(plant_->num_positions());
    VectorX<Expression> q = ExtractQFromAbstractCSpaceVariable(s);
    Vector3<Expression> p_AA = x.template segment<3>(plant_->num_positions()),
                        p_BB = x.template tail<3>();
    Vector3<Expression> p_WA, p_WB;
    symbolic_plant_->SetPositions(symbolic_context_.get(), q);
    symbolic_plant_->CalcPointsPositions(*symbolic_context_, frameA, p_AA,
                                         symbolic_plant_->world_frame(), &p_WA);
    symbolic_plant_->CalcPointsPositions(*symbolic_context_, frameB, p_BB,
                                         symbolic_plant_->world_frame(), &p_WB);
    *y = p_WA - p_WB;
  }
};

// Takes q, p_AA, and p_BB and enforces that p_WA == p_WB.
class SamePointConstraint : public SamePointConstraintAbstractCSpace {
 public:
  SamePointConstraint(const MultibodyPlant<double>* plant,
                      const Context<double>& context)
      : SamePointConstraintAbstractCSpace(plant, context) {}

  VectorXd ExtractAbstractCSpaceVariableFromQ(
      const Ref<const VectorXd>& q) const override {
    return q;
  }

  VectorXd ExtractQFromAbstractCSpaceVariable(
      const Ref<const VectorXd>& cspace_variable) const override {
    return cspace_variable;
  }
  AutoDiffVecXd ExtractQFromAbstractCSpaceVariable(
      const Ref<const AutoDiffVecXd>& cspace_variable) const override {
    return cspace_variable;
  }
  VectorX<symbolic::Expression> ExtractQFromAbstractCSpaceVariable(
      const Ref<const VectorX<symbolic::Expression>>& cspace_variable)
      const override {
    return cspace_variable;
  }

 protected:
  inline void Compute_Js_v_WF(const Ref<const AutoDiffVecXd>&,
                              const Ref<const MatrixXd>& Jq_v_WF,
                              EigenPtr<MatrixXd> Js_v_WF) const override {
    *Js_v_WF = Eigen::Map<const Eigen::MatrixXd>(Jq_v_WF.data(), Jq_v_WF.rows(),
                                                 Jq_v_WF.cols());
  }
};

// takes s, p_AA,and p_BB and enforces that p_WA == p_WB
class SamePointConstraintRational : public SamePointConstraintAbstractCSpace {
 public:
  SamePointConstraintRational(const multibody::RationalForwardKinematics*
                                  rational_forward_kinematics_ptr,
                              const Ref<const VectorXd>& q_star,
                              const Context<double>& context)
      : SamePointConstraintAbstractCSpace(
            &rational_forward_kinematics_ptr->plant(), context),
        rational_forward_kinematics_ptr_(rational_forward_kinematics_ptr),
        q_star_(q_star) {}

  VectorXd ExtractAbstractCSpaceVariableFromQ(
      const Ref<const VectorXd>& q) const override {
    return rational_forward_kinematics_ptr_->ComputeSValue(q, q_star_);
  }
  VectorXd ExtractQFromAbstractCSpaceVariable(
      const Ref<const VectorXd>& cspace_variable) const override {
    return rational_forward_kinematics_ptr_->ComputeQValue(cspace_variable,
                                                           q_star_);
  }
  AutoDiffVecXd ExtractQFromAbstractCSpaceVariable(
      const Ref<const AutoDiffVecXd>& cspace_variable) const override {
    return rational_forward_kinematics_ptr_->ComputeQValue(cspace_variable,
                                                           q_star_);
  }
  VectorX<symbolic::Expression> ExtractQFromAbstractCSpaceVariable(
      const Ref<const VectorX<symbolic::Expression>>& cspace_variable)
      const override {
    return rational_forward_kinematics_ptr_->ComputeQValue(cspace_variable,
                                                           q_star_);
  }

 protected:
  inline void Compute_Js_v_WF(const Ref<const AutoDiffVecXd>& s,
                              const Ref<const MatrixXd>& Jq_v_WF,
                              EigenPtr<MatrixXd> Js_v_WF) const override {
    // If q = g(s) then qdot = dg/ds * ds/dt and so dqdot/ds_dot = dg/ds.
    // To get dg/ds, we compute q after setting the gradient of s to the
    // identity.
    AutoDiffVecXd s_identity_derivative = math::InitializeAutoDiff(
        math::ExtractValue(s), Eigen::MatrixXd::Identity(s.rows(), s.rows()));
    AutoDiffVecXd q = ExtractQFromAbstractCSpaceVariable(s_identity_derivative);
    Eigen::MatrixXd dqdot_dsdot =
        MatrixXd::Identity(plant_->num_positions(), plant_->num_positions());
    // since sᵢ = f(qᵢ) then we only need to operate on the diagonal.
    for (int i = 0; i < static_cast<int>(dqdot_dsdot.rows()); ++i) {
      dqdot_dsdot(i, i) = q(i).derivatives()(i);
    }
    *Js_v_WF = Jq_v_WF * dqdot_dsdot;
  }

 private:
  const multibody::RationalForwardKinematics* rational_forward_kinematics_ptr_;
  const VectorXd q_star_;
};

// Defines a MathematicalProgram to solve the problem
// min_q (q-d) CᵀC (q-d)
// s.t. setA in frameA and setB in frameB are in collision in q.
//      Aq ≤ b.
// where C, d are the matrix and center from the hyperellipsoid E.
//
// The class design supports repeated solutions of the (nearly) identical
// problem from different initial guesses.
class ClosestCollisionProgram {
 public:
  ClosestCollisionProgram(
      std::shared_ptr<SamePointConstraintAbstractCSpace> same_point_constraint,
      const multibody::Frame<double>& frameA,
      const multibody::Frame<double>& frameB, const ConvexSet& setA,
      const ConvexSet& setB, const Hyperellipsoid& E,
      const Ref<const Eigen::MatrixXd>& A, const Ref<const VectorXd>& b) {
    q_ = prog_.NewContinuousVariables(A.cols(), "q");

    P_constraint_ = prog_.AddLinearConstraint(
        A,
        VectorXd::Constant(b.size(), -std::numeric_limits<double>::infinity()),
        b, q_);

    // Scale the objective so the eigenvalues are close to 1, using
    // scale*lambda_min = 1/scale*lambda_max.
    const MatrixXd Asq = E.A().transpose() * E.A();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Asq);
    const double scale = 1.0 / std::sqrt(es.eigenvalues().maxCoeff() *
                                         es.eigenvalues().minCoeff());
    prog_.AddQuadraticErrorCost(scale * Asq, E.center(), q_);

    auto p_AA = prog_.NewContinuousVariables<3>("p_AA");
    auto p_BB = prog_.NewContinuousVariables<3>("p_BB");
    setA.AddPointInSetConstraints(&prog_, p_AA);
    setB.AddPointInSetConstraints(&prog_, p_BB);

    same_point_constraint->set_frameA(&frameA);
    same_point_constraint->set_frameB(&frameB);
    prog_.AddConstraint(same_point_constraint, {q_, p_AA, p_BB});

    // Help nonlinear optimizers (e.g. SNOPT) avoid trivial local minima at the
    // origin.
    prog_.SetInitialGuess(p_AA, Vector3d::Constant(.01));
    prog_.SetInitialGuess(p_BB, Vector3d::Constant(.01));
  }

  void UpdatePolytope(const Ref<const Eigen::MatrixXd>& A,
                      const Ref<const VectorXd>& b) {
    P_constraint_->evaluator()->UpdateCoefficients(
        A,
        VectorXd::Constant(b.size(), -std::numeric_limits<double>::infinity()),
        b);
  }

  // Returns true iff a collision is found.
  // Sets `closest` to an optimizing solution q*, if a solution is found.
  bool Solve(const solvers::SolverInterface& solver,
             const Eigen::Ref<const VectorXd>& q_guess, VectorXd* closest) {
    prog_.SetInitialGuess(q_, q_guess);
    solvers::MathematicalProgramResult result;
    solver.Solve(prog_, std::nullopt, std::nullopt, &result);
    if (result.is_success()) {
      *closest = result.GetSolution(q_);
      return true;
    }
    return false;
  }

    bool Solve(const solvers::SolverInterface& solver,
             const Eigen::Ref<const VectorXd>& q_guess, 
             const std::optional<solvers::SolverOptions>& solver_options,
             VectorXd* closest) {
    prog_.SetInitialGuess(q_, q_guess);
    solvers::MathematicalProgramResult result;
    solver.Solve(prog_, std::nullopt, solver_options, &result);
    if (result.is_success()) {
      *closest = result.GetSolution(q_);
      return true;
    }
    return false;
  }

 private:
  MathematicalProgram prog_;
  solvers::VectorXDecisionVariable q_;
  std::optional<Binding<solvers::LinearConstraint>> P_constraint_{};
};

// Takes a constraint bound to another mathematical program and defines a new
// constraint that is the negation of one index and one (lower/upper) bound.
class CounterExampleConstraint : public Constraint {
 public:
  DRAKE_NO_COPY_NO_MOVE_NO_ASSIGN(CounterExampleConstraint)

  explicit CounterExampleConstraint(const MathematicalProgram* prog)
      : Constraint(1, prog->num_vars(),
                   Vector1d(-std::numeric_limits<double>::infinity()),
                   Vector1d::Constant(-kSolverConstraintTolerance - 1e-14)),
        prog_{prog} {
    DRAKE_DEMAND(prog != nullptr);
  }

  ~CounterExampleConstraint() = default;

  // Sets the actual constraint to be falsified, overwriting any previously set
  // constraints. The Binding<Constraint> must remain valid for the lifetime of
  // this object (or until a new Binding<Constraint> is set).
  void set(const Binding<Constraint>* binding_with_constraint_to_be_falsified,
           int index, bool falsify_lower_bound) {
    DRAKE_DEMAND(binding_with_constraint_to_be_falsified != nullptr);
    const int N =
        binding_with_constraint_to_be_falsified->evaluator()->num_constraints();
    DRAKE_DEMAND(index >= 0 && index < N);
    binding_ = binding_with_constraint_to_be_falsified;
    index_ = index;
    falsify_lower_bound_ = falsify_lower_bound;
  }

 private:
  void DoEval(const Ref<const VectorXd>& x, VectorXd* y) const override {
    DRAKE_DEMAND(binding_ != nullptr);
    const double val = prog_->EvalBinding(*binding_, x)[index_];
    if (falsify_lower_bound_) {
      // val - lb <= -kSolverConstraintTolerance < 0.
      (*y)[0] = val - binding_->evaluator()->lower_bound()[index_];
    } else {
      // ub - val <= -kSolverConstraintTolerance < 0.
      (*y)[0] = binding_->evaluator()->upper_bound()[index_] - val;
    }
  }

  void DoEval(const Ref<const AutoDiffVecXd>& x,
              AutoDiffVecXd* y) const override {
    DRAKE_DEMAND(binding_ != nullptr);
    const AutoDiffXd val = prog_->EvalBinding(*binding_, x)[index_];
    if (falsify_lower_bound_) {
      // val - lb <= -kSolverConstraintTolerance < 0.
      (*y)[0] = val - binding_->evaluator()->lower_bound()[index_];
    } else {
      // ub - val <= -kSolverConstraintTolerance < 0.
      (*y)[0] = binding_->evaluator()->upper_bound()[index_] - val;
    }
  }

  void DoEval(const Ref<const VectorX<symbolic::Variable>>&,
              VectorX<symbolic::Expression>*) const override {
    // MathematicalProgram::EvalBinding doesn't support symbolic, and we
    // shouldn't get here.
    throw std::logic_error(
        "CounterExampleConstraint doesn't support DoEval for symbolic.");
  }

  const MathematicalProgram* prog_{};
  const Binding<Constraint>* binding_{};
  int index_{0};
  bool falsify_lower_bound_{true};

  // To find a counter-example for a constraints,
  //  g(x) ≤ ub,
  // we need to ask the solver to find
  //  g(x) + kSolverConstraintTolerance > ub,
  // which we implement as
  //  g(x) + kSolverConstraintTolerance ≥ ub + eps.
  // The variable is static so that it is initialized by the time it is accessed
  // in the initializer list of the constructor.
  // TODO(russt): We need a more robust way to get this from the solver. This
  // value works for SNOPT and is reasonable for most solvers.
  static constexpr double kSolverConstraintTolerance{1e-6};
};

// Defines a MathematicalProgram to solve the problem
// min_q (q-d)*CᵀC(q-d)
// s.t. counter-example-constraint
//      Aq ≤ b.
// where C, d are the matrix and center from the hyperellipsoid E.
//
// The class design supports repeated solutions of the (nearly) identical
// problem from different initial guesses.
class CounterExampleProgram {
 public:
  CounterExampleProgram(
      std::shared_ptr<CounterExampleConstraint> counter_example_constraint,
      const Hyperellipsoid& E, const Ref<const Eigen::MatrixXd>& A,
      const Ref<const VectorXd>& b) {
    q_ = prog_.NewContinuousVariables(A.cols(), "q");

    P_constraint_ = prog_.AddLinearConstraint(
        A,
        VectorXd::Constant(b.size(), -std::numeric_limits<double>::infinity()),
        b, q_);
    // Scale the objective so the eigenvalues are close to 1, using
    // scale*lambda_min = 1/scale*lambda_max.
    const MatrixXd Asq = E.A().transpose() * E.A();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Asq);
    const double scale = 1.0 / std::sqrt(es.eigenvalues().maxCoeff() *
                                         es.eigenvalues().minCoeff());
    prog_.AddQuadraticErrorCost(scale * Asq, E.center(), q_);

    prog_.AddConstraint(counter_example_constraint, q_);
  }

  void UpdatePolytope(const Ref<const Eigen::MatrixXd>& A,
                      const Ref<const VectorXd>& b) {
    P_constraint_->evaluator()->UpdateCoefficients(
        A,
        VectorXd::Constant(b.size(), -std::numeric_limits<double>::infinity()),
        b);
  }

  // Returns true iff a counter-example is found.
  // Sets `closest` to an optimizing solution q*, if a solution is found.
  bool Solve(const solvers::SolverInterface& solver,
             const Eigen::Ref<const VectorXd>& q_guess, 
             const std::optional<solvers::SolverOptions>& solver_options,
             VectorXd* closest) {
    prog_.SetInitialGuess(q_, q_guess);
    solvers::MathematicalProgramResult result;
    solver.Solve(prog_, std::nullopt, solver_options, &result);
    if (result.is_success()) {
      *closest = result.GetSolution(q_);
      return true;
    }
    return false;
  }

 private:
  MathematicalProgram prog_;
  solvers::VectorXDecisionVariable q_;
  std::optional<Binding<solvers::LinearConstraint>> P_constraint_{};
};

// Add the tangent to the (scaled) ellipsoid at @p point as a
// constraint.
void AddTangentToPolytope(
    const Hyperellipsoid& E, const Eigen::Ref<const VectorXd>& point,
    double configuration_space_margin,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>* A,
    Eigen::VectorXd* b, int* num_constraints) {
  while (*num_constraints >= A->rows()) {
    // Increase pre-allocated polytope size.
    A->conservativeResize(A->rows() * 2, A->cols());
    b->conservativeResize(b->rows() * 2);
  }

  A->row(*num_constraints) =
      (E.A().transpose() * E.A() * (point - E.center())).normalized();
  (*b)[*num_constraints] =
      A->row(*num_constraints) * point - configuration_space_margin;
  if (A->row(*num_constraints) * E.center() > (*b)[*num_constraints]) {
    throw std::logic_error(
        "The current center of the IRIS region is within "
        "options.configuration_space_margin of being infeasible.  Check your "
        "sample point and/or any additional constraints you've passed in via "
        "the options. The configuration space surrounding the sample point "
        "must have an interior.");
  }
  *num_constraints += 1;
}

void MakeGuessFeasible(const HPolyhedron& P, VectorXd* guess) {
  const auto& A = P.A();
  const auto& b = P.b();
  const int M = A.rows();
  // Add kEps below because we want to be strictly on the correct side of the
  // inequality (and robust to floating point errors).
  const double kEps = 1e-14;
  // Try projecting the guess onto any violated constraints, one-by-one.
  for (int i = M - 1; i >= 0; --i) {
    if (A.row(i) * *guess - b[i] > 0) {
      // guess = argmin_x ||x - guess||^2 s.t. A.row(i) * x = b(i).
      *guess -=
          A.row(i).normalized().transpose() * (A.row(i) * *guess - b[i] + kEps);
    }
  }
  // If this causes a different constraint to be violated, then just return the
  // Chebyshev center.
  if (!P.PointInSet(*guess)) {
    // Note: This can throw if the set becomes empty. But it should not happen;
    // we check that the seed is feasible in AddTangentToPolytope.
    *guess = P.ChebyshevCenter();
  }
}

struct GeometryPairWithDistance {
  GeometryId geomA;
  GeometryId geomB;
  double distance;

  GeometryPairWithDistance(GeometryId gA, GeometryId gB, double dist)
      : geomA(gA), geomB(gB), distance(dist) {}

  bool operator<(const GeometryPairWithDistance& other) const {
    return distance < other.distance;
  }
};

bool CheckTerminate(const IrisOptions& options, const HPolyhedron& P,
                    const std::string& error_msg, const std::string& info_msg,
                    const bool is_initial_region) {
  if (options.termination_func && options.termination_func(P)) {
    if (is_initial_region) {
      throw std::runtime_error(error_msg);
    }
    log()->info(info_msg);
    return true;
  }
  return false;
}

//}  // namespace

//  Method for actually running Iris. We assume that the HPolyhedron in
//  P_ptr has been set to a finite bounding box defining the joint limits of the
//  plant.
//  @param P_ptr must have its ambient dimension be the same size as the number
//  of joints in the plant.
//  @param E_ptr must have its ambient dimension be the same size as the nubmer
//  of joints in the plant.
//  @return
HPolyhedron _DoIris_(const MultibodyPlant<double>& plant,
                     const Context<double>& context, const IrisOptions& options,
                     const VectorXd& sample,
                     const std::shared_ptr<SamePointConstraintAbstractCSpace>&
                         same_point_constraint) {
  // Check the inputs.
  plant.ValidateContext(context);
  const int nq = plant.num_positions();
  const VectorXd q_sample =
      same_point_constraint->ExtractQFromAbstractCSpaceVariable(sample);
  const int nc = static_cast<int>(options.configuration_obstacles.size());
  // Note: We require finite joint limits to define the bounding box for the
  // IRIS algorithm.
  DRAKE_DEMAND(plant.GetPositionLowerLimits().array().isFinite().all());
  DRAKE_DEMAND(plant.GetPositionUpperLimits().array().isFinite().all());
  DRAKE_DEMAND(options.num_collision_infeasible_samples >= 0);
  for (int i = 0; i < nc; ++i) {
    DRAKE_DEMAND(options.configuration_obstacles[i]->ambient_dimension() == nq);
    if (options.configuration_obstacles[i]->PointInSet(sample)) {
      throw std::runtime_error(
          fmt::format("The seed point is in configuration obstacle {}", i));
    }
  }

  if (options.prog_with_additional_constraints) {
    DRAKE_DEMAND(options.prog_with_additional_constraints->num_vars() == nq);
    DRAKE_DEMAND(options.num_additional_constraint_infeasible_samples >= 0);
  }

  // Make the polytope and ellipsoid.
  HPolyhedron P = HPolyhedron::MakeBox(
      same_point_constraint->ExtractAbstractCSpaceVariableFromQ(
          plant.GetPositionLowerLimits()),
      same_point_constraint->ExtractAbstractCSpaceVariableFromQ(
          plant.GetPositionUpperLimits()));
  DRAKE_DEMAND(P.A().rows() == 2 * nq);

  if (options.bounding_region) {
    DRAKE_DEMAND(options.bounding_region->ambient_dimension() == nq);
    P = P.Intersection(*options.bounding_region);
  }

  const double kEpsilonEllipsoid = 1e-2;
  Hyperellipsoid E = options.starting_ellipse.value_or(
      Hyperellipsoid::MakeHypersphere(kEpsilonEllipsoid, sample));

  // Make all of the convex sets and supporting quantities.
  auto query_object =
      plant.get_geometry_query_input_port().Eval<QueryObject<double>>(context);
  const SceneGraphInspector<double>& inspector = query_object.inspector();
  IrisConvexSetMaker maker(query_object, inspector.world_frame_id());
  std::unordered_map<GeometryId, copyable_unique_ptr<ConvexSet>> sets{};
  std::unordered_map<GeometryId, const multibody::Frame<double>*> frames{};
  const std::vector<GeometryId> geom_ids =
      inspector.GetAllGeometryIds(Role::kProximity);
  copyable_unique_ptr<ConvexSet> temp_set;
  for (GeometryId geom_id : geom_ids) {
    // Make all sets in the local geometry frame.
    FrameId frame_id = inspector.GetFrameId(geom_id);
    maker.set_reference_frame(frame_id);
    maker.set_geometry_id(geom_id);
    inspector.GetShape(geom_id).Reify(&maker, &temp_set);
    sets.emplace(geom_id, std::move(temp_set));
    frames.emplace(geom_id, &plant.GetBodyFromFrameId(frame_id)->body_frame());
  }

  auto pairs = inspector.GetCollisionCandidates();
  const int n = static_cast<int>(pairs.size());
  std::map<std::pair<GeometryId, GeometryId>, std::vector<VectorXd>>
      counter_examples;

  // As a surrogate for the true objective, the pairs are sorted by the distance
  // between each collision pair from the seed point configuration. This could
  // improve computation times and produce regions with fewer faces.
  std::vector<GeometryPairWithDistance> sorted_pairs;
  for (const auto& [geomA, geomB] : pairs) {
    const double distance =
        query_object.ComputeSignedDistancePairClosestPoints(geomA, geomB)
            .distance;
    if (distance < 0.0) {
      throw std::runtime_error(
          fmt::format("The seed point is in collision; geometry {} is in "
                      "collision with geometry {}",
                      inspector.GetName(geomA), inspector.GetName(geomB)));
    }
    sorted_pairs.emplace_back(geomA, geomB, distance);
  }
  std::sort(sorted_pairs.begin(), sorted_pairs.end());

  // On each iteration, we will build the collision-free polytope represented as
  // {x | A * x <= b}.  Here we pre-allocate matrices with a generous maximum
  // size.
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A(
      P.A().rows() + 2 * n + nc, nq);
  VectorXd b(P.A().rows() + 2 * n + nc);
  A.topRows(P.A().rows()) = P.A();
  b.head(P.A().rows()) = P.b();
  int num_initial_constraints = P.A().rows();

  std::shared_ptr<CounterExampleConstraint> counter_example_constraint{};
  std::unique_ptr<CounterExampleProgram> counter_example_prog{};
  std::vector<Binding<Constraint>> additional_constraint_bindings{};
  if (options.prog_with_additional_constraints) {
    counter_example_constraint = std::make_shared<CounterExampleConstraint>(
        options.prog_with_additional_constraints);
    additional_constraint_bindings =
        options.prog_with_additional_constraints->GetAllConstraints();
    // Fail fast if the seed point is infeasible.
    {
      if (!options.prog_with_additional_constraints->CheckSatisfied(
              additional_constraint_bindings, sample)) {
        throw std::runtime_error(
            "options.prog_with_additional_constraints is infeasible at the "
            "seed point. The seed point must be feasible.");
      }
    }
    // Handle bounding box and linear constraints as a special case (extracting
    // them from the additional_constraint_bindings).
    auto AddConstraint = [&](const MatrixXd& new_A, const VectorXd& new_b,
                             const solvers::VectorXDecisionVariable& vars) {
      while (num_initial_constraints + new_A.rows() >= A.rows()) {
        // Increase pre-allocated polytope size.
        A.conservativeResize(A.rows() * 2, A.cols());
        b.conservativeResize(b.rows() * 2);
      }
      for (int i = 0; i < new_b.rows(); ++i) {
        if (!std::isinf(new_b[i])) {
          A.row(num_initial_constraints).setZero();
          for (int j = 0; j < vars.rows(); ++j) {
            const int index = options.prog_with_additional_constraints
                                  ->FindDecisionVariableIndex(vars[j]);
            A(num_initial_constraints, index) = new_A(i, j);
          }
          b[num_initial_constraints++] = new_b[i];
        }
      }
    };
    auto HandleLinearConstraints = [&](const auto& bindings) {
      for (const auto& binding : bindings) {
        AddConstraint(binding.evaluator()->get_sparse_A(),
                      binding.evaluator()->upper_bound(), binding.variables());
        AddConstraint(-binding.evaluator()->get_sparse_A(),
                      -binding.evaluator()->lower_bound(), binding.variables());
        auto pos = std::find(additional_constraint_bindings.begin(),
                             additional_constraint_bindings.end(), binding);
        DRAKE_ASSERT(pos != additional_constraint_bindings.end());
        additional_constraint_bindings.erase(pos);
      }
    };
    HandleLinearConstraints(
        options.prog_with_additional_constraints->bounding_box_constraints());
    HandleLinearConstraints(
        options.prog_with_additional_constraints->linear_constraints());
    counter_example_prog = std::make_unique<CounterExampleProgram>(
        counter_example_constraint, E, A.topRows(num_initial_constraints),
        b.head(num_initial_constraints));

    P = HPolyhedron(A.topRows(num_initial_constraints),
                    b.head(num_initial_constraints));
  }

  if (options.termination_func && options.termination_func(P)) {
    throw std::runtime_error(
        "IrisInConfigurationSpace: The options.termination_func() returned "
        "true for the initial region (defined by the linear constraints in "
        "prog_with_additional_constraints and bounding_region arguments).  "
        "Please check the implementation of your termination_func.");
  }

  DRAKE_THROW_UNLESS(P.PointInSet(sample, 1e-12));

  double best_volume = E.Volume();
  int iteration = 0;
  VectorXd closest(nq);
  RandomGenerator generator(options.random_seed);
  std::vector<std::pair<double, int>> scaling(nc);
  MatrixXd closest_points(nq, nc);

  auto solver = solvers::MakeFirstAvailableSolver(
      {solvers::SnoptSolver::id(), solvers::IpoptSolver::id()});

  VectorXd guess = sample;

  // For debugging visualization.
  Vector3d point_to_draw = Vector3d::Zero();
  int num_points_drawn = 0;
  bool do_debugging_visualization = options.meshcat && nq <= 3;

  const std::string seed_point_error_msg =
      "IrisInConfigurationSpace: require_sample_point_is_contained is true but "
      "the seed point exited the initial region. Does the provided "
      "options.starting_ellipse not contain the seed point?";
  const std::string seed_point_msg =
      "IrisInConfigurationSpace: terminating iterations because the seed point "
      "is no longer in the region.";
  const std::string termination_error_msg =
      "IrisInConfigurationSpace: the termination function returned false on "
      "the computation of the initial region. Are the provided "
      "options.starting_ellipse and options.termination_func compatible?";
  const std::string termination_msg =
      "IrisInConfigurationSpace: terminating iterations because "
      "options.termination_func returned false.";

  while (true) {
    log()->info("IrisInConfigurationSpace iteration {}", iteration);
    int num_constraints = num_initial_constraints;
    HPolyhedron P_candidate = HPolyhedron(A.topRows(num_initial_constraints),
                                          b.head(num_initial_constraints));
    DRAKE_ASSERT(best_volume > 0);
    // Find separating hyperplanes

    // Add constraints from configuration space obstacles to reduce the domain
    // for later optimization.
    if (options.configuration_obstacles.size() > 0) {
      const ConvexSets& obstacles = options.configuration_obstacles;
      for (int i = 0; i < nc; ++i) {
        const auto touch = E.MinimumUniformScalingToTouch(*obstacles[i]);
        scaling[i].first = touch.first;
        scaling[i].second = i;
        closest_points.col(i) = touch.second;
      }
      std::sort(scaling.begin(), scaling.end());

      for (int i = 0; i < nc; ++i) {
        // Only add a constraint if this obstacle still has overlap with the set
        // that has been constructed so far on this iteration.
        if (HPolyhedron(A.topRows(num_constraints), b.head(num_constraints))
                .IntersectsWith(*obstacles[scaling[i].second])) {
          const VectorXd point = closest_points.col(scaling[i].second);
          AddTangentToPolytope(E, point, 0.0, &A, &b, &num_constraints);
          if (options.require_sample_point_is_contained) {
            const bool seed_point_requirement =
                A.row(num_constraints - 1) * sample <= b(num_constraints - 1);
            if (!seed_point_requirement) {
              if (iteration == 0) {
                throw std::runtime_error(seed_point_error_msg);
              }
              log()->info(seed_point_msg);
              return P;
            }
          }
          if (CheckTerminate(options,
                             HPolyhedron(A.topRows(num_constraints),
                                         b.head(num_constraints)),
                             termination_error_msg, termination_msg,
                             iteration == 0)) {
            return P;
          }
        }
      }

      P_candidate =
          HPolyhedron(A.topRows(num_constraints), b.head(num_constraints));
      MakeGuessFeasible(P_candidate, &guess);
    }

    // Use the fast nonlinear optimizer until it fails
    // num_collision_infeasible_samples consecutive times.
    for (const auto& pair_w_distance : sorted_pairs) {
      std::pair<GeometryId, GeometryId> geom_pair(pair_w_distance.geomA,
                                                  pair_w_distance.geomB);
      int consecutive_failures = 0;
      ClosestCollisionProgram prog(
          same_point_constraint, *frames.at(pair_w_distance.geomA),
          *frames.at(pair_w_distance.geomB), *sets.at(pair_w_distance.geomA),
          *sets.at(pair_w_distance.geomB), E, A.topRows(num_constraints),
          b.head(num_constraints));
      std::vector<VectorXd> prev_counter_examples =
          std::move(counter_examples[geom_pair]);
      // Sort by the current ellipsoid metric.
      std::sort(prev_counter_examples.begin(), prev_counter_examples.end(),
                [&E](const VectorXd& x, const VectorXd& y) {
                  return (E.A() * x - E.center()).squaredNorm() <
                         (E.A() * y - E.center()).squaredNorm();
                });
      std::vector<VectorXd> new_counter_examples;
      int counter_example_searches_for_this_pair = 0;
      bool warned_many_searches = false;
      while (consecutive_failures < options.num_collision_infeasible_samples) {
        // First use previous counter-examples for this pair as the seeds.
        if (counter_example_searches_for_this_pair <
            ssize(prev_counter_examples)) {
          guess = prev_counter_examples[counter_example_searches_for_this_pair];
        } else {
          MakeGuessFeasible(P_candidate, &guess);
          guess = P_candidate.UniformSample(&generator, guess,
                                            options.mixing_steps);
        }
        ++counter_example_searches_for_this_pair;
        if (do_debugging_visualization) {
          ++num_points_drawn;
          point_to_draw.head(nq) = guess;
          std::string path = fmt::format("iteration{:02}/{:03}/guess",
                                         iteration, num_points_drawn);
          options.meshcat->SetObject(path, Sphere(0.01),
                                     geometry::Rgba(0.1, 0.1, 0.1, 1.0));
          options.meshcat->SetTransform(path,
                                        RigidTransform<double>(point_to_draw));
        }
        if (prog.Solve(*solver, guess, options.solver_options, &closest)) {
          if (do_debugging_visualization) {
            point_to_draw.head(nq) = closest;
            std::string path = fmt::format("iteration{:02}/{:03}/found",
                                           iteration, num_points_drawn);
            options.meshcat->SetObject(path, Sphere(0.01),
                                       geometry::Rgba(0.8, 0.1, 0.8, 1.0));
            options.meshcat->SetTransform(
                path, RigidTransform<double>(point_to_draw));
          }
          consecutive_failures = 0;
          new_counter_examples.emplace_back(closest);
          AddTangentToPolytope(E, closest, options.configuration_space_margin,
                               &A, &b, &num_constraints);
          P_candidate =
              HPolyhedron(A.topRows(num_constraints), b.head(num_constraints));
          MakeGuessFeasible(P_candidate, &guess);
          if (options.require_sample_point_is_contained) {
            const bool seed_point_requirement =
                A.row(num_constraints - 1) * sample <= b(num_constraints - 1);
            if (!seed_point_requirement) {
              if (iteration == 0) {
                throw std::runtime_error(seed_point_error_msg);
              }
              log()->info(seed_point_msg);
              return P;
            }
          }
          if (CheckTerminate(options, P_candidate, termination_error_msg,
                             termination_msg, iteration == 0)) {
            return P;
          }
          prog.UpdatePolytope(A.topRows(num_constraints),
                              b.head(num_constraints));
        } else {
          if (do_debugging_visualization) {
            point_to_draw.head(nq) = closest;
            std::string path = fmt::format("iteration{:02}/{:03}/closest",
                                           iteration, num_points_drawn);
            options.meshcat->SetObject(path, Sphere(0.01),
                                       geometry::Rgba(0.1, 0.8, 0.8, 1.0));
            options.meshcat->SetTransform(
                path, RigidTransform<double>(point_to_draw));
          }
          if (counter_example_searches_for_this_pair >
              ssize(counter_examples[geom_pair])) {
            // Only count the failures once we start the random guesses.
            ++consecutive_failures;
          }
        }
        if (!warned_many_searches &&
            counter_example_searches_for_this_pair -
                    ssize(counter_examples[geom_pair]) >=
                10 * options.num_collision_infeasible_samples) {
          warned_many_searches = true;
          log()->info(
              " Checking {} against {} has already required {} counter-example "
              "searches; still searching...",
              inspector.GetName(pair_w_distance.geomA),
              inspector.GetName(pair_w_distance.geomB),
              counter_example_searches_for_this_pair);
        }
      }
      counter_examples[geom_pair] = std::move(new_counter_examples);
      if (warned_many_searches) {
        log()->info(
            " Finished checking {} against {} after {} counter-example "
            "searches.",
            inspector.GetName(pair_w_distance.geomA),
            inspector.GetName(pair_w_distance.geomB),
            counter_example_searches_for_this_pair);
      }
    }

    if (options.prog_with_additional_constraints) {
      counter_example_prog->UpdatePolytope(A.topRows(num_constraints),
                                           b.head(num_constraints));
      for (const auto& binding : additional_constraint_bindings) {
        for (int index = 0; index < binding.evaluator()->num_constraints();
             ++index) {
          for (bool falsify_lower_bound : {true, false}) {
            int consecutive_failures = 0;
            if (falsify_lower_bound &&
                std::isinf(binding.evaluator()->lower_bound()[index])) {
              continue;
            }
            if (!falsify_lower_bound &&
                std::isinf(binding.evaluator()->upper_bound()[index])) {
              continue;
            }
            counter_example_constraint->set(&binding, index,
                                            falsify_lower_bound);
            while (consecutive_failures <
                   options.num_additional_constraint_infeasible_samples) {
              if (counter_example_prog->Solve(
                      *solver, guess, options.solver_options, &closest)) {
                consecutive_failures = 0;
                AddTangentToPolytope(E, closest,
                                     options.configuration_space_margin, &A, &b,
                                     &num_constraints);
                P_candidate = HPolyhedron(A.topRows(num_constraints),
                                          b.head(num_constraints));
                MakeGuessFeasible(P_candidate, &guess);
                if (options.require_sample_point_is_contained) {
                  const bool seed_point_requirement =
                      A.row(num_constraints - 1) * sample <=
                      b(num_constraints - 1);
                  if (!seed_point_requirement) {
                    if (iteration == 0) {
                      throw std::runtime_error(seed_point_error_msg);
                    }
                    log()->info(seed_point_msg);
                    return P;
                  }
                }
                if (CheckTerminate(options, P_candidate, termination_error_msg,
                                   termination_msg, iteration == 0)) {
                  return P;
                }
                counter_example_prog->UpdatePolytope(A.topRows(num_constraints),
                                                     b.head(num_constraints));
              } else {
                ++consecutive_failures;
              }
              guess = P_candidate.UniformSample(&generator, guess,
                                                options.mixing_steps);
            }
          }
        }
      }
    }

    P = HPolyhedron(A.topRows(num_constraints), b.head(num_constraints));

    iteration++;
    if (iteration >= options.iteration_limit) {
      log()->info(
          "IrisInConfigurationSpace: Terminating because the iteration limit "
          "{} has been reached.",
          options.iteration_limit);
      break;
    }

    E = P.MaximumVolumeInscribedEllipsoid();
    const double volume = E.Volume();
    const double delta_volume = volume - best_volume;
    if (delta_volume <= options.termination_threshold) {
      log()->info(
          "IrisInConfigurationSpace: Terminating because the hyperellipsoid "
          "volume change {} is below the threshold {}.",
          delta_volume, options.termination_threshold);
      break;
    } else if (delta_volume / best_volume <=
               options.relative_termination_threshold) {
      log()->info(
          "IrisInConfigurationSpace: Terminating because the hyperellipsoid "
          "relative volume change {} is below the threshold {}.",
          delta_volume / best_volume, options.relative_termination_threshold);
      break;
    }
    best_volume = volume;
  }
  return P;
}
}  // namespace

HPolyhedron IrisInConfigurationSpace(const MultibodyPlant<double>& plant,
                                     const Context<double>& context,
                                     const IrisOptions& options) {
  plant.ValidateContext(context);
  const VectorXd sample = plant.GetPositions(context);
  auto same_point_constraint =
      std::make_shared<SamePointConstraint>(&plant, context);
  return _DoIris_(plant, context, options, sample, same_point_constraint);
}

HPolyhedron IrisInRationalConfigurationSpace(
    const multibody::MultibodyPlant<double>& plant,
    const systems::Context<double>& context,
    const Eigen::Ref<const Eigen::VectorXd>& q_star,
    const IrisOptions& options) {
  plant.ValidateContext(context);
  multibody::RationalForwardKinematics rational_forward_kinematics{&plant};

  const VectorXd q_sample = plant.GetPositions(context);
  const VectorXd sample =
      rational_forward_kinematics.ComputeSValue(q_sample, q_star);
  auto same_point_constraint = std::make_shared<SamePointConstraintRational>(
      &rational_forward_kinematics, q_star, context);
  return _DoIris_(plant, context, options, sample, same_point_constraint);
}

void SetEdgeContainmentTerminationCondition(
    IrisOptions* options, const Eigen::Ref<const Eigen::VectorXd>& x_1,
    const Eigen::Ref<const Eigen::VectorXd>& x_2, const double epsilon,
    const double tol) {
  const auto ab = AffineBall::MakeAffineBallFromLineSegment(x_1, x_2, epsilon);
  const auto hyperellipsoid = Hyperellipsoid(ab);
  options->starting_ellipse = hyperellipsoid;
  options->termination_func = [=](const HPolyhedron& polytope) {
    return !polytope.PointInSet(x_1, tol) || !polytope.PointInSet(x_2, tol);
  };
}

}  // namespace optimization
}  // namespace geometry
}  // namespace drake
