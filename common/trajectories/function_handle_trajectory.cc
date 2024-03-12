#include "drake/common/trajectories/function_handle_trajectory.h"

#include "drake/math/jacobian.h"

#include <utility>

namespace drake {
namespace trajectories {

template <typename T>
FunctionHandleTrajectory<T>::FunctionHandleTrajectory(
    std::function<MatrixX<T>(const T&)> func, int rows, int cols,
    double start_time, double end_time)
    : func_(std::move(func)),
      rows_(rows),
      cols_(cols),
      start_time_(start_time),
      end_time_(end_time) {
  DRAKE_THROW_UNLESS(func_ != nullptr);
  DRAKE_THROW_UNLESS(rows >= 0);
  DRAKE_THROW_UNLESS(cols >= 0);
  DRAKE_THROW_UNLESS(start_time <= end_time);
  DRAKE_ASSERT(func(start_time).rows() == rows);
  DRAKE_ASSERT(func(start_time).cols() == cols);
}

template <typename T>
FunctionHandleTrajectory<T>::~FunctionHandleTrajectory() = default;

template <typename T>
std::unique_ptr<Trajectory<T>> FunctionHandleTrajectory<T>::Clone() const {
  using Self = FunctionHandleTrajectory<T>;
  return std::unique_ptr<Self>(
      new Self(func_, rows_, cols_, start_time_, end_time_));
}

template <typename T>
MatrixX<T> FunctionHandleTrajectory<T>::value(const T& t) const {
  return func_(t);
}

template <typename T>
MatrixX<T> FunctionHandleTrajectory<T>::DoEvalDerivative(
    const T& t, int derivative_order) const {
  const double eps = 1e-3;
  // https://en.wikipedia.org/wiki/Five-point_stencil
  // T t0 = std::max(start_time(), t - eps);
  // T t1 = std::min(end_time(), t + eps);
  // T dt = t1 - t0;
  // DRAKE_THROW_UNLESS(dt > 0);
  if (t == start_time() || t == end_time()) {
    return Eigen::MatrixXd::Zero(rows_, cols_);
  }
  if (derivative_order == 1) {
    // return math::jacobian(func_, Eigen::Vector<double, 1>(static_cast<double>(t)));
    // return math::jacobian(func_, Eigen::Vector<T, 1>(t));
    // return math::jacobian(func_, t);
    return (value(t + eps) - value(t - eps)) / (2 * eps);
    // return (value(t1) - value(t0)) / (dt);
  //   return (-value(t + 2 * eps) + 8 * value(t + eps) - 8 * value(t - eps) + value(t - 2 * eps)) / (12 * eps);
  // } else if (derivative_order == 2) {
  //   return (-value(t + 2 * eps) + 16 * value(t + eps) - 30 * value(t) + 16 * value(t - eps) - value(t - 2 * eps)) / (12 * eps * eps);
  // } else if (derivative_order == 2) {
  //   return math::hessian(func_, Eigen::Vector<T, 1>(t));
  // } else if (derivative_order == 2) {
  //   return (value(t - eps) - (2 * value(t)) + value(t + eps)) / (eps * eps);
  } else {
    return (DoEvalDerivative(t + eps, derivative_order - 1) - DoEvalDerivative(t - eps, derivative_order - 1)) / (2 * eps);
    // return (DoEvalDerivative(t1, derivative_order - 1) - DoEvalDerivative(t0, derivative_order - 1)) / (dt);
  }
}

}  // namespace trajectories
}  // namespace drake

DRAKE_DEFINE_CLASS_TEMPLATE_INSTANTIATIONS_ON_DEFAULT_SCALARS(
    class drake::trajectories::FunctionHandleTrajectory)
