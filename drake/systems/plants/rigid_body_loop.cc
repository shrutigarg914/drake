#include "drake/systems/plants/rigid_body_loop.h"

std::ostream& operator<<(std::ostream& os, const RigidBodyLoop& obj) {
  os << "loop connects pt "
     << obj.frameA->transform_to_body.matrix().topRightCorner(3, 1).transpose()
     << " on " << obj.frameA->body->name_ << " to pt "
     << obj.frameB->transform_to_body.matrix().topRightCorner(3, 1).transpose()
     << " on " << obj.frameB->body->name_ << std::endl;
  return os;
}
