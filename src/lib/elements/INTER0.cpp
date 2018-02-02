#include "elements/INTER0.h"

namespace nla3d {

ElementINTER0::ElementINTER0 () {
  type = ElementType::INTER0;
}

void ElementINTER0::pre () {
  for (uint16 i = 0; i < getNNodes(); i++) {
    storage->addNodeDof(getNodeNumber(i), {Dof::UX, Dof::UY, Dof::UZ});
  }
}

void ElementINTER0::buildK() {
  Eigen::MatrixXd Ke(6,6);
  Eigen::MatrixXd K(3,3);

  K << k , 0., 0.,
       0., k , 0.,
       0., 0., k ;

  Ke <<  K , -K,
        -K,   K;

  assembleK(Ke, {Dof::UX, Dof::UY, Dof::UZ});
}

void ElementINTER0::update () {
}

bool  ElementINTER0::getVector(math::Vec<3>* vector, vectorQuery query, uint16 gp, const double scale) {
  /*
  switch (query) {
    case vectorQuery::FLUX:
      *vector += flux*scale;
      return true;
    case vectorQuery::GRADT:
      *vector += flux*(scale/k);
      return true;
  } 
  */
  return false;
}

} //namespace nla3d
