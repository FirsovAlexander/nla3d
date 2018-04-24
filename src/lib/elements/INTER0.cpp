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

  Eigen::MatrixXd T(3,6);
  
  Ke.setZero();
  K.setZero();
  T.setZero();
  
  T(0,0) = x[0]/x.length();
  T(1,1) = x[1]/x.length();
  T(2,2) = x[2]/x.length();
  T(0,3) = T(0,0);
  T(1,4) = T(1,1);
  T(2,5) = T(2,2);
  
  /*
  T(0,0) = x[0]/x.length();
  T(0,1) = x[1]/x.length();
  T(0,2) = x[2]/x.length();
  T(1,3) = T(0,0);
  T(1,4) = T(0,1);
  T(1,5) = T(0,2); 
  
  K << 1.0, -1.0,
      -1.0,  1.0;
  K *= k;
  Ke.triangularView<Eigen::Upper>() = T.transpose() * K * T;
  */
  
  
  K << x[0]/x.length()*k , 0., 0.,
       0., x[1]/x.length()*k , 0.,
       0., 0., x[2]/x.length()*k ;

  Ke <<  K , -K,
        -K,   K;

  Ke.triangularView<Eigen::Upper>() = Ke;

  assembleK(Ke, {Dof::UX, Dof::UY, Dof::UZ});
}

void ElementINTER0::update () {
  Eigen::VectorXd U(6);
  for (uint16 i = 0; i < getNNodes(); i++) {
    U(i*3 + 0) = storage->getNodeDofSolution(getNodeNumber(i), Dof::UX);
    U(i*3 + 1) = storage->getNodeDofSolution(getNodeNumber(i), Dof::UY);
    U(i*3 + 2) = storage->getNodeDofSolution(getNodeNumber(i), Dof::UZ);
  }
  deltaPos[0] = (U(0)+U(3));
  deltaPos[1] = (U(1)+U(4));
  deltaPos[2] = (U(2)+U(5));
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
