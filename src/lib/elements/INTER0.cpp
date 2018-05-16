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

  Eigen::MatrixXd T(3,3);
  
  Ke.setZero();
  K.setZero();

  //Произвольный линенйно независмый базис
  math::Vec<3> s1(n[0],n[1]+1.,n[2]+1.);
  //ортогонализаця базиса
  math::Vec<3> proj_s1_n = n*((s1[0]*n[0]+s1[1]*n[1]+s1[2]*n[2])/n.qlength());
  s1 = s1 - proj_s1_n;
  math::Vec<3> s2(n[1]*s1[2]-n[2]*s1[1],n[2]*s1[0]-n[0]*s1[2],n[0]*s1[1]-n[1]*s1[0]);
  
  n = n*(1./n.length());
  s1 = s1*(1./s1.length());
  s2 = s2*(1./s2.length());

  T << s1[0], s2[0], n[0],
       s1[1], s2[1], n[1],
       s1[2], s2[2], n[2];

  Eigen::MatrixXd invT(3,3);
  invT = T.inverse();

  K << ks , 0, 0,
       0,  ks, 0,
       0,  0, kn;

  K = invT.transpose()*K*invT;

  Ke << K ,  -K,
       -K,    K;

  LOG(INFO) << Ke;

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
