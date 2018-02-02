// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "elements/element.h"

namespace nla3d {

class ElementINTER0 : public ElementVERTEX {
public:
  ElementINTER0 ();

  void pre();

  void buildK();

  void update();

  // strength
  double k = 0.0;

  //postproc procedures
  bool getVector(math::Vec<3>* vector, vectorQuery code, uint16 gp, const double scale);
};

} //namespace nla3d
