// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#pragma once
#include "elements/element.h"

namespace nla3d {

class ElementINTER3 : public Element {
public:
  ElementINTER3();

  ElementINTER3& operator= (const ElementINTER3& from) {
    Element::operator= (from);
    return *this;
  }

  void pre();

  void buildK();
  void make_D(math::MatSym<3>& D);
  math::Mat<3,18> make_B(uint16 np);
  void make_T(math::Mat<3,3>& T);

  void makeJacob();

  void update();

  inline uint16 getNNodes() {
    return 6;
  }

  inline uint16 getDim(){
    return 3;
  }

  uint16 nOfIntPoints();

  double intWeight(uint16 np);
  double intL1(uint16 np);
  double intL2(uint16 np);
  double intL3(uint16 np);

  // strength
  double k = 0.0;

  uint16 i_int = 0; // index of integration scheme
  double det = 0.; // determinant of Jacob matrix

  //postproc procedures
  bool getVector(math::Vec<3>* vector, vectorQuery code, uint16 gp, const double scale);
};


// structures for convenient keeping of quadrature constants
struct TrianglePt {
  double L1;
  double L2;
  double L3;
  double W;
};

// gauss quadrature for 3D trinangle from (0) to (1) in L-coordinates
// 1st order
static const TrianglePt _triangle_o1[] = {
  {1./3., 1./3., 1./3., 1./2.}
};

// 2nd order
static const TrianglePt _triangle_o2[] = {
  {1./2.,     0.,     1./2.,  1./6.},
  {1./2.,     1./2.,  0.,     1./6.},
  {0.,        1./2.,  1./2.,  1./6.}
};

// array of number of quadrature points in integration scheme
static const uint16 _np_triangle[] = {
  sizeof(_triangle_o1) / sizeof(TrianglePt),
  sizeof(_triangle_o2) / sizeof(TrianglePt)
};

inline uint16 ElementINTER3::nOfIntPoints(){
  return _np_triangle[i_int];
}

static const TrianglePt* _table_triangle[] = {
  _triangle_o1,
  _triangle_o2
};

inline double ElementINTER3::intWeight(uint16 np){
  return _table_triangle[i_int][np].W * det;
}

inline double ElementINTER3::intL1(uint16 np){
  return _table_triangle[i_int][np].L1;
}

inline double ElementINTER3::intL2(uint16 np){
  return _table_triangle[i_int][np].L2;
}

inline double ElementINTER3::intL3(uint16 np){
  return _table_triangle[i_int][np].L3;
}

} //namespace nla3d
