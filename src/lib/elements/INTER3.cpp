#include "elements/INTER3.h"

namespace nla3d {

ElementINTER3::ElementINTER3 () {
  shape = ElementShape::TRIANGLE;
  nodes = new uint32[getNNodes()];
  type = ElementType::INTER3;
}

void ElementINTER3::pre () {
  for (uint16 i = 0; i < getNNodes(); i++) {
    storage->addNodeDof(getNodeNumber(i), {Dof::UX, Dof::UY, Dof::UZ});
  }
}

void ElementINTER3::buildK() {
  math::MatSym<18> Ke;
  Ke.zero();
  math::MatSym<3> D;
  D.zero();

  make_D(D);

  makeJacob();

  // build Ke
  double dWt; //Gaussian quadrature weight
  for (uint16 np=0; np < nOfIntPoints(); np++) {
    dWt = intWeight(np);
    math::Mat<3,18> matB = make_B(np);
    matBTDBprod(matB, D, dWt, Ke);
  }// loop over integration points

  assembleK(Ke, {Dof::UX, Dof::UY, Dof::UZ});
}

void ElementINTER3::update () {


}

void ElementINTER3::makeJacob(){
  //Находим координаты узлов в локальной декартовой системе координат треугольника (t1, t2)
  //Начало координт - первый узел треугольника
  math::Vec<2> locX1(0,0);
  //Вектор s1 направлен по одной из сторон
  math::Vec<3> t1 = storage->getNode(getNodeNumber(1)).pos-storage->getNode(getNodeNumber(0)).pos;
  math::Vec<3> t2 = storage->getNode(getNodeNumber(2)).pos-storage->getNode(getNodeNumber(0)).pos;
  math::Vec<2> locX2(t1.length(),0);
  //Для координат третьего узла нужна ортогонолизация. Ищем угол треугольника при начале координат
  double mult = t1[0]*t2[0] + t1[1]*t2[1] + t1[2]*t2[2];
  double angle = acos(mult/t1.length()/t2.length());
  math::Vec<2> locX3(t2.length()*cos(angle),t2.length()*sin(angle));

  math::Mat<2,2> J(locX1[0]-locX3[0],locX1[1]-locX3[1],
                   locX2[0]-locX3[0],locX2[1]-locX3[1]);
  det = J.det(); //Якобиан перехода между L координатами и локальными декартовыми
}

void ElementINTER3::make_D(math::MatSym<3>& D){
  D.comp(0,0) = k;
  D.comp(1,1) = k;
  D.comp(2,2) = k;
}

math::Mat<3,18> ElementINTER3::make_B(uint16 np){
  math::Mat<3,18> B;
  B.zero();
  //Функции формы впервого порядка на треугольнике
  B[0][0] = intL1(np);
  B[1][1]  = intL1(np);
  B[2][2]  = intL1(np);
  B[0][3]  = intL2(np);
  B[1][4]  = intL2(np);
  B[2][5]  = intL2(np);
  B[0][6]  = intL3(np);
  B[1][7]  = intL3(np);
  B[2][8]  = intL3(np);
  B[0][9]  = -intL1(np);
  B[1][10] = -intL1(np);
  B[2][11] = -intL1(np);
  B[0][12] = -intL2(np);
  B[1][13] = -intL2(np);
  B[2][14] = -intL2(np);
  B[0][15] = -intL3(np);
  B[1][16] = -intL3(np);
  B[2][17] = -intL3(np);

  //матрица поворота в глобальную декартову СК
  B = make_T()*B;

  return B;
}

math::Mat<3,3> ElementINTER3::make_T(){
  //Востанавливаем локальный базис s1,s2,n
  //s1 совпадает с одной из сторон
  math::Vec<3> s1 = storage->getNode(getNodeNumber(1)).pos - storage->getNode(getNodeNumber(0)).pos;
  math::Vec<3> t2 = storage->getNode(getNodeNumber(2)).pos - storage->getNode(getNodeNumber(0)).pos;
  //Востанавливаем нормаль как векторное произведение двух вектров в плоскости треугольника
  math::Vec<3> n(s1[1]*t2[2]-s1[2]*t2[1],s1[2]*t2[0]-s1[0]*t2[2],s1[0]*t2[1]-s1[1]*t2[0]);
  //Востанавливаем s2 как векторное произведение n x s1
  math::Vec<3> s2(n[1]*s1[2]-n[2]*s1[1],n[2]*s1[0]-n[0]*s1[2],n[0]*s1[1]-n[1]*s1[0]);
  
  //Матрица поворота от глоб к локальной ск
  math::Mat<3,3> T ( s1[0],s2[0],n[0],
                    s1[1],s2[1],n[1],
                    s1[2],s2[2],n[2]);
  return T;
}

bool  ElementINTER3::getVector(math::Vec<3>* vector, vectorQuery query, uint16 gp, const double scale) {
  return false;
}

} //namespace nla3d
