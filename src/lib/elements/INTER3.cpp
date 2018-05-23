#include "elements/INTER3.h"

namespace nla3d {

ElementINTER3::ElementINTER3 () {
  shape = ElementShape::WEDGE;
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
     for (uint16 npj=0; npj < nOfIntPoints(); npj++) {
      //dWt = intWeight(np);
      dWt = radoIntWeight(np,npj);
      math::Mat<3,18> matB = make_B(np,npj);
      matBTDBprod(matB, D, dWt, Ke);
    }
  }// loop over integration points

  math::Mat<18,18> T = make_T();
  math::MatSym<18> KeT; 
  KeT.zero();

  matBTDBprod(T, Ke, 1., KeT);

  LOG(DEBUG) << "Ke = " << KeT.toMat(); 

  assembleK(KeT, {Dof::UX, Dof::UY, Dof::UZ});
}

void ElementINTER3::update () {
  math::Vec<18> U;
  U.zero();
  for (uint16 i = 0; i < getNNodes(); i++) {
    U[i*3 + 0] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UX);
    U[i*3 + 1] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UY);
    U[i*3 + 2] = storage->getNodeDofSolution(getNodeNumber(i), Dof::UZ);
  }
  math::Mat<18,18> T = make_T();
  strains.zero();
  for (uint16 np=0; np < nOfIntPoints(); np++) {
    for (uint16 npj=0; npj < nOfIntPoints(); npj++) {
      //double dWt = intWeight(np);
      double dWt = radoIntWeight(np,npj);
      math::Mat<3,18> matB = make_B(np,npj);
      matB = matB*T;
      matBVprod(matB, U, 1., strains);
    }
  }

  math::MatSym<3> D;
  D.zero();
  make_D(D);
  stress = D.toMat()*strains;

  //strains = make_T()*strains;
}

void ElementINTER3::makeJacob(){
  //Находим координаты узлов в локальной декартовой системе координат треугольника (t1, t2)
  //Начало координт - первый узел треугольника
  math::Vec<2> locX1(0.,0.);
  //Вектор s1 направлен по одной из сторон
  math::Vec<3> t1 = storage->getNode(getNodeNumber(1)).pos-storage->getNode(getNodeNumber(0)).pos;
  math::Vec<3> t2 = storage->getNode(getNodeNumber(2)).pos-storage->getNode(getNodeNumber(0)).pos;

  math::Vec<2> locX2(t1.length(),0.);
  //Для координат третьего узла нужна ортогонолизация. Ищем угол треугольника при начале координат
  double mult = t1[0]*t2[0] + t1[1]*t2[1] + t1[2]*t2[2];
  double angle = acos(mult/t1.length()/t2.length());
  math::Vec<2> locX3(t2.length()*cos(angle),t2.length()*sin(angle));

  math::Mat<2,2> J(locX1[0]-locX3[0],locX1[1]-locX3[1],
                   locX2[0]-locX3[0],locX2[1]-locX3[1]);
  math::Mat<2,2> invJ = J.inv(J.det());

  det = abs(invJ.det()); //Якобиан перехода между L координатами и локальными декартовыми

  LOG(INFO) << "det = " << det;
}

void ElementINTER3::make_D(math::MatSym<3>& D){
  D.comp(0,0) = ks;
  D.comp(1,1) = ks;
  D.comp(2,2) = kn;
} 

math::Mat<3,18> ElementINTER3::make_B(uint16 np, uint16 npj){
  math::Mat<3,18> B;
  B.zero();
  //Функции формы впервого порядка на треугольнике
  /*
  B[0][0] = -intL1(np);
  B[1][1]  = -intL1(np);
  B[2][2]  = -intL1(np);
  B[0][3]  = -intL2(np);
  B[1][4]  = -intL2(np);
  B[2][5]  = -intL2(np);
  B[0][6]  = -intL3(np);
  B[1][7]  = -intL3(np);
  B[2][8]  = -intL3(np);
  B[0][9]  = intL1(np);
  B[1][10] = intL1(np);
  B[2][11] = intL1(np);
  B[0][12] = intL2(np);
  B[1][13] = intL2(np);
  B[2][14] = intL2(np);
  B[0][15] = intL3(np);
  B[1][16] = intL3(np);
  B[2][17] = intL3(np);
  */
  B[0][0] = radoIntL1(np);
  B[1][1]  = radoIntL1(np);
  B[2][2]  = radoIntL1(np);
  B[0][3]  = radoIntL2(np,npj);
  B[1][4]  = radoIntL2(np,npj);
  B[2][5]  = radoIntL2(np,npj);
  B[0][6]  = radoIntL3(np,npj);
  B[1][7]  = radoIntL3(np,npj);
  B[2][8]  = radoIntL3(np,npj);
  B[0][9]  = -radoIntL1(np);
  B[1][10] = -radoIntL1(np);
  B[2][11] = -radoIntL1(np);
  B[0][12] = -radoIntL2(np,npj);
  B[1][13] = -radoIntL2(np,npj);
  B[2][14] = -radoIntL2(np,npj);
  B[0][15] = -radoIntL3(np,npj);
  B[1][16] = -radoIntL3(np,npj);
  B[2][17] = -radoIntL3(np,npj);
  return B;
}

math::Mat<18,18> ElementINTER3::make_T(){
  //Востанавливаем локальный базис s1,s2,n
  //s1 совпадает с одной из сторон
  math::Vec<3> s1 = storage->getNode(getNodeNumber(1)).pos - storage->getNode(getNodeNumber(0)).pos;
  math::Vec<3> t2 = storage->getNode(getNodeNumber(2)).pos - storage->getNode(getNodeNumber(0)).pos;
  //Востанавливаем нормаль как векторное произведение двух вектров в плоскости треугольника
  math::Vec<3> n(s1[1]*t2[2]-s1[2]*t2[1],s1[2]*t2[0]-s1[0]*t2[2],s1[0]*t2[1]-s1[1]*t2[0]);
  //Востанавливаем s2 как векторное произведение n x s1
  math::Vec<3> s2(n[1]*s1[2]-n[2]*s1[1],n[2]*s1[0]-n[0]*s1[2],n[0]*s1[1]-n[1]*s1[0]);

  n = n*(1./n.length());
  s1 = s1*(1./s1.length());
  s2 = s2*(1./s2.length());
  //Матрица поворота от глоб к локальной ск
  math::Mat<3,3> Tn (s1[0],s2[0],n[0],
                    s1[1],s2[1],n[1],
                    s1[2],s2[2],n[2]);
  math::Mat<3,3> invT = Tn.inv(Tn.det());
  //Tn = invT;
  math::Mat<18,18> T;
  T.zero();

  T[0][0] = Tn[0][0];
  T[0][1] = Tn[0][1];
  T[0][2] = Tn[0][2];

  T[1][0] = Tn[1][0];
  T[1][1] = Tn[1][1];
  T[1][2] = Tn[1][2];

  T[2][0] = Tn[2][0];
  T[2][1] = Tn[2][1];
  T[2][2] = Tn[2][2];
  //
  T[3][3] = Tn[0][0];
  T[3][4] = Tn[0][1];
  T[3][5] = Tn[0][2];

  T[4][3] = Tn[1][0];
  T[4][4] = Tn[1][1];
  T[4][5] = Tn[1][2];

  T[5][3] = Tn[2][0];
  T[5][4] = Tn[2][1];
  T[5][5] = Tn[2][2];
  //
  T[6][6] = Tn[0][0];
  T[6][7] = Tn[0][1];
  T[6][8] = Tn[0][2];

  T[7][6] = Tn[1][0];
  T[7][7] = Tn[1][1];
  T[7][8] = Tn[1][2];

  T[8][6] = Tn[2][0];
  T[8][7] = Tn[2][1];
  T[8][8] = Tn[2][2];
  //
  T[9][9] = Tn[0][0];
  T[9][10] = Tn[0][1];
  T[9][11] = Tn[0][2];

  T[10][9] = Tn[1][0];
  T[10][10] = Tn[1][1];
  T[10][11] = Tn[1][2];

  T[11][9] = Tn[2][0];
  T[11][10] = Tn[2][1];
  T[11][11] = Tn[2][2];
    //
  T[12][12] = Tn[0][0];
  T[12][13] = Tn[0][1];
  T[12][14] = Tn[0][2];

  T[13][12] = Tn[1][0];
  T[13][13] = Tn[1][1];
  T[13][14] = Tn[1][2];

  T[14][12] = Tn[2][0];
  T[14][13] = Tn[2][1];
  T[14][14] = Tn[2][2];
    //
  T[15][15] = Tn[0][0];
  T[15][16] = Tn[0][1];
  T[15][17] = Tn[0][2];

  T[16][15] = Tn[1][0];
  T[16][16] = Tn[1][1];
  T[16][17] = Tn[1][2];

  T[17][15] = Tn[2][0];
  T[17][16] = Tn[2][1];
  T[17][17] = Tn[2][2];

  return T;
}

bool  ElementINTER3::getVector(math::Vec<3>* vector, vectorQuery query, uint16 gp, const double scale) {
  return false;
}

bool ElementINTER3::getTensor(math::MatSym<3>* tensor, tensorQuery query, uint16 gp, const double scale){
  if (query == tensorQuery::C){
      tensor->comp(0,0) += strains[0];
      tensor->comp(1,1) += strains[1];
      tensor->comp(2,2) += strains[2];
      tensor->comp(0,1) += 0.;
      tensor->comp(1,2) += 0.;
      tensor->comp(0,2) += 0.;
      return true;
  }
  if (query == tensorQuery::E){
    tensor->comp(0,0) += stress[0];
    tensor->comp(1,1) += stress[1];
    tensor->comp(2,2) += stress[2];
    tensor->comp(0,1) += 0.;
    tensor->comp(1,2) += 0.;
    tensor->comp(0,2) += 0.;
    return true;
  }
}


} //namespace nla3d
