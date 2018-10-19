#include "elements/TETRA0.h"

using namespace std;

namespace nla3d {

ElementTETRA0::ElementTETRA0 () {
  type = ElementType::TETRA0;
  rotmat.resize(3,3);
  rotmat << 1., 0., 0.,
            0., 1., 0.,
            0., 0., 1.;
}

void ElementTETRA0::pre () {
  for (uint16 i = 0; i < getNNodes(); i++) {
    storage->addNodeDof(getNodeNumber(i), {Dof::UX, Dof::UY, Dof::UZ});
  }
}

// here stiffness matrix is built
void ElementTETRA0::buildK() {
  Eigen::MatrixXd matS(4,4);
  matS.setZero();
  matS<< 1. , storage->getNode(getNodeNumber(0)).pos[0] , storage->getNode(getNodeNumber(0)).pos[1] , storage->getNode(getNodeNumber(0)).pos[2] ,
          1. , storage->getNode(getNodeNumber(1)).pos[0] , storage->getNode(getNodeNumber(1)).pos[1] , storage->getNode(getNodeNumber(1)).pos[2] ,
          1. , storage->getNode(getNodeNumber(2)).pos[0] , storage->getNode(getNodeNumber(2)).pos[1] , storage->getNode(getNodeNumber(2)).pos[2] ,
          1. , storage->getNode(getNodeNumber(3)).pos[0] , storage->getNode(getNodeNumber(3)).pos[1] , storage->getNode(getNodeNumber(3)).pos[2];

  vol = matS.determinant()/6.;
  // Ke will store element stiffness matrix in global coordinates
  Eigen::MatrixXd matKe(12,12);
  matKe.setZero();

  // matB is strain matrix
  Eigen::MatrixXd matB(6,12);

  // matC is 3d elastic  matrix
  Eigen::MatrixXd matC(6,6);

  // fill here matC
  makeC(matC);
  // fill here matB
  makeB(matB);

  matKe=vol*matB.transpose()*matC*matB;

  if ((alpha != 0. && T != 0.) || strains.qlength() != 0. || stress.qlength() != 0.){
    //node forces calculations
    Eigen::VectorXd Fe(12);
    Fe.setZero();

    //from initial strains
    Eigen::VectorXd strainsE(6);
    strainsE = Eigen::Map<Eigen::VectorXd>(strains.ptr(),6);

    //mechanical initial stress (replace initial from strains)
    if (stress.qlength() != 0.){
      strainsE.setZero();

      Eigen::VectorXd stressE(6);
      stressE = Eigen::Map<Eigen::VectorXd>(stress.ptr(),6);
      Eigen::MatrixXd matP;
      matP = matC.inverse();
      strainsE = matP*stressE;
    }
    //termal initial strains
    if (alpha != 0. && T != 0.){
      Eigen::VectorXd tStrains(6);
      tStrains << alpha*T,alpha*T,alpha*T,0.,0.,0.;
      strainsE = strainsE+tStrains;
    }

    Fe = (-1)*vol*matB.transpose()*matC*strainsE;
    assembleK(matKe,Fe,{Dof::UX, Dof::UY, Dof::UZ});
  }
  else{
    assembleK(matKe, {Dof::UX, Dof::UY, Dof::UZ});
  }
}

// after solution it's handy to calculate stresses, strains and other stuff in elements.
void ElementTETRA0::update () {
  // matB is strain matrix
  Eigen::MatrixXd matB(6,12);

  // matC is 3d elastic  matrix
  Eigen::MatrixXd matC(6,6);

  // fill here matC
  makeC(matC);
  // fill here matB
  makeB(matB);

  // get nodal solutions from storage
  Eigen::VectorXd U(12);
  for (uint16 i = 0; i < getNNodes(); i++) {
    U(i*3 + 0) = storage->getNodeDofSolution(getNodeNumber(i), Dof::UX);
    U(i*3 + 1) = storage->getNodeDofSolution(getNodeNumber(i), Dof::UY);
    U(i*3 + 2) = storage->getNodeDofSolution(getNodeNumber(i), Dof::UZ);
  }
  
  // restore strains
  // initial strains
  Eigen::VectorXd strainsE0(6);
  strainsE0 = Eigen::Map<Eigen::VectorXd>(strains.ptr(),6);

  Eigen::VectorXd strainsE(6);
  strainsE.setZero();
  strainsE = (-1.)*matB*U;

  //calc term strains
  if (T > 0.){
    Eigen::VectorXd tStrains(6);
    tStrains << alpha*T, alpha*T, alpha*T, 0., 0., 0.;
    strainsE = strainsE-tStrains;
  }
  strainsE = strainsE - strainsE0;
  Eigen::Map<Eigen::VectorXd>( strains.ptr(), 6) = strainsE;

  Eigen::VectorXd stressE(6);
  stressE.setZero();
  stressE = matC*strainsE;
  Eigen::Map<Eigen::VectorXd>( stress.ptr(), 6) = stressE;
}

void ElementTETRA0::makeB(Eigen::MatrixXd &B)
{
    math::Mat<6,12> matB;
    matB.zero();
    double *B_L = matB.ptr();
    double b[4], c[4], d[4];
    //Eigen::MatrixXd mb(3,3), mc(3,3), md(3,3);
    int x=0, y = 1, z=2;

    double x12 = storage->getNode(getNodeNumber(0)).pos[x] - storage->getNode(getNodeNumber(1)).pos[x];
    double x13 = storage->getNode(getNodeNumber(0)).pos[x] - storage->getNode(getNodeNumber(2)).pos[x];
    double x14 = storage->getNode(getNodeNumber(0)).pos[x] - storage->getNode(getNodeNumber(3)).pos[x];
    double x23 = storage->getNode(getNodeNumber(1)).pos[x] - storage->getNode(getNodeNumber(2)).pos[x];
    double x24 = storage->getNode(getNodeNumber(1)).pos[x] - storage->getNode(getNodeNumber(3)).pos[x];
    double x34 = storage->getNode(getNodeNumber(2)).pos[x] - storage->getNode(getNodeNumber(3)).pos[x];

    double x21 = -1.*x12;
    double x31 = -1.*x13;
    double x32 = -1.*x23;
    double x42 = -1.*x24;
    double x43 = -1.*x34;

    double y12 = storage->getNode(getNodeNumber(0)).pos[y] - storage->getNode(getNodeNumber(1)).pos[y];
    double y13 = storage->getNode(getNodeNumber(0)).pos[y] - storage->getNode(getNodeNumber(2)).pos[y];
    double y14 = storage->getNode(getNodeNumber(0)).pos[y] - storage->getNode(getNodeNumber(3)).pos[y];
    double y23 = storage->getNode(getNodeNumber(1)).pos[y] - storage->getNode(getNodeNumber(2)).pos[y];
    double y24 = storage->getNode(getNodeNumber(1)).pos[y] - storage->getNode(getNodeNumber(3)).pos[y];
    double y34 = storage->getNode(getNodeNumber(2)).pos[y] - storage->getNode(getNodeNumber(3)).pos[y];

    double y21 = -1.*y12;
    double y31 = -1.*y13;
    double y32 = -1.*y23;
    double y42 = -1.*y24;
    double y43 = -1.*y34;

    double z12 = storage->getNode(getNodeNumber(0)).pos[z] - storage->getNode(getNodeNumber(1)).pos[z];
    double z13 = storage->getNode(getNodeNumber(0)).pos[z] - storage->getNode(getNodeNumber(2)).pos[z];
    double z14 = storage->getNode(getNodeNumber(0)).pos[z] - storage->getNode(getNodeNumber(3)).pos[z];
    double z23 = storage->getNode(getNodeNumber(1)).pos[z] - storage->getNode(getNodeNumber(2)).pos[z];
    double z24 = storage->getNode(getNodeNumber(1)).pos[z] - storage->getNode(getNodeNumber(3)).pos[z];
    double z34 = storage->getNode(getNodeNumber(2)).pos[z] - storage->getNode(getNodeNumber(3)).pos[z];

    double z21 = -1.*z12;
    double z31 = -1.*z13;
    double z32 = -1.*z23;
    double z42 = -1.*z24;
    double z43 = -1.*z34;


    b[0] = y42*z32 - y32*z42;
    b[1] = y31*z43 - y34*z13;
    b[2] = y24*z14 - y14*z24;
    b[3] = y13*z21 - y12*z31;

    c[0] = x32*z42-x42*z32;
    c[1] = x43*z31-x13*z34;
    c[2] = x14*z24-x24*z14;
    c[3] = x21*z13-x31*z12;

    d[0] = x42*y32-x32*y42;
    d[1] = x31*y43-x34*y13;
    d[2] = x24*y14-x14*y24;
    d[3] = x13*y21-x12*y31;

    const double A = -1./6./vol;
    for (int i = 0; i < 4; i++){
      B_L[0 + 3*i] = b[i]*A;
      B_L[13 + 3*i] = c[i]*A;
      B_L[26 + 3*i] = d[i]*A;
      B_L[36 + 3*i] = c[i]*A;
      B_L[37 + 3*i] = b[i]*A;
      B_L[49 + 3*i] = d[i]*A;
      B_L[50 + 3*i] = c[i]*A;
      B_L[60 + 3*i] = d[i]*A;
      B_L[62 + 3*i] = b[i]*A;
    }
    B = Eigen::Map<Eigen::MatrixXd>(matB.transpose().ptr(),6,12);
}

void ElementTETRA0::makeT (Eigen::MatrixXd &T){
  //From Lekhnitskiy, The theory of anysotropic body elasticity[in Russian]
  T.setZero();

  double l1 = rotmat(0,0);
  double l2 = rotmat(1,0);
  double l3 = rotmat(2,0);
  double m1 = rotmat(0,1);
  double m2 = rotmat(1,1);
  double m3 = rotmat(2,1);
  double n1 = rotmat(0,2);
  double n2 = rotmat(1,2);
  double n3 = rotmat(2,2);

  T <<    l1*l1, m1*m1, n1*n1, 2.*l1*m1, 2.*m1*n1, 2.*n1*l1,
          l2*l2, m2*m2, n2*n2, 2.*l2*m2, 2.*m2*n2, 2.*n2*l2,
          l3*l3, m3*m3, n3*n3, 2.*l3*m3, 2.*m3*n3, 2.*n3*l3,
          l1*l2, m1*m2, n1*n2, l1*m2+m1*l2, m1*n2+m2*n1, n1*l2 + n2*l1,
          l2*l3, m2*m3, n2*n3, l2*m3+m2*l3, m2*n3+m3*n2, n2*l3+l2*n3,
          l3*l1, m3*m1, n3*n1, l3*m1+m3*l1, m3*n1+m1*n3, n3*l1+n1*l3;
}


void ElementTETRA0::makeC (Eigen::MatrixXd &C) {
  if (anisotropy == 0){
    const double A = E/((1.+my)*(1.-2.*my));
    C << (1.-my)*A , my*A,      my*A,       0.,           0.,           0.,
         my*A,       (1.-my)*A, my*A,       0.,           0.,           0.,
         my*A,        my*A,     (1.-my)*A,  0.,           0.,           0.,
         0.,          0.,       0.,         (1./2.-my)*A, 0.,           0.,
         0.,          0.,       0.,         0,            (1./2.-my)*A, 0.,
         0.,          0.,       0.,         0,            0.,           (1./2.-my)*A;
  }
  else if (anisotropy == 1){
    Eigen::MatrixXd P(6,6);
    P <<    1./EX,       -myXY/EY, -myXZ/EZ,  0.,         0.,          0.,
            -myXY/EX,    1./EY,    -myYZ/EZ,  0.,         0.,          0.,
            -myXZ/EX,    -myYZ/EY,  1./EZ,    0.,         0.,          0.,
            0.,    0.,       0.,              1./GXY,     0.,          0.,
            0.,    0.,       0.,              0.,         1./GYZ,      0.,
            0.,    0.,       0.,              0.,         0.,          1./GXZ;
    
    Eigen::MatrixXd T(6,6);
    makeT(T);
    // Elasticity matrix in local cs
    Eigen::MatrixXd Cloc(6,6);
    Cloc = P.inverse();
    C = T*Cloc*T.transpose();
  }
  else
    LOG(FATAL) << "Wrong anisotropy type!";
}

bool ElementTETRA0::getScalar(double* scalar, scalarQuery query, uint16 gp, const double scale) {
  if (query == scalarQuery::VOL){
     *scalar += vol;
    return true;
  }
  return false;
}

bool  ElementTETRA0::getTensor(math::MatSym<3>* tensor, tensorQuery query, uint16 gp, const double scale) {
  if (query == tensorQuery::C){
      tensor->comp(0,0) += strains[0];
      tensor->comp(1,1) += strains[1];
      tensor->comp(2,2) += strains[2];
      tensor->comp(0,1) += strains[3];
      tensor->comp(1,2) += strains[4];
      tensor->comp(0,2) += strains[5];
      return true;
  }
  if (query == tensorQuery::E){
    tensor->comp(0,0) += stress[0];
    tensor->comp(1,1) += stress[1];
    tensor->comp(2,2) += stress[2];
    tensor->comp(0,1) += stress[3];
    tensor->comp(1,2) += stress[4];
    tensor->comp(0,2) += stress[5];
    return true;
  }
  if (query == tensorQuery::TSTRAIN){
    tensor->comp(0,0) += alpha*T;
    tensor->comp(1,1) += alpha*T;
    tensor->comp(2,2) += alpha*T;
    tensor->comp(0,1) += 0.;
    tensor->comp(1,2) += 0.;
    tensor->comp(0,2) += 0.;
    return true;
  }
  
  return false;
}
} //namespace nla3d
