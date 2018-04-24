#include "sys.h"
#include "FEStorage.h"
#include "VtkProcessor.h"
#include "FESolver.h"
#include "elements/TETRA0.h"
#include "elements/INTER0.h"
#include "elements/INTER3.h"


using namespace nla3d;

int main (int argc, char* argv[]) {
  const uint32 numberOfNodes = 8;
  double nodeTable[numberOfNodes][3] = {
                          {0.0, 0.0, 0.0},
                          {1.,0., 0.},
                          {0., -0.8660254037844387, 0.49999999999999994},
                          {0., 0.49999999999999994, 0.8660254037844387},

                          {0.0, 0.0, 0.0},
                          {1.,0., 0.},
                          {0., -0.8660254037844387, 0.49999999999999994},
                          {0.,-0.49999999999999994, -0.8660254037844387}};
  double kn = 1e10;
  double ks = 1e20;
  /*

  ElementTETRA0* el1 = new ElementTETRA0();
  el1->E = 1.e9;
  el1->my = 0.3;
  el1->getNodeNumber(0) = 3;
  el1->getNodeNumber(1) = 2;
  el1->getNodeNumber(2) = 1;
  el1->getNodeNumber(3) = 4;
  
  ElementTETRA0* el2 = new ElementTETRA0();
  el2->E = 1.e9;
  el2->my = 0.3;
  el2->getNodeNumber(0) = 7;
  el2->getNodeNumber(1) = 5;
  el2->getNodeNumber(2) = 6;
  el2->getNodeNumber(3) = 8;
  
  ElementINTER0* inter1 = new ElementINTER0();
  inter1->kn = kn;
  inter1->ks = kn;
  inter1->getNodeNumber(0) = 1;
  inter1->getNodeNumber(1) = 5;

  ElementINTER0* inter2 = new ElementINTER0();
  inter2->kn = kn;
  inter2->ks = ks;
  inter2->getNodeNumber(0) = 2;
  inter2->getNodeNumber(1) = 6;

  ElementINTER0* inter3 = new ElementINTER0();
  inter3->kn = kn;
  inter3->ks = ks;
  inter3->getNodeNumber(0) = 3;
  inter3->getNodeNumber(1) = 7;
  
  /*1D interface solution */
  /*
  FEStorage storage;
  // Create and add nodes into FEStorage
  for (uint32 i = 1; i <= numberOfNodes; i++) {
    Node* no = new Node;
    no->pos[0] = nodeTable[i-1][0];
    no->pos[1] = nodeTable[i-1][1];
    no->pos[2] = nodeTable[i-1][2];
    storage.addNode(no);
  }
  storage.addElement(el1);
  storage.addElement(el2);
  storage.addElement(inter1);
  storage.addElement(inter2);
  storage.addElement(inter3);
  
  LinearFESolver solver;
  solver.addFix(8, Dof::UX);
  solver.addFix(8, Dof::UY);
  solver.addFix(8, Dof::UZ);

  solver.addFix(4, Dof::UX);
  solver.addFix(4, Dof::UY);
  solver.addFix(4, Dof::UZ, 0.001);

  solver.attachFEStorage(&storage);

  VtkProcessor* vtk = new VtkProcessor(&storage, "inter");
  vtk->writeAllResults();
  solver.addPostProcessor(vtk);
  solver.solve();

  /*2D interface solution*/
  
  FEStorage storage2;
  // Create and add nodes into FEStorage
  for (uint32 i = 1; i <= numberOfNodes; i++) {
    Node* no = new Node;
    no->pos[0] = nodeTable[i-1][0];
    no->pos[1] = nodeTable[i-1][1];
    no->pos[2] = nodeTable[i-1][2];
    storage2.addNode(no);
  }

  ElementTETRA0* el3 = new ElementTETRA0();
  el3->E = 1.e9;
  el3->my = 0.3;
  el3->getNodeNumber(0) = 3;
  el3->getNodeNumber(1) = 2;
  el3->getNodeNumber(2) = 1;
  el3->getNodeNumber(3) = 4;
  
  ElementTETRA0* el4 = new ElementTETRA0();
  el4->E = 1.e9;
  el4->my = 0.3;
  el4->getNodeNumber(0) = 7;
  el4->getNodeNumber(1) = 5;
  el4->getNodeNumber(2) = 6;
  el4->getNodeNumber(3) = 8;

  ElementINTER3* inter4 = new ElementINTER3();
  inter4->kn = kn;
  inter4->ks = ks;
  inter4->getNodeNumber(0) = 1;
  inter4->getNodeNumber(1) = 2;
  inter4->getNodeNumber(2) = 3;
  inter4->getNodeNumber(3) = 5;
  inter4->getNodeNumber(4) = 6;
  inter4->getNodeNumber(5) = 7;

  storage2.addElement(el3);
  storage2.addElement(el4);
  storage2.addElement(inter4);
  
  LinearFESolver solver2;
  solver2.addFix(8, Dof::UX);
  solver2.addFix(8, Dof::UY);
  solver2.addFix(8, Dof::UZ);

  solver2.addFix(4, Dof::UX);
  solver2.addFix(4, Dof::UY, 0.001*0.49999999999999994);
  solver2.addFix(4, Dof::UZ, 0.001*0.8660254037844387);
  
  solver2.attachFEStorage(&storage2);

  VtkProcessor* vtk2 = new VtkProcessor(&storage2, "inter2");
  vtk2->writeAllResults();
  solver2.addPostProcessor(vtk2);
  solver2.solve();
  
}