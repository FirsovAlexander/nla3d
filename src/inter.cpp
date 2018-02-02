#include "sys.h"
#include "FEStorage.h"
#include "VtkProcessor.h"
#include "FESolver.h"
#include "elements/TETRA0.h"
#include "elements/INTER0.h"


using namespace nla3d;

int main (int argc, char* argv[]) {
  const uint32 numberOfNodes = 8;
  double nodeTable[numberOfNodes][3] = {
                          {0.0, 0.0, 0.0},
                          {1.,0., 0.},
                          {0., -1., 0.},
                          {0., 0., 1.},

                          {0.0, 0.0, 0.0},
                          {1.,0., 0.},
                          {0., -1., 0.},
                          {0., 0., -1.}};
  FEStorage storage;
  LinearFESolver solver;
  
  // Create and add nodes into FEStorage
  for (uint32 i = 1; i <= numberOfNodes; i++) {
    Node* no = new Node;
    no->pos[0] = nodeTable[i-1][0];
    no->pos[1] = nodeTable[i-1][1];
    no->pos[2] = nodeTable[i-1][2];
    storage.addNode(no);
  }

  ElementTETRA0* el1 = new ElementTETRA0();
  el1->E = 1.e9;
  el1->my = 0.3;
  el1->getNodeNumber(0) = 3;
  el1->getNodeNumber(1) = 2;
  el1->getNodeNumber(2) = 1;
  el1->getNodeNumber(3) = 4;
  storage.addElement(el1);

  ElementTETRA0* el2 = new ElementTETRA0();
  el2->E = 1.e9;
  el2->my = 0.3;
  el2->getNodeNumber(0) = 7;
  el2->getNodeNumber(1) = 5;
  el2->getNodeNumber(2) = 6;
  el2->getNodeNumber(3) = 8;
  storage.addElement(el2);

  ElementINTER0* inter1 = new ElementINTER0();
  inter1->k = 1e8;
  inter1->getNodeNumber(0) = 1;
  inter1->getNodeNumber(1) = 5;
  storage.addElement(inter1);

  ElementINTER0* inter2 = new ElementINTER0();
  inter2->k = 1e8;
  inter2->getNodeNumber(0) = 2;
  inter2->getNodeNumber(1) = 6;
  storage.addElement(inter2);

  ElementINTER0* inter3 = new ElementINTER0();
  inter3->k = 1e8;
  inter3->getNodeNumber(0) = 3;
  inter3->getNodeNumber(1) = 7;
  storage.addElement(inter3);

  solver.addFix(4, Dof::UX);
  solver.addFix(4, Dof::UY);
  solver.addFix(4, Dof::UZ);

  solver.addFix(8, Dof::UX);
  solver.addFix(8, Dof::UY);
  solver.addFix(8, Dof::UZ, -0.001);

  solver.attachFEStorage(&storage);

  VtkProcessor* vtk = new VtkProcessor(&storage, "inter");
  vtk->writeAllResults();
  solver.addPostProcessor(vtk);

  solver.solve();
}