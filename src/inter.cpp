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
                          {0., -1., 0.},
                          {0., 0., 1.},

                          {0.0, 0.0, 0.0},
                          {1.,0., 0.},
                          {0., -1., 0.},
                          {0., 0., -1.}};
  double kn = 1e8;
  double ks = 1e8;

  math::Vec<3> loc = {0.,0.,1.};

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
  inter1->ks = ks;
  inter1->n = loc;
  inter1->getNodeNumber(0) = 1;
  inter1->getNodeNumber(1) = 5;

  ElementINTER0* inter11 = new ElementINTER0();
  inter11->kn = kn;
  inter11->ks = ks;
  inter11->n = loc;
  inter11->getNodeNumber(0) = 1;
  inter11->getNodeNumber(1) = 6;

  ElementINTER0* inter12 = new ElementINTER0();
  inter12->kn = kn;
  inter12->ks = ks;
  inter12->n = loc;
  inter12->getNodeNumber(0) = 1;
  inter12->getNodeNumber(1) = 7;

  ElementINTER0* inter2 = new ElementINTER0();
  inter2->kn = kn;
  inter2->ks = kn;
  inter2->n = loc;
  inter2->getNodeNumber(0) = 2;
  inter2->getNodeNumber(1) = 6;

  ElementINTER0* inter21 = new ElementINTER0();
  inter21->kn = kn;
  inter21->ks = kn;
  inter21->n = loc;
  inter21->getNodeNumber(0) = 2;
  inter21->getNodeNumber(1) = 5;

  ElementINTER0* inter22 = new ElementINTER0();
  inter22->kn = kn;
  inter22->ks = kn;
  inter22->n = loc;
  inter22->getNodeNumber(0) = 2;
  inter22->getNodeNumber(1) = 7;

  ElementINTER0* inter3 = new ElementINTER0();
  inter3->kn = kn;
  inter3->ks = kn;
  inter3->n = loc;
  inter3->getNodeNumber(0) = 3;
  inter3->getNodeNumber(1) = 7;

  ElementINTER0* inter31 = new ElementINTER0();
  inter31->kn = kn;
  inter31->ks = kn;
  inter31->n = loc;
  inter31->getNodeNumber(0) = 3;
  inter31->getNodeNumber(1) = 5;

  ElementINTER0* inter32 = new ElementINTER0();
  inter32->kn = kn;
  inter32->ks = kn;
  inter32->n = loc;
  inter32->getNodeNumber(0) = 3;
  inter32->getNodeNumber(1) = 6;
  
  
  
  /*1D interface solution */
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
  storage.addElement(inter11);
  storage.addElement(inter12);
  storage.addElement(inter21);
  storage.addElement(inter22);
  storage.addElement(inter31);
  storage.addElement(inter32);
  
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

  math::PARDISO_equationSolver eqSolver = math::PARDISO_equationSolver();
  solver.attachEquationSolver(&eqSolver);
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
  solver2.addFix(4, Dof::UY);
  solver2.addFix(4, Dof::UZ, 0.001);
  
  solver2.attachFEStorage(&storage2);

  VtkProcessor* vtk2 = new VtkProcessor(&storage2, "inter2");
  vtk2->writeAllResults();
  solver2.addPostProcessor(vtk2);

  math::PARDISO_equationSolver eqSolver2 = math::PARDISO_equationSolver();
  solver2.attachEquationSolver(&eqSolver2);  
  solver2.solve();
  
}