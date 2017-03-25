
#include "sys.h"
#include "FEStorage.h"
#include "FESolver.h"
#include "VtkProcessor.h"
#include "ReactionProcessor.h"
#include "materials/MaterialFactory.h"
#include "FEReaders.h"
#include "elements/TETRA0.h"

using namespace nla3d;

int main(){
	
	//чтение сетки
	FEStorage storage;
	
	MeshData md;
	std::string modelFilename = "../test/TETRA/tetra2_coarse.dat";
	if (!readCdbFile(modelFilename, md)) {
		LOG(FATAL) << "Can't read FE info from " << modelFilename << "file. exiting..";
	}
	md.compressNumbers();
	LinearFESolver solver;
	// add nodes
	auto sind = storage.createNodes(md.nodesNumbers.size());
	auto ind = md.nodesNumbers;
	for (uint32 i = 0; i < sind.size(); i++) {
		storage.getNode(sind[i]).pos = md.nodesPos[i];
	}
	// add elements
	for (int isolid = 1; isolid < 4;isolid++){
		ind = md.getCellsByAttribute("TYPE", isolid);
		sind = storage.createElements(ind.size(), ElementType::TETRA0);
		for (uint32 i = 0; i < sind.size(); i++) {
			ElementTETRA0& el = (ElementTETRA0&)storage.getElement(sind[i]);
			if (isolid == 1)
				el.E = 2.e11;
			else
				el.E = 2.e11;
			el.my = 0.3;
			el.getNodeNumber(0) = md.cellNodes[ind[i]][0];
			el.getNodeNumber(1) = md.cellNodes[ind[i]][1];
			el.getNodeNumber(2) = md.cellNodes[ind[i]][2];
			el.getNodeNumber(3) = md.cellNodes[ind[i]][4];
		}
	}

	for (auto v : md.feComps["X0"].list) {
		solver.addFix(v, Dof::UX);
		solver.addFix(v, Dof::UY);
		solver.addFix(v, Dof::UZ);
	}
	for (auto v : md.feComps["X1"].list) {
		solver.addFix(v, Dof::UY, 0.001);
	}
	
	solver.attachFEStorage (&storage);

	math::PARDISO_equationSolver eqSolver = math::PARDISO_equationSolver();
    solver.attachEquationSolver(&eqSolver);

	VtkProcessor* vtk = new VtkProcessor (&storage, "tetra");
    solver.addPostProcessor(vtk);
    vtk->writeAllResults();

    solver.solve();
}