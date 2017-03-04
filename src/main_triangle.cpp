
#include "sys.h"
#include "FEStorage.h"
#include "FESolver.h"
#include "VtkProcessor.h"
#include "ReactionProcessor.h"
#include "materials/MaterialFactory.h"
#include "FEReaders.h"
#include "elements/TRIANGLE4.h"

using namespace nla3d;

int main(){
	//чтение сетки
	FEStorage storage;
	
	MeshData md;
	std::string modelFilename = "../test/TRIANGLE/triangle.cdb";
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
	ind = md.getCellsByAttribute("TYPE", 1);
	sind = storage.createElements(ind.size(), ElementType::TRIANGLE4);
	for (uint32 i = 0; i < sind.size(); i++) {
		ElementTRIANGLE4& el = (ElementTRIANGLE4&)storage.getElement(sind[i]);
		el.E = 2.e11;
		el.my = 0.3;
		el.state = PlaneState::Strain;
		for (uint16 j = 0; j < el.getNNodes(); j++) {
			el.getNodeNumber(j) = md.cellNodes[ind[i]][j];
		}
	}

	for (auto v : md.feComps["FIX"].list) {
		solver.addFix(v, Dof::UX);
		solver.addFix(v, Dof::UY);
	}
	for (auto v : md.feComps["DISP"].list) {
		solver.addFix(v, Dof::UX, 0.01);
	}
	
	solver.attachFEStorage (&storage);

	VtkProcessor* vtk = new VtkProcessor (&storage, "jobname");
    solver.addPostProcessor(vtk);

    solver.solve();
}