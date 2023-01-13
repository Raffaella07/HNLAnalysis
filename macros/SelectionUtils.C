#include <Math/VectorUtil.h>
#include <ROOT/RVec.hxx>
#include "Math/Vector4D.h"

using namespace ROOT::VecOps;
const RVec<float> den_sel(const RVec<float> &sel_vec,const RVec<float> &pi_pt,const RVec<float> &pi_eta,const RVec<float> &trgmu_pt,const RVec<float> &trgmu_eta,const RVec<float> &selLep_pt,const RVec<float> &selLep_eta ){
		RVec<float> den;
		std::cout << sel_vec[0] <<std::endl;
		den = sel_vec[ (pi_pt > 0.7 && abs(pi_eta) < 2.5 && trgmu_pt >7  && abs(trgmu_eta)<2 && selLep_pt>1.5 && abs(selLep_eta)<2.5) ];	
		std::cout << den[0] <<std::endl;
		return den;	
}
