// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 12.09.2011 (m,d,y)

#include "registry/registry.h"
#include "common/ug_config.h"
#include "common/error.h"
#include "lib_disc/domain.h"
#include "lib_grid/lib_grid.h"

#include "vm2ug.h"
#include "bg.h"

#include <string>

using namespace std;

int SayHello()
{
	UG_LOG("HELLO From MembranePotentialMapping!\n");
	return 0;
}

extern "C" UG_API void InitUGPlugin(ug::bridge::Registry* reg, string parentGroup)
{
	string grp(parentGroup); grp.append("Neuro/");

	reg->add_function("SayHelloMPM", &SayHello, grp);

   typedef vug::Vm2uG<std::string> TVm2uG;
reg->add_class_<TVm2uG>("MembranePotentialMapper", grp)
       .add_constructor<void (*)(std::string, std::string, bool)>("Timestep-File#File-Extension#NodesCanChange")
       .add_method("get_potential", &TVm2uG::get_potential, grp)
       .add_method("build_tree", &TVm2uG::buildTree, grp);
   
   typedef bg::BG TBG;
   reg->add_class_<TBG>("BorgGraham", grp)
       .add_constructor()
       .add_method("install_can_gates", &TBG::install_can_gates, grp)
       .add_method("get_current", &TBG::timestepping_of_gates_and_calc_current, grp)
       .add_method("calc_current_at_start", &TBG::calc_current_at_start, grp)
       .add_method("get_Neumann_Flux", &TBG::get_Neumann_Flux, grp);
}
