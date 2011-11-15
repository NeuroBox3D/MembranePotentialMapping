// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 12.09.2011 (m,d,y)

#include "registry/registry.h"
#include "common/ug_config.h"
#include "common/error.h"
#include "lib_disc/domain.h"
#include "lib_grid/lib_grid.h"

#include "vm2ug.hh"

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

}
