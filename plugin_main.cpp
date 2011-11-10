// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 12.09.2011 (m,d,y)

#include "registry/registry.h"
#include "common/ug_config.h"
#include "common/error.h"
#include "lib_disc/domain.h"
#include "lib_grid/lib_grid.h"

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
}
