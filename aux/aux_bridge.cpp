/* mandatory includes */
#include "edge_utilities.h"
#include "lib_grid/algorithms/subset_util.h"
#include "lib_grid/tools/subset_handler_interface.h"

/* begin namespace ug */
namespace ug {
	/* begin namespace mpm */
	namespace membrane_potential_mapping {
		/* begin namespace aux */
		namespace aux {
			/* wrapper for EdgeSum (\see edge_utilities.h) */
			template <typename TDomain>
			static number EdgeSum(TDomain& dom, ISubsetHandler& sh, int si, int lvl)
			{
				typename TDomain::position_accessor_type& aaPos = dom.position_accessor();
				UG_ASSERT(TDomain::position_type::Size <= 3, "too many coordinates.");

				return EdgeSum(sh, si, lvl, aaPos);
			}

			/* wrapper for EdgeSumSq (\see edge_utilities.h) */
		    template <typename TDomain>
		    static number EdgeSumSq(TDomain& dom, ISubsetHandler& sh, int si, int lvl)
		    {
			   typename TDomain::position_accessor_type& aaPos = dom.position_accessor();
		 	   UG_ASSERT(TDomain::position_type::Size <= 3, "too many coordinates.");

    		   return EdgeSumSq(sh, si, lvl, aaPos);
		    }
		/* end namespace aux */
		}
	/* end namespace mpm */
	}
/* end namespace ug*/
}
