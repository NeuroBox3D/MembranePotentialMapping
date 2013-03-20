/*! mandatory includes */
#include "edge_utilities.h"
#include <lib_grid/algorithms/geom_obj_util/edge_util.h>
#include <lib_grid/grid/geometric_base_objects.h>

/*! begin namespace ug */
namespace ug {
	/*! begin namespace mpm */
	namespace membrane_potential_mapping {
		/*! begin namespace aux */
		namespace aux {

			template <typename TEdgeIterator, typename TAAPosVRT>
			number EdgeSum(TEdgeIterator eBegin, TEdgeIterator eEnd, TAAPosVRT& aaPos)
			{
				number circumference = 0.;

				if(eBegin == eEnd)
					return circumference;

				for(; eBegin != eEnd; ++eBegin) {
					EdgeBase* cEdge = *eBegin;
					circumference += EdgeLength(cEdge, aaPos);
				}
			}

			template <typename TEdgeIterator, typename TAAPosVRT>
			number EdgeSumSq(TEdgeIterator eBegin, TEdgeIterator eEnd, TAAPosVRT& aaPos)
			{
				number circumference = 0.;

				if(eBegin == eEnd)
					return circumference;

				for(; eBegin != eEnd; ++eBegin) {
					EdgeBase* cEdge = *eBegin;
					circumference += EdgeLength(cEdge, aaPos);
				}
			}


			template <typename TAAPosVRT>
			number EdgeSum(ISubsetHandler& sh, int si, int lvl, TAAPosVRT& aaPos)
			{
				number circumference = 0.;
				GeometricObjectCollection goc = sh.get_geometric_objects_in_subset(si);

				if (goc.num<Edge>(lvl) == 0) {
					UG_WARNING("WARNING: Given subset doesn't contain any edges on provided level " << lvl << std::endl);
				} else {
					for (EdgeBaseIterator eIt = goc.edges_begin(lvl); eIt != goc.edges_end(lvl); ++eIt) {
						EdgeBase* curEdge = *eIt;
						circumference += EdgeLength(curEdge, aaPos);
					}
				}

				return circumference;
			}

			template <typename TAAPosVRT>
			number EdgeSumSq(ISubsetHandler& sh, int si, int lvl, TAAPosVRT& aaPos)
			{
				number circumference = 0.;
				GeometricObjectCollection goc = sh.get_geometric_objects_in_subset(si);

				if (goc.num<Edge>(lvl) == 0) {
					UG_WARNING("WARNING: Given subset doesn't contain any edges on provided level " << lvl << std::endl);
				} else {
					for (EdgeBaseIterator eIt = goc.edges_begin(lvl); eIt != goc.edges_end(lvl); ++eIt) {
						EdgeBase* curEdge = *eIt;
						circumference += EdgeLength(curEdge, aaPos);
					}
				}

				return circumference;
			}
		/*! end namespace aux */
		}
	/*! end namespace mpm */
	}
/*! end namespace ug */
}
