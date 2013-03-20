/*!
 *
 * \file circumference.h
 * \brief auxiliary functions to calculate circumference
 *
 * \author Stephan Grein
 * \date Tue, 19 Mar 2013
 *
 */

#ifndef __H__UG_MEMBRANE_POTENTIAL_MAPPING__AUX__EDGE_UTILITIES__
#define __H__UG_MEMBRANE_POTENTIAL_MAPPING__AUX__EDGE_UTILITIES__

#include "lib_grid/algorithms/subset_util.h"
#include "lib_grid/tools/subset_handler_interface.h"

/*! begin namespace ug */
namespace ug {
	/*! begin namespace mpm */
	namespace membrane_potential_mapping {
		/*! begin namespace aux */
		namespace aux {

			/*!
			 *
			 * \fn EdgeSum
			 * \brief iterator for calculation of the linear edge sum for an edge subset
			 *
			 * \param[in] eBegin begin of an edge iterator
			 * \param[in] eEnd end of an edge iterator
			 * \param[in] aaPos the attachment
			 *
			 * \return \c the linear edge sum as a number
			 *
			 */
			template <typename TEdgeIterator, typename TAAPosVRT>
			number EdgeSum(TEdgeIterator eBegin, TEdgeIterator eEnd, TAAPosVRT& aaPos);

			/*!
			 *
			 * \fn EdgeSumSq
			 * \brief iterator for calculation of the squared edge sum for an edge subset
			 *
			 * \param[in] eBegin begin of an edge iterator
			 * \param[in] eEnd end of an edge iterator
			 * \param[in] aaPos the attachment
			 *
			 * \return \c the squared edge sum as a number
			 *
			 */
			template <typename TEdgeIterator, typename TAAPosVRT>
			number EdgeSum(TEdgeIterator eBegin, TEdgeIterator eEnd, TAAPosVRT& aaPos);

			/*!
			 *
			 * \fn EdgeSum
			 * \brief calculates the linear edge sum for a given subset si on level lvl
			 *
			 * \param[in] sh subset handler
			 * \param[in] si subset index
			 * \param[in] lvl grid level
			 * \param[in] aaPos the attachment
			 *
			 * \return \c the linear edge sum as a number
			 *
			 */
			template <typename TAAPosVRT>
			number EdgeSum(ISubsetHandler& sh, int si, int lvl, TAAPosVRT& aaPos);

			/*!
			 *
			 * \fn EdgeSumSq
			 * \brief calculates the squared edge sum for a given subset si on level lvl
			 *
			 * \param[in] sh subset handler
			 * \param[in] si subset index
			 * \param[in] lvl grid level
			 * \param[in] aaPos the attachment
			 *
			 * \return \c the squared edge sum as a number
			 *
			 */
			template <typename TAAPosVRT>
			number EdgeSumSq(ISubsetHandler& sh, int si, int lvl, TAAPosVRT& aaPos);

		/*! end namespace aux */
		}
	/*! end namespace mpm */
	}
/*! end namespace ug */
}

#include "edge_utilities_impl.hpp"

#endif //__H__UG_MEMBRANE_POTENTIAL_MAPPING__AUX__EDGE_UTILITIES__
