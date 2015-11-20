/*!
 * \file vm2ug_rework.h
 * \brief vm2ug rework
 *
 * \author Stephan Grein <stephan.grein@gcsc.uni-frankfurt.de>
 * \date 19/06/2015
 */

/// guard
#ifndef  __H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG_REWORK__
#define  __H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG_REWORK__

/// includes
#include "kdtree/kd_tree.h"


/*! \defgroup mpm_plugin Membrane Potential Mapping plugin
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
	namespace membrane_potential_mapping {
		template <size_t dim, typename M> class Mapper {
		private:
			//// non-static members
			kd_tree<dim, M> m_kdtree;
			size_t m_delimIndex;

			/// common typedefs
			typedef typename std::vector<std::pair<MathVector<dim, number>, M> >::const_iterator CITVPMNM;
			typedef typename std::vector<std::pair<std::vector<number>, M> >::const_iterator CITVPVNM;

		public:
			/*!
			 * \brief ctor
			 */
			Mapper() : m_kdtree(kd_tree<dim, M>()), m_delimIndex(0) {

			}

			/*!
			 * \brief build an empty tree
			 */
			void build_tree();

			/*!
			 * \brief build populated tree from given vector of pairs
			 * \param[in] points to be used for the tree construction
			 */
			void build_tree(const std::vector<std::pair<MathVector<dim, number>, M> >& points);

			/*!
			 * \brief build populated tree from given file
			 * \param[in] filename where the points are stored (with meta data)
			 * \param[in] delimiter in the file to separate values
			 */
			void build_tree(const std::string& filename, const std::string& delim);

			/*!
			 * \brief build populated tree from given vector of pairs
			 * \param[in] points to be used for the tree construction
			 */
			void build_tree(const std::vector<std::pair<std::vector<number>, M> >& points);

			/*!
			 * \brief add a single node
			 * \param[in] node a single point with meta data
			 */
			void add_node(const std::pair<MathVector<dim, number>, M>& node);

			/*!
			 * \brief add a single node
			 * \param[in] node a single point with meta data
			 */
			void add_node(const std::pair<std::vector<number>, M>& node);

			/*!
			 * \brief add a single node with a value
			 * \param[in] node a single point with meta data
			 * \param[in] value the meta data
			 */
			inline void add_node(const std::vector<number>&, const M& value);

			/*!
			 * \brief query the tree for the data of the very nearest neighbor
			 * \param[in] query coordinates of a given point
			 */
			M get_data_from_nearest_neighbor(const MathVector<dim, number>& query) const;

			/*!
			 * \brief query the tree for the data of the very nearest neighbor
			 * \param[in] query coordinates of a given point
			 */
			M get_data_from_nearest_neighbor(const std::vector<number>& query) const;

		};
	} /// end namespace mpm
} /// end namespace ug
//<! \}

#include "vm2ug_rework_impl.h"

#endif /// __H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG_REWORK__
