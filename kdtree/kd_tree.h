/*!
 * \file kd_tree.h
 *
 *  Created on: Apr 21, 2015
 *      Author: stephan
 */

/// guard
#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__KD_TREE__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__KD_TREE__

/// includes
#include "kd_node.h"
#include "common/math/ugmath.h"

/* \defgroup mpm_plugin Membrane Potential Mapping plugin
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
	namespace membrane_potential_mapping {
		/*!
		 * \brief pretty simplistic KD Tree
		 * Note the worst case runtime for building the static kd tree,
		 * is dominated by the utilized sorting algorithm and here O(n * log(n))
		 */
		template <size_t dim, typename M>
		class kd_tree {
		private:
			int m_visited;
        	std::vector<kd_node<dim, M> > nodes;
        	kd_node<dim, M>* m_pRoot;
    		kd_node<dim, M>* m_pWps;

		public:
			/*!
			 * \brief ctor
			 */
			kd_tree();

			/*!
			 * \brief dtor
			 */
			~kd_tree();

		private:
			/*!
			 * \brief euclidean distance in dim-d
			 * \param[in] a kd_node
			 * \param[in] b kd_node
			 * \param[out] M distance
			 */
			M dist(kd_node<dim, M>* a, kd_node<dim, M> *b) const;

			/*!
			 * \brief swap nodes
			 * \param[in] x kd_node
			 * \param[in] y kd_node
			 */
			void swap(kd_node<dim, M>* x, kd_node<dim, M>* y);

			/*!
			 * \brief find median by median-of-medians, a quickselect method
			 * Note the runtime and memory complexity:
			 * 		- worst case runtime: O(n) and
			 * 		- worst case memory: O(1)
			 *
			 * \param[in] start of kd_nodes
			 * \param[in] end of kd_nodes
			 * \param[in] idx index
			 * \param[out] kd_node representing median
			 */
			kd_node<dim, M>* find_median(kd_node<dim, M> *start, kd_node<dim, M> *end, size_t idx);

			/*!
			 * \brief generates a kd tree from a list of kd nodes
			 * \param[in] t list of kd_nodes
			 * \param[in] len length of kd_nodes
			 * \param[in] i index
			 * \param[out] kd_node root
			 */
			kd_node<dim, M>* make_tree(kd_node<dim, M>* t, size_t len, size_t i);

			/*!
			 * \brief get nearest neighbor in tree wrt to a query point (euclidean distance is used)
			 * \param[in] root kd_tree root
			 * \param[in] i index
			 * \param[out] best kd_nodes's best node
			 * \param[out] best_dist associated distance to best node (see above)
			 */
			void nearest(kd_node<dim, M> *root, kd_node<dim, M> *nd, size_t i, kd_node<dim, M> **best, M *best_dist);

       public:
			/*!
			 * \brief add a node w meta data
			 * \param[in] vec node to be added
			 * \param[in] m meta data
			 */
        	void add_node_meta(const MathVector<dim>& vec, number m);

        	/*!
        	 * \brief add a node w/o meta data
			 * \param[in] vec node to be added
        	 */
        	void add_node(const MathVector<dim>& vec);

           	/*!
        	 * \brief query the tree with a given vector for the nearest neighbor and return the attached metadata
        	 */
        	bool build_tree();

        	/*!
        	 * \brief query the tree for the nearest neighbor
        	 * \param[in] vec query node
        	 * \param[out] M meta data
        	 */
        	M query(const MathVector<dim, M>& vec);
         };
	} // namespace membrane_potential_mapping
} // namespace ug
//<! \}

#include "kd_tree_impl.h"

#endif // __H__UG__MEMBRANE_POTENTIAL_MAPPING__KD_TREE__
