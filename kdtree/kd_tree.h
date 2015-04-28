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

/* \defgroup mpm_plugin Membrane Potential Mapping plugin
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
	namespace membrane_potential_mapping {
		/*!
		 * \brief pretty simplistic KD Tree
		 */
		template <size_t dim, typename M>
		class kd_tree {
		private:
			int m_visited;
        	std::vector<kd_node<dim, M> > nodes;
        	kd_node<dim, M>* root;
    		kd_node<dim, M>* wps;

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
			 */
			M dist(kd_node<dim, M>* a, kd_node<dim, M> *b) const;

			/*!
			 * \brief swap nodes
			 */
			void swap(kd_node<dim, M>* x, kd_node<dim, M>* y);

			/*!
			 * \brief find median
			 */
			kd_node<dim, M>* find_median(kd_node<dim, M> *start, kd_node<dim, M> *end, int idx);

			/*!
			 * \brief generates a kd tree from a list of kd nodes
			 */
			kd_node<dim, M>* make_tree(kd_node<dim, M>* t, int len, int i);

			/*!
			 * \brief get nearest neighbor in tree wrt to a query point (euclidean distance is used)
			 */
			void nearest(kd_node<dim, M> *root, kd_node<dim, M> *nd, int i, kd_node<dim, M> **best, M *best_dist);

       public:
			/*!
			 * \brief add a node w meta data
			 */
        	void add_node_meta(const MathVector<dim>& vec, number m);

        	/*!
        	 * \brief add a node w/o meta data
        	 */
        	void add_node(const MathVector<dim>& vec);

           	/*!
        	 * \brief query the tree with a given vector for the nearest neighbor and return the attached metadata
        	 */
        	bool build_tree();

        	/*!
        	 * \brief query the tree for the nearest neighbor
        	 */
        	M query(const MathVector<dim, M>& vec);
         };
	} // namespace membrane_potential_mapping
} // namespace ug
//<! \}

#include "kd_tree_impl.h"

#endif // __H__UG__MEMBRANE_POTENTIAL_MAPPING__KD_TREE__
