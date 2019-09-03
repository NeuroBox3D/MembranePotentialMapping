/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2015-04-28
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/// guard
#ifndef UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__KDTREE__KD_TREE_H
#define UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__KDTREE__KD_TREE_H

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
		//int m_visited;
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
		void nearest(kd_node<dim, M> *root, kd_node<dim, M> *nd, size_t i, kd_node<dim, M> **best, M *best_dist) const;

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
		 * @brief print tree
		 */
		void print_tree(kd_node<dim, M>* root, size_t lvl) const;

		/*!
		 * \brief query the tree for the nearest neighbor
		 * \param[in] vec query node
		 * \param[out] M meta data
		 */
		M query(const MathVector<dim, M>& vec) const;
};

} // namespace membrane_potential_mapping
} // namespace ug
//<! \}

#include "kd_tree_impl.h"

#endif // UG__PLUGINS__MEMBRANE_POTENTIAL_MAPPING__KDTREE__KD_TREE_H
