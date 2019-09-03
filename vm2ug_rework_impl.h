/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Author: Stephan Grein
 * Creation date: 2015-06-19
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

/// includes
#include <sstream>
#include <fstream>

#include "boost/lexical_cast.hpp"
#include "vm2ug.h"


namespace ug {
namespace membrane_potential_mapping {


//////////////////////////////////////////////////////////
/// build_tree
//////////////////////////////////////////////////////////
template <size_t dim, typename M>
void Mapper<dim, M>::build_tree() {
	m_kdtree.build_tree();
}

//////////////////////////////////////////////////////////
/// build_tree
//////////////////////////////////////////////////////////
template <size_t dim, typename M>
void Mapper<dim, M>::build_tree(const std::vector<std::pair<MathVector<dim, number>, M> >& points) {
	m_kdtree = kd_tree<dim, M>();
	for (CITVPMNM it = points.begin(); it != points.end(); ++it) {
		std::pair<MathVector<dim, number>, M> pair = *it;
		m_kdtree.add_node_meta(pair.first, pair.second);
	}
	m_kdtree.build_tree();
}

//////////////////////////////////////////////////////////
/// build_tree
//////////////////////////////////////////////////////////
template <size_t dim, typename M>
void Mapper<dim, M>::build_tree(const std::string& filename) {
	m_kdtree = kd_tree<dim, M>();

	std::ifstream file(filename.c_str());
	UG_COND_THROW(!file, "Could not open data file: " + filename);

	std::string line;
	size_t lineNo = 0;
	while (std::getline(file, line)) {
		++lineNo;

		/// read data
		MathVector<dim, number> node(0.0);
		std::istringstream iss(line);

		if (! (iss >> node[0]))
			break;

		for (size_t i = 1; i < dim; ++i)
			if ((! (iss >> node[i])))
				UG_THROW("Reading coordinate " << i << " did not succeed on line " << lineNo << ".");

		M meta;
		if (! (iss >> meta))
			UG_THROW("Reading meta data did not succeed on line " << lineNo << ".");

		/// add the node
		m_kdtree.add_node_meta(node, meta);
	}

	/// finally build the tree
	if (!m_kdtree.build_tree())
	   UG_THROW("KD tree could not be built.");
}

//////////////////////////////////////////////////////////
/// build_tree
//////////////////////////////////////////////////////////
template <size_t dim, typename M>
void Mapper<dim, M>::build_tree(const std::vector<std::pair<std::vector<number>, M> >& points) {
	m_kdtree = kd_tree<dim, M>();
	for (CITVPVNM it = points.begin(); it != points.end(); ++it) {
		std::pair<std::vector<number>, M> pair = *it;

		MathVector<dim, number> coords;
		for (size_t i = 0; i < dim; i++) {
			coords[i] = pair.first[0];
		}

		m_kdtree.add_node_meta(coords, pair.second);
	}
	m_kdtree.build_tree();
}

//////////////////////////////////////////////////////////
/// add_node
//////////////////////////////////////////////////////////
template <size_t dim, typename M>
void Mapper<dim, M>::add_node(const std::pair<MathVector<dim, number>, M>& node) {
	m_kdtree.add_node_meta(node.first, node.second);
}

//////////////////////////////////////////////////////////
/// add_node
//////////////////////////////////////////////////////////
template <size_t dim, typename M>
void Mapper<dim, M>::add_node(const std::pair<std::vector<number>, M>& node) {
	MathVector<dim, number> coords;

	for (size_t i = 0; i < dim; i++) {
		coords[i] = node.first[0];
	}

	m_kdtree.add_node_meta(coords, node.second);
}

//////////////////////////////////////////////////////////
/// add_node
//////////////////////////////////////////////////////////
template <size_t dim, typename M>
void Mapper<dim, M>::add_node(const std::vector<number>& node, const M& data) {
	add_node(std::make_pair(node, data));
}

//////////////////////////////////////////////////////////
/// get_data_from_nearest_neighbor
//////////////////////////////////////////////////////////
template <size_t dim, typename M>
M Mapper<dim, M>::get_data_from_nearest_neighbor(const MathVector<dim, number>& query) const {
	return m_kdtree.query(query);
}

//////////////////////////////////////////////////////////
/// get_data_from_nearest_neighbor
//////////////////////////////////////////////////////////
template <size_t dim, typename M>
M Mapper<dim, M>::get_data_from_nearest_neighbor(const std::vector<number>& query) const {
	MathVector<dim, number> coords;
	for (size_t i = 0; i < dim; i++) {
		coords[i] = query[i];
	}
	return m_kdtree.query(coords);
}

} // end namespace mpm
} // end namespace ug

