/*!
 * \file vm2ug_rework_impl.h
 * \brief implementation for vm2ug rework
 *
 * \author Stephan Grein <stephan.grein@gcsc.uni-frankfurt.de>
 * \date 19/06/2015
 */

/// includes
#include <vector>
#include <utility>
#include <sstream>
#include <fstream>

#include "common/math/ugmath.h"

#include "vm2ug_rework.h"


/// using directives
using namespace ug::membrane_potential_mapping;

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
void Mapper<dim, M>::build_tree(const std::string& filename, const char& delim) {
	   std::ifstream file(filename.c_str());
	   UG_COND_THROW(!file, "Could not open data file: " + filename);

	   std::string line;
	   while (std::getline(file, line, delim)) {
		   MathVector<dim, number> node;
		   M meta;

	       std::istringstream iss(line);
	       for (size_t i = 0; i < dim; i++) {
	    	   iss >> node[i];
	       }

	       iss >> meta;
	       m_kdtree.add_node_meta(node, meta);
	   }
}

//////////////////////////////////////////////////////////
/// build_tree
//////////////////////////////////////////////////////////
template <size_t dim, typename M>
void Mapper<dim, M>::build_tree(const std::vector<std::pair<std::vector<number>, M> >& points) {
	for (CITVPVNM it = points.begin(); it != points.end(); ++it) {
		std::pair<std::vector<number>, M> pair = *it;

		MathVector<dim, number> coords;
		for (int i = 0; i < dim; i++) {
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
	for (int i = 0; i < dim; i++) {
		coords[i] = node.first[0];
	}

	m_kdtree.add_node_meta(coords, node.second);
}

//////////////////////////////////////////////////////////
/// get_data_from_nearest_neighbor
//////////////////////////////////////////////////////////
template <size_t dim, typename M>
M Mapper<dim, M>::get_data_from_nearest_neighbor(const MathVector<dim, number>& query) {
	return m_kdtree.query(query);
}

//////////////////////////////////////////////////////////
/// get_data_from_nearest_neighbor
//////////////////////////////////////////////////////////
template <size_t dim, typename M>
M Mapper<dim, M>::get_data_from_nearest_neighbor(const std::vector<number>& query) {
	MathVector<dim, number> coords;
	for (int i = 0; i < dim; i++) {
		coords[i] = query[i];
	}
	return m_kdtree.query(coords);
}
