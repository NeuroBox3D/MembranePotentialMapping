/*!
 * \file vm2ug_rework_impl.h
 * \brief implementation for vm2ug rework
 *
 * \author Stephan Grein <stephan.grein@gcsc.uni-frankfurt.de>
 * \date 19/06/2015
 */

/// includes
#include <sstream>
#include <fstream>

#include "vm2ug_rework.h"

#include "boost/lexical_cast.hpp"


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

