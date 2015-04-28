/*!
 * \file kdtree.h
 *
 *  Created on: Apr 21, 2015
 *      Author: stephan
 */

/// guard
#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__KDTREE__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__KDTREE__

/// includes
#include <cstdlib>
#include <common/log.h>

////////////////////////////////////////////////////////////////////////
/// debug ids
////////////////////////////////////////////////////////////////////////
ug::DebugID MPM_KDTREE("MPM_DID.kdtree");

/* \defgroup mpm_plugin Membrane Potential Mapping plugin
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
	namespace membrane_potential_mapping {
		/*!
		 * \brief simple kd node with meta attachment
		 */
		template <size_t dim, typename M>
        struct kd_node {
           M m_coords[dim];
           M m_meta;
           kd_node<dim, M>* left;
           kd_node<dim, M>* right;

           /*!
            * \brief pretty print point coordinates with meta
            * \function operator<<
            */
           friend std::ostream& operator<<(std::ostream &out, const kd_node& node) {
        	   for (size_t i = 0; i < dim-1; i++) {
        		   out << node.m_coords[i] << ", ";
        	   }
        	   out << node.m_coords[dim] << "(" << *node.m_meta << ")";
        	   return out;
           }
        };

		/*!
		 * \brief pretty simplistic KD Tree
		 */
		template <size_t dim, typename M>
		class KDTree {
		private:
			int m_visited;
        	std::vector<kd_node<dim, M> > nodes;
        	kd_node<dim, M>* root;
    		kd_node<dim, M>* wps;

		public:
			/*!
			 * \brief ctor
			 */
			KDTree() : m_visited(0), root(NULL), wps(NULL) {

			}

			/*!
			 * \brief dtor
			 */
			~KDTree() {

			}

		private:
			/*!
			 * \brief euclidean distance in dim-d
			 */
			M dist(kd_node<dim, M>* a, kd_node<dim, M> *b) const {
				M t, d = 0;
                size_t dimension = dim;
                while (dimension--) {
                        t = a->m_coords[dimension] - b->m_coords[dimension];
                        d += t * t;
                }
                return d;
			}

			/*!
			 * \brief swap nodes
			 */
			void swap(kd_node<dim, M>* x, kd_node<dim, M>* y) {
                M temp[dim];
                memcpy(temp,  x->m_coords, sizeof(temp));
                memcpy(x->m_coords, y->m_coords, sizeof(temp));
                memcpy(y->m_coords, temp,  sizeof(temp));
			}

			/*!
			 * \brief find median
			 */
			kd_node<dim, M>* find_median(kd_node<dim, M> *start, kd_node<dim, M> *end, int idx) {
                if (end <= start) {
                        return NULL;
                }

                if (end == start + 1) {
                        return start;
                }

                kd_node<dim, M> *p, *store, *md = start + (end - start) / 2;
                M pivot;

                while (1) {
                        pivot = md->m_coords[idx];

                        swap(md, end - 1);

                        for (store = p = start; p < end; p++) {
                                if (p->m_coords[idx] < pivot) {
                                        if (p != store) {
                                                swap(p, store);
                                        }
                                        store++;
                                }
                        }

                        swap(store, end - 1);

                        if (store->m_coords[idx] == md->m_coords[idx]) {
                                return md;
                        }

                        if (store > md)	{
                        	end = store;
                        } else {
                        	start = store;
                        }

                }
                return NULL;
			}

			/*!
			 * \brief generates a kd tree from a list of kd nodes
			 */
			kd_node<dim, M>* make_tree(kd_node<dim, M>* t, int len, int i) {
                kd_node<dim, M>* n;

                if (!len) return 0;

                if ((n = find_median(t, t + len, i))) {
                        i = (i + 1) % dim;
                        n->left  = make_tree(t, n - t, i);
                        n->right = make_tree(n + 1, t + len - (n + 1), i);
                }
                return n;
			}

			/*!
			 * \brief generates a kd tree from a vector of kd nodes
			 */
			kd_node<dim, M>* make_tree(const std::vector<kd_node<dim, M> >& nodes, int len, int i) {
				kd_node<dim, M> snodes[nodes.size()];
				int count = 0;
				for (typename std::vector<kd_node<dim, M> >::const_iterator cit = nodes.begin(); cit != nodes.end(); cit++) {
					snodes[count] = *cit;
					count++;
				}
				return make_tree(snodes, len, i);
			}

			/*!
			 * \brief get nearest neighbor in tree wrt to a query point (euclidean distance is used)
			 */
			void nearest(kd_node<dim, M> *root, kd_node<dim, M> *nd, int i,
                     kd_node<dim, M> **best, M *best_dist) {
                M d;
                M dx;
                M dx_squared;

                if (!root) {
                	return;
                }

                d = dist(root, nd);
                dx = root->m_coords[i] - nd->m_coords[i];
                dx_squared = dx * dx;

                m_visited++;

                if (!*best || d < *best_dist) {
                        *best_dist = d;
                        *best = root;
                }

                if (!*best_dist) {
                	return;
                }

                if (++i >= dim) {
                	i = 0;
                }

                nearest(dx > 0 ? root->left : root->right, nd, i, best, best_dist);

                if (dx_squared >= *best_dist) {
                	return;
                }

                nearest(dx > 0 ? root->right : root->left, nd, i, best, best_dist);
        	}

		public:
			/*!
			 * \brief add a node w meta data
			 */
        	void add_node_meta(const MathVector<dim>& vec, M m) {
        		kd_node<dim, M> node;

        		for (int i = 0; i < dim; i++) {
        			node.m_coords[i] = vec[i];
        		}

        		node.m_meta = m;

        		this->nodes.push_back(node);
        	}

        	/*!
        	 * \brief add a node w/o meta data
        	 */
        	void add_node(const MathVector<dim>& vec) {
        		kd_node<dim, M> node;

        		for (int i = 0; i < dim; i++) {
        			node.m_coords[i] = vec[i];
        		}

        		this->nodes.push_back(node);
        	}

           	/*!
        	 * \brief query the tree with a given vector for the nearest neighbor and return the attached metadata
        	 */
        	bool build_tree() {
        		kd_node<dim, M> wp[this->nodes.size()];

       			this->wps = (kd_node<dim, M>*)malloc(sizeof(kd_node<dim, M>)*nodes.size());

        		for (int i = 0; i < nodes.size(); i++) {
        			wp[i] = nodes[i];
        			this->wps[i] = wp[i];
        		}

       			root = make_tree(this->wps, nodes.size(), 0);
       			if (root) return true;
       			else return false;
        	}

        	/*!
        	 * \brief query the tree for the nearest neighbor
        	 */
        	M query(const MathVector<dim, M>& vec) {
        		/// meta value
        		M meta = 0;

        		/// check if tree has been build
        		if (!root) {
        			UG_DLOGN(MPM_KDTREE, 0, "KDTree not build for dim=" << dim);
        		} else {
        			/// build query
        			kd_node<dim, M> query;
        			for (int i = 0; i < dim; i++)
        				query.m_coords[i] = vec[i];

        			/// query tree
        			kd_node<dim, M>* found = NULL;
        			M best_dist = 0;
        			nearest(root, &query, 0, &found, &best_dist);

        			/// debug output
        			UG_DLOGN(MPM_KDTREE, 0, "Best distance for query [" << query << "]: " << std::sqrt(best_dist));
        			UG_DLOGN(MPM_KDTREE, 0, "Found point [" << found << "]";)
        			meta = found->m_meta;
        		}

        		/// return the meta value
        		return meta;
        	}
         };
	} // namespace synapse_provider
} // namespace ug
//<! \}
#endif // __H__UG__MEMBRANE_POTENTIAL_MAPPING__KDTREE__
