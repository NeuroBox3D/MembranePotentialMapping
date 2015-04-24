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

/* \defgroup mpm_plugin Membrane Potential Mapping plugin
 * \ingroup plugins_experimental
 * \{
 */
namespace ug {
	namespace membrane_potential_mapping {
		/*!
		 * \brief simple kd node
		 */
		template <size_t dim, typename M>
        struct kd_node {
           number m_coords[dim];
           kd_node<dim, M>* left;
           kd_node<dim, M>* right;
           M* m_meta;
        };

		/*!
		 * \brief pretty simplistic KD Tree
		 */
		template <size_t dim, typename M>
		class KDTree {
		private:
        	std::vector<kd_node<dim, M> >nodes;
			int m_visited;

		public:
			/*!
			 * \brief ctor
			 */
			KDTree() : m_visited(0) {

			}

			/*!
			 * \brief dtor
			 */
			~KDTree() {

			}

			/*!
			 * \brief euclidean distance in dim-d
			 */
			number dist(kd_node<dim, M>* a, kd_node<dim, M> *b) const {
				number t, d = 0;
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
                number temp[dim];
                memcpy(temp,  x->m_coords, sizeof(temp));
                memcpy(x->m_coords, y->m_coords, sizeof(temp));
                memcpy(y->m_coords, temp,  sizeof(temp));
			}


			/*!
			 * \brief find median
			 */
			struct kd_node<dim, M>* find_median(kd_node<dim, M> *start,
				   kd_node<dim, M> *end, int idx) {
                if (end <= start) {
                        return NULL;
                }

                if (end == start + 1) {
                        return start;
                }

                kd_node<dim, M> *p, *store, *md = start + (end - start) / 2;
                number pivot;

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
			 * \brief generate a kd_tree from a list of kd_node structs
			 */
			kd_node<dim, M>* make_tree(kd_node<dim, M> *t, int len, int i) {
                kd_node<dim, M>* n;

                if (!len) return 0;

                if ((n = find_median(t, t + len, i))) {
                        i = (i + 1) % dim;
                        n->left  = make_tree(t, n - t, i);
                        n->right = make_tree(n + 1, t + len - (n + 1), i);
                }
                return n;
			}

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
                     kd_node<dim, M> **best, number *best_dist) {
                number d;
                number dx;
                number dx_squared;

                if (!root) {
                	return;
                }

                d = dist(root, nd);
                dx = root->m_coords[i] - nd->m_coords[i];
                dx_squared = std::pow(d, 2);

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

#ifdef UG_DIM_3
        	void add_node_with_meta(const MathVector<dim>& vec, M* m) {
        		/// add a node
        		struct kd_node<dim, M> node;

        		for (int i = 0; i < dim-1; i++) {
        			node.m_coords[0] = vec[i];
        		}

        		node.m_meta = m;

        		this->nodes.push_back(node);
        	}

        	bool build_tree() {
        		struct kd_node<dim, M>* root;
        		root = make_tree(this->nodes, sizeof(this->nodes) / sizeof(this->nodes[0]), 0);
        		return (root != NULL);
        	}

        	M query_node(const MathVector<dim>& vec) {
        		/// query the tree
        		struct kd_node<dim, M> query;
        		for (int i = 0; i < dim-1; i++) {
        			query.m_coords[0] = vec[i];
        		}

        		struct kd_node<dim, M> *root, *found;

        		number best_dist;
        		root = make_tree(this->nodes, sizeof(this->nodes) / sizeof(this->nodes[0]), 0);
        		nearest(root, &query, 0, &found, &best_dist);

        		/// and return vm
        		return *found->m_meta;
        	}
#endif
        };
	} // namespace synapse_provider
} // namespace ug
//<! \}
#endif // __H__UG__MEMBRANE_POTENTIAL_MAPPING__KDTREE__
