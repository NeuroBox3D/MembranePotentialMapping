/*!
 * \file kd_tree_impl.h
 *
 *  Created on: Apr 28, 2015
 *      Author: stephan
 */

/// guard
#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__KD_TREE_IMPL__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__KD_TREE_IMPL__

/// includes
#include <cstdlib>
#include <common/log.h>
#include "kd_node.h"
#include "kd_tree.h"

namespace ug {
	namespace membrane_potential_mapping {
	    /////////////////////////////////////////////////////////
		/// ctor
	    /////////////////////////////////////////////////////////
		template <size_t dim, typename M>
		kd_tree<dim, M>::kd_tree() : m_visited(0), m_pRoot(NULL), m_pWps(NULL) {
		}

	    /////////////////////////////////////////////////////////
		/// dtor
	    /////////////////////////////////////////////////////////
		template <size_t dim, typename M>
		kd_tree<dim, M>::~kd_tree() {
			delete m_pRoot;
			delete m_pWps;
		}

	    /////////////////////////////////////////////////////////
		/// dist
	    /////////////////////////////////////////////////////////
		template <size_t dim, typename M>
		M kd_tree<dim, M>::dist(kd_node<dim, M>* a, kd_node<dim, M> *b) const {
			M t, d = 0;
            size_t dimension = dim;
            while (dimension--) {
                    t = a->m_coords[dimension] - b->m_coords[dimension];
                    d += t * t;
            }
            return d;
		}


	    /////////////////////////////////////////////////////////
		/// swap
	    /////////////////////////////////////////////////////////
		template <size_t dim, typename M>
		void kd_tree<dim, M>::swap(kd_node<dim, M>* x, kd_node<dim, M>* y) {
                M temp[dim];
                memcpy(temp,  x->m_coords, sizeof(temp));
                memcpy(x->m_coords, y->m_coords, sizeof(temp));
                memcpy(y->m_coords, temp,  sizeof(temp));
		}

	    /////////////////////////////////////////////////////////
		/// find_median
	    /////////////////////////////////////////////////////////
		template <size_t dim, typename M>
		kd_node<dim, M>* kd_tree<dim, M>::find_median(kd_node<dim, M> *start, kd_node<dim, M> *end, int idx) {
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

	    /////////////////////////////////////////////////////////
		/// make_tree
	    /////////////////////////////////////////////////////////
		template <size_t dim, typename M>
		kd_node<dim, M>* kd_tree<dim, M>::make_tree(kd_node<dim, M>* t, int len, int i) {
           kd_node<dim, M>* n;
           if (!len) return 0;

           if ((n = find_median(t, t + len, i))) {
        	   i = (i + 1) % dim;
               n->left  = make_tree(t, n - t, i);
               n->right = make_tree(n + 1, t + len - (n + 1), i);
            }

            return n;
		}

	    /////////////////////////////////////////////////////////
		/// make_tree
	    /////////////////////////////////////////////////////////
		template <size_t dim, typename M>
		void kd_tree<dim, M>::nearest(kd_node<dim, M> *root, kd_node<dim, M> *nd, int i, kd_node<dim, M> **best, M *best_dist) {
           M d;
           M dx;
           M dx_squared;

           if (!root) return;

           d = dist(root, nd);
           dx = root->m_coords[i] - nd->m_coords[i];
           dx_squared = dx * dx;

           m_visited++;

           if (!*best || d < *best_dist) {
        	   *best_dist = d;
        	   *best = root;
           }

           if (!*best_dist) return;

           if (++i >= dim) i = 0;

           nearest(dx > 0 ? root->left : root->right, nd, i, best, best_dist);
           if (dx_squared >= *best_dist) return;
           nearest(dx > 0 ? root->right : root->left, nd, i, best, best_dist);
        }

	    /////////////////////////////////////////////////////////
		/// add_node_meta
	    /////////////////////////////////////////////////////////
		template <size_t dim, typename M>
        void kd_tree<dim, M>::add_node_meta(const MathVector<dim>& vec, number m) {
        	kd_node<dim, M> node;
        	for (size_t i = 0; i < dim; i++) node.m_coords[i] = vec[i];
        	node.m_meta = m;
        	this->nodes.push_back(node);
        }

	    /////////////////////////////////////////////////////////
		/// add_node
	    /////////////////////////////////////////////////////////
		template <size_t dim, typename M>
        void kd_tree<dim, M>::add_node(const MathVector<dim>& vec) {
        	kd_node<dim, M> node;
        	for (size_t i = 0; i < dim; i++) node.m_coords[i] = vec[i];
        	this->nodes.push_back(node);
		}

	    /////////////////////////////////////////////////////////
		/// build_tree
	    /////////////////////////////////////////////////////////
		template <size_t dim, typename M>
        bool kd_tree<dim, M>::build_tree() {
       		/// allocate nodes
       		size_t numNodes = nodes.size();
       		this->m_pWps = (kd_node<dim, M>*)malloc(sizeof(kd_node<dim, M>)*numNodes);
       		for (size_t i = 0; i < numNodes; i++) this->m_pWps[i] = nodes[i];
       		this->nodes.clear();

       		/// make the tree and indicate success or not
       		m_pRoot = make_tree(this->m_pWps, numNodes, 0);
       		if (m_pRoot) return true;
       		else return false;
       	}

	    /////////////////////////////////////////////////////////
		/// query
	    /////////////////////////////////////////////////////////
		template <size_t dim, typename M>
        M kd_tree<dim, M>::query(const MathVector<dim, M>& vec) {
        	/// meta value
        	M meta = 0;

        	/// check if tree has been build
        	if (!m_pRoot) {
        		UG_DLOGN(MPM_KDTREE, 0, "KDTree not build for dim=" << dim);
        	} else {
        		/// build query
        		kd_node<dim, M> query;
        		for (size_t i = 0; i < dim; i++) query.m_coords[i] = vec[i];

        		/// query tree
        		kd_node<dim, M>* found = NULL;
        		M best_dist = 0;
        		nearest(m_pRoot, &query, 0, &found, &best_dist);

        		/// debug output
        		UG_DLOGN(MPM_KDTREE, 0, "Best distance for query [" << query << "]: " << std::sqrt(best_dist));
        		UG_DLOGN(MPM_KDTREE, 0, "Found point [" << found << "]";)
        		meta = found->m_meta;
        	}

        	/// return the meta value
        	return meta;
        }
	} // namespace membrane_potential_mapping
} // namespace ug
//<! \}

#endif /// __H__UG__MEMBRANE_POTENTIAL_MAPPING__KD_TREE_IMPL__
