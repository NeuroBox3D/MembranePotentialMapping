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
			int m_visited;
        	std::vector<kd_node<dim, M> > nodes;
        	kd_node<dim, M>* root;
    		kd_node<dim, number>* wps;

		public:
			/*!
			 * \brief ctor
			 */
			KDTree() : m_visited(0), root(NULL) {

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
                     kd_node<dim, M> **best, number *best_dist) {
                number d;
                number dx;
                number dx_squared;

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
			 * \brief add a vector with meta data to the tree
			 */
        	void add_vec_with_meta(const MathVector<dim>& vec, M* m) {
        		struct kd_node<dim, M> node;

        		for (int i = 0; i < dim; i++) {
        			node.m_coords[0] = vec[i];
        		}

        		node.m_meta = m;

        		this->nodes.push_back(node);
        	}


        	/*!
        	 * \brief builds the tree from the stored nodes
        	 */
        	bool build_tree() {
        		root = make_tree(this->nodes, sizeof(this->nodes) / sizeof(this->nodes[0]), 0);
        		return (root != NULL);
        	}

        	/*!
        	 * \brief query the tree with a given vector for the nearest neighbor and return the attached metadata
        	 */
        	M* build_tree(const MathVector<dim>& vec) {
        		struct kd_node<dim, M> query;
        		for (int i = 0; i < dim; i++) {
        			std::cout << "there" << std::endl;
        			query.m_coords[0] = vec[i];
        		}

        		std::cout << "Here" << std::endl;

              		int i;
        		/// if wps go out of scope we have memory corruption -> aka seg faults, store tree within the class
        			kd_node<dim, number> wp[this->nodes.size()];
        			kd_node<dim, number> a, b, q;
        		for (int i = 0; i < dim; i++) {
   //     			a.m_coords[i] = 0;
 //       			b.m_coords[i] = 1;
        			q.m_coords[i] = vec[i];
        		}
 //       		wp[0] = a;
   //     		wp[1] = b;
       			this->wps = (kd_node<dim, number>*)malloc(sizeof(kd_node<dim, number>)*nodes.size());
        			for (int i = 0; i < nodes.size(); i++) {
        				wp[i] = nodes[i];
        				this->wps[i] = wp[i];
        			}

       		//	this->wps[0] = a;
       			//this->wps[1] = b;

        			kd_node<dim, number> *found;
        			double best_dist;

 //       			root = make_tree(this->wps, sizeof(this->wps) / sizeof(this->wps[1]), 0);
       			root = make_tree(this->wps, nodes.size(), 0);
        			found = 0;
        			nearest(root, &q, 0, &found, &best_dist);
        			std::cout << "best_distance: " << best_dist;

        			found = 0;
        			nearest(root, &q, 0, &found, &best_dist);
        			std::cout << "best_distance: " << best_dist;

        			found = 0;
        			nearest(root, &q, 0, &found, &best_dist);
        			std::cout << "best_distance: " << best_dist;

        			query_vec_only(vec);

        	}

        	void query_vec_only(const MathVector<dim, number>& vec) {
        		struct kd_node<dim, M> query;
        		        		for (int i = 0; i < dim; i++) {
        		        			std::cout << "there" << std::endl;
        		        			query.m_coords[i] = vec[i];
        		        		}
        			kd_node<dim, number> *found;
        			double best_dist = 0;
        		found = 0;
        		if (!root) {
        			std::cout << "root not available!" << std::endl;
        		}
        			nearest(root, &query, 0, &found, &best_dist);
        			std::cout << "best_distance: " << best_dist;
        	}
         };
	} // namespace synapse_provider
} // namespace ug
//<! \}
#endif // __H__UG__MEMBRANE_POTENTIAL_MAPPING__KDTREE__
