/*!
 * \file vm2ug.h
 * \brief header for the membrane_potential_mapping
 * \addtogroup mpm_plugin
 * \see docs for additional information
 *
 * \author Stephan Grein
 * \date July, 2011
 *
 * Notes:
 * TODO: Migrate to the UG internal kd-tree for dismissing the ANNLIB
 * TODO: Code optimization	and cleanup (cctor calls, typename vs. class for templates)
 * TODO: introduce template parameter for DIM in mapper -- done!
 * TODO: code cleanup
 */

#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG__


#include <vector>
#include <string>
#include <istream>
#include <ostream>
#include <cmath>

#include "common/common.h"
#include "common/util/smart_pointer.h"

#include "kdtree/kd_tree.h"
// FIXME: below fails if EMBEDDED_PLUGINGS=ON -> why?
// comment: hmm, works fine for me (mbreit)
#include "ANN/ANN.h"

#include "common_typedefs.h"
#include "mvec.h"
#include "transformator.h"

namespace ug {
namespace membrane_potential_mapping {



// forward declarations
template <class T> class Vm2uG;
class sPoint;
class uGPoint;

template <size_t dim, typename M> class kd_tree;

template <class T> std::ostream& operator<< (std::ostream& output, const Vm2uG<T>& p);
std::ostream& operator<< (std::ostream& output, const sPoint& p);
std::ostream& operator<< (std::ostream& output, const uGPoint& p);



/*!
 * \brief Vm2uG represents a class to perform the k nearest neighbor search with a kd tree and can return the neighbors encapsulates in classes
 * The template parameter T encodes the type used for specification of time steps (int, number, std::string).
 *
 * \addtogroup mpm_plugin
 *
 * \see docs for additional informations
 */
template <class T> class Vm2uG {

	friend std::ostream& operator<< <>(std::ostream& output, const Vm2uG<T>& p);

	public:
		/*!
		 * \brief main constructor
		 *
		 * \param[in] dataFileBaseName basename of the input files (timesteps of the Neuron simulation)
		 * \param[in] dataFileExt extension of the input files (i.e. .csv)
		 * \param[in] promise_ specify if ordering of the datapoints in the input files can change (if true, tree needs to be rebuild less often)
		 */
		Vm2uG(const std::string& dataFileBaseName=std::string("timesteps/"), const std::string& dataFileExt=std::string(".csv"), bool promise_ = false);

		/*!
		 * \brief enhanced constructor
		 *
		 * \param[in] dataFileBaseName basename of the input files (timesteps of the Neuron simulation)
		 * \param[in] dim_file dimension of the datapoints in the input file (default: 3)
		 * \param[in] maxPts maximum number of datapoints in the input files (default: 10000)
		 * \param[in] eps approximation factor for nearest neighbor search (default: 0)
		 * \param[in] k search for k nearest neighbors (default: 1)
		 */
		Vm2uG(const std::string& dataFileBaseName, size_t dim_file, size_t maxPts, number eps, size_t k);


#ifdef MPMNEURON
		/*!
		 * \brief another constructor for neuron usage
		 */

		Vm2uG(SmartPtr<Transformator> transformator) {
			m_transformator = transformator;
			isTreeBuilt = false;
			promise = true;
			dim_file = 3;
			timestep = 0;
			maxPts = 100000;
			k = 1;
			eps = 0.0;
			buildTree();
		}
#endif

		/*!
		 * \brief default destructor
		 */
		~Vm2uG();


		/*!
		 * \brief return the associated membrane potential of the first nearest neighbor at Cartesian coordinates x,y,z at timestep t: special case DIM=3.
		 *
		 * \param[in] x coordinate x
		 * \param[in] y coordinate y
		 * \param[in] z coordinate z
		 * \param[in] t timestep
		 *
		 * @deprecated This method is deprecated and might be removed in future. Use
		 * 			   template <size_t dim>
		 * 			   number get_potential(const MathVector<dim>& coords, const T& t)
		 * 			   instead.
		 * \return \c number the membrane potential
		 */
		number get_potential(number x, number y, number z, const T& t) {
			MathVector<3> node;
			node[1] = x; node[2] = y; node[3] = z;
			return get_potential(node, t);
		}

		/*!
		 * \brief Returns nearest neighbor potential to query point coordinates.
		 * Returns the associated membrane potential of the nearest neighbor
		 * at coordinates coords and time step t.
		 *
		 * \param[in] coords coordinates
		 * \param[in] t time step
		 *
		 * \return \c number the membrane potential
		 */
		template <size_t dim>
		number get_potential(const MathVector<dim>& coords, const T& t)
		{
			return vm_t(t, coords).getVm();
		}

#ifdef MPMNEURON
		number get_potential(number x, number y, number z) {
			number node[] = {x, y, z};
			return vm_t(node);
		}

		number get_potential(const number coords[])
		{
			return vm_t(coords);
		}
#endif



		/*!
		 * \brief return the associated membrane potential of the first nearest neighbor at cartesian coordinates x,y,z at timestep t: special case DIM=3.
		 *
		 * \param[in] x coordinate x
		 * \param[in] y coordinate y
		 * \param[in] z coordinate z
		 * \param[in] t timestep
		 * \param[in] cutoff
		 * \param[in] k
		 *
		 * @deprecated This method is deprecated and might be removed in future. Use
		 * 			   template <size_t dim>
		 * 			   number get_potential_lin(const MathVector<dim>& coords, const T& t, const number cutoff, const int k)
		 * 			   instead.
		 * \return \c number the membrane potential (linearly interpolated)
		 */
		number get_potential_lin(number x, number y, number z, const T& t, number cutoff, size_t k=2) {
			MathVector<3> node;
			node[1] = x; node[2] = y; node[3] = z;
			return get_potential_lin(node, t, cutoff, k);
		}

		/*!
		 * \brief Returns linearly interpolated k-NN potential to query point coordinates.
		 * Returns the linearly interpolated associated membrane potential of the k nearest neighbors
		 * at coordinates coords and time step t.
		 *
		 * \param[in] coords 	coordinates
		 * \param[in] t 		timestep
		 * \param[in] cutoff
		 * \param[in] k
		 *
		 * \return \c number the membrane potential (linearly interpolated)
		 */
		template <size_t dim>
		number get_potential_lin(const MathVector<dim>& coords, const T& t, number cutoff, size_t k=2)
		{
			return interp_lin_vms(t, coords, cutoff, k);
		}

		/*!
		 * \brief return the associated membrane potential of the first nearest neighbor at cartesian coordinates x,y,z at timestep t: special case DIM=3.
		 *
		 * \param[in] x coordinate x
		 * \param[in] y coordinate y
		 * \param[in] z coordinate z
		 * \param[in] t timestep
		 * \param[in] cutoff
		 * \param[in] k
		 *
		 * @deprecated This method is deprecated and might be removed in future. Use
		 * 			   template <size_t dim>
		 * 			   number get_potential_bilin(const MathVector<dim>& coords, const T& t, const number cutoff, const int k)
		 * 			   instead.
		 *
		 * \return \c number the membrane potential (bilinearly interpolated)
		 */
		number get_potential_bilin(number x, number y, number z, const T& t, number cutoff, size_t k=4) {
			MathVector<3> node;
			node[1] = x; node[2] = y; node[3] = z;
			return get_potential_bilin(node, t, cutoff, k);
		}

		/*!
		 * \brief Returns bi-linearly interpolated k-NN potential to query point coordinates.
		 * Returns the bi-linearly interpolated associated membrane potential of the k nearest neighbors
		 * at coordinates coords and time step t.
		 *
		 * \param[in] coords 	coordinates
		 * \param[in] t 		timestep
		 * \param[in] cutoff
		 * \param[in] k
		 *
		 * \return \c number the membrane potential (bi-linearly interpolated)
		 */
		template <size_t dim>
		number get_potential_bilin(const MathVector<dim>& coords, const T& t, number cutoff, size_t k=2)
		{
			return interp_bilin_vms(t, coords, cutoff, k);
		}

#ifdef MPMNEURON
		void buildTree();
#endif

		/*!
		 * \brief builds the tree with timestep
		 *
		 * \param[in] timestep the timestep for build-up of the kd-tree
		 * \return \c void
		 */
		void buildTree(const T& timestep);


		/*!
		 * \brief get the nearest neighbor and return the associated membrane potential and distance to the query point
		 *
		 * \param[in] timestep the timestep for which a membrane potential shall be searched
		 * \param[in] node the query point (cartesian coordinate of grid point from UG)
		 *
		 * \see uGPoint
		 * \see sPoint
		 *
		 * \return \c uGPoint an grid point from UG with additional attributes, e. g. membrane potential, distance to query point, etc.
		 */

	private:
		template <size_t dim>
		uGPoint vm_t(const T& timestep, const MathVector<dim>&);

#ifdef MPMNEURON
		number vm_t(number node[]);
#endif

		/*!
		 * \brief the same as vm_t except supplied is a list of query points; parallelized with OMP.
		 *
		 * \param[in] timestep the timestep for which a mmebrane potential shall be searched
		 * \param[in] nodes the query points (cartesian coordinates for grid points from UG)
		 *
		 * \return \c std::vector<uGPoint> grid points from UG with additional attributes, e. g. membrane potential, distance to query points, etc.
		 *
		 * \see QDIM
		 * \see vm_t
		 * \see sPoint
		 * \see uGPoint
		 *
		 */
		template <size_t dim>
		std::vector<uGPoint> vm_t_many(const T& timestep, const std::vector<MathVector<dim> >& nodes);


		/*!
		 * \brief the same as vm_t but get the k nearest neighbors of one query points
		 *
		 * \param[in] timestep the time step for which a membrane potential shall be searched
		 * \param[in] node the query point (Cartesian coordinate of grid point from UG)
		 * \param[in] k how many nearest neighbors should be searched
		 *
		 * \return \c uGPoint an grid point from UG with additional attributes, e. g. membrane potentials, distance to query points, etc.
		 * \see vm_t
		 * \see uGPoint
		 * \see sPoint
		 */
		template <size_t dim>
		uGPoint vm_t_k(const T& timestep, const MathVector<dim>& node, const size_t k);

		/*!
		 * \brief the same as vm_t_many but get the k nearest neighbors of many query points
		 *
		 * \param[in] timestep the time step for which a membrane potential shall be searched
		 * \param[in] nodes the query points (Cartesian coordinates for grid points from UG)
		 * \param[in] k how many nearest neighbors should be searched
		 *
		 * \return \c std::vector<uGPoint> grid points from UG with additional attributes, e. g. membrane potentials, distance to query points, etc.
		 *
		 * \see QDIM
		 * \see vm_t_many
		 * \see sPoint
		 * \see uGPoint
		 */
		template <size_t dim>
		std::vector<uGPoint> vm_t_many_k(const T& timestep, const std::vector<MathVector<dim> >& nodes, size_t k);

		/*!
		 * \brief interpolate linearly the membrane potentials
		 *
		 * \param[in] node the query point (uG)
		 * \param[in] timestep for which the nearest neighbor should be computed
		 * \param[in] cutoff maximum distance for nearest neighbor
		 * \param[in] k how many nearest neighbors should be considered for linear interpolation
		 * \return \c number the membrane potential (potentially linearly interpolated)
		 */
		template <size_t dim>
		number interp_lin_vms(const T& timestep, const MathVector<dim>& coords, number cutoff, size_t k);

		/*!
		 * \brief interpolate bilinearly the membrane potentials
		 *
		 * \param[in] node the query point (uG)
		 * \param[in] timestep for which the nearest neighbor should be computed
		 * \param[in] cutoff maximum distance for nearest neighbor
		 * \param[in] k how many nearest neighbors should be considered for bilinear interpolation
		 *
		 * \return \c number the membrane potential (potentially bilinearly interpolated)
		 */
		template <size_t dim>
		number interp_bilin_vms(const T& timestep, const MathVector<dim>& coords, number cutoff, size_t k);

	public:
		// setters
		void setK(size_t k);
		void setDim(size_t dim_file);
		void setTimestep(const T& timestep);
		void setMaxPts(size_t maxPts);
		void setEps(number eps);
		void setPromise(bool promise);
		void setdataFileBaseName(const std::string& dataFileBaseName);
		void setdataFileExt(const std::string& dataFileExt);

		// getters
		bool treeBuilt() { return isTreeBuilt; }

	protected:
		bool promise;
		bool isTreeBuilt;

		size_t k;
		size_t dim_file;

		size_t maxPts;

		long timestep;	// hash value encoding the current time step
		T ts_oType;		// time step in original type

		number eps;

		std::string dataFileBaseName;
		std::string dataFileExt;

	private:
		size_t qPts;
		size_t nPts;

		std::istream* queryIn;
		std::istream* dataIn;

		std::vector<ANNidxArray> annQueryPtsIdxArray;
		std::vector<ANNdistArray> annQueryPtsDistsArray;

		ANNpointArray dataPts;
		ANNpoint queryPt;

		ANNidxArray nnIdx;
		ANNdistArray dists;
		ANNkd_tree* kdTree;

		ANNpointArray queryPts;

		static std::ifstream dataStream;

#ifdef MPMNEURON
		SmartPtr<Transformator> m_transformator;
#endif

		/*!
		 * \brief prints a data or query point with its cartesian coordinates
		 *
		 * \param[in] out where the output should go to (default: std::cout)
		 * \param[in] p an ANN point
		 * \param[in] newline if true prints newline otherwise newline is omitted
		 *
		 * \return \c void
		 */
		inline void printPt(std::ostream &out, const ANNpoint& p, bool newline = true) const;

		/*!
		 * \brief reads a data or query point with its Cartesian coordinates
		 *
		 * \param[in] node a point with its coordinates
		 *
		 * \return \c void
		 */
		template <size_t dim>
		inline void readPt(const MathVector<dim>& node);

		/*!
		 * \brief reads  multiple data or query points with their Cartesian coordinates
		 *
		 * \param[in] nodes points with their coordinates
		 *
		 * \return \c void
		 */
		template <size_t dim>
		inline void readPts(const std::vector<MathVector<dim> >& nodes);

		/*!
		 * \brief rebuilds the kd-tree iff the timestep has changed
		 *
		 * \param[in] timestep which timestep is demanded
		 *
		 * \return \c void
		 */
		void rebuildTree(const T& timestep);

#ifdef MPMNEURON
		void rebuildTree();
#endif

		/*!
		 * \brief generates a hash for the current timestep
		 *
		 * \param[in] timestep which timestep is demanded
		 *
		 * \return \c number the hashcode (runtime: log(n))
		 */
		long genHash(const T& timestep);

		/*!
		 * \brief check if two numbers can be considered equal (with accuracy < eps)
		 *
		 * \param[in] a first number
		 * \param[in] b second number
		 *
		 * \return \c bool, true if considered equal
		*/
		bool areSame(number a, number b);
};




/*!
 * \brief sPoint represents a generalized n-dimensional point which has an associated membrane potential
 * \addtogroup mpm_plugin
 *
 */
class sPoint {

	friend std::ostream& operator<< (std::ostream& output, const sPoint& p);

	public:
		/*!
		 * \brief main constructor
		 *
		 * \param[in] coordinates coordinates of a hoc point from Neuron
		 * \param[in] Vm the associated membrane potential of the hoc point
		 * \param[in] dist the distance to the nearest neighbor (nearest grid point) in uG
		 * \param[in] index the index of the point w. r. t. the datafile in which the point is included
		 * \param[in] timestep specifies the timestep in which the hoc point should be evaluated
		 *
		 * \see uGPoint
		 * \see Vm2uG
		 */
		sPoint(const std::vector<number>& coordinates, number Vm, number dist, int index, long timestep);

		/*!
		 * \brief default constructor
		 */
		sPoint();

		/*!
		 * brief default destructor
		 */
		~sPoint();

		// getters
		number getVm() const;
		number getDist() const;
		number getIndex() const;
		long getTimestep() const;
		std::vector<number> getCoordinates() const;

	protected:
		int index;

		long timestep;

		number Vm;
		number dist;

		std::vector<number> coordinates; /// coordinates extracted from a .hoc timestep file
};




/*!
 * \brief uGPoint represents a class (composites sPoint) to store the k nearest neighbors (as sPoints) of a query point (UG grid point).
 * \addtogroup mpm_plugin
 */
class uGPoint {

	friend std::ostream& operator<< (std::ostream& output, const uGPoint& p);

	public:
		/*!
		 * \brief main constructor
		 *
		 * \param[in] coordinates coordinates of a query point (a grid point in uG)
		 * \param[in] nearestNeighbors the k nearest neighbors (sPoints) of the query point
		 *
		 * \see sPoint
		 * \see Vm2uG
		 */
		uGPoint(const std::vector<number>& coordinates, const std::vector<sPoint>& nearestNeighbors);

		/*!
		 * \brief default constructor
		 */
		uGPoint();

		/*!
		 * \brief default destructor
		 */
		~uGPoint();

		// getters
		std::vector<sPoint> getNearestNeighbors() const;
		std::vector<number> getCoordinates() const;
		number getVm();
		number getDist();

	protected:
		std::vector<number> coordinates; /// cartesian coordinates from an UG grid point
		std::vector<sPoint> nearestNeighbors;
};


} // end namespace mpm
} // end namespace ug


// include vm2ug implementation
#include "vm2ug_impl.h"

#endif // __H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG__
