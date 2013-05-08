/*!
 * \file vm2ug.h
 * \brief header for the membrane_potential_mapping
 * \addtogroup mpm_plugin
 * \see docs for additional information
 *
 * \author Stephan Grein
 * \date July, 2011
 *
 * Notes: Migrate to the UG internal kd-tree
 * TODO: Code optimization	and cleanup (cctor calls, typename vs. class for templates)
 * TOOD: introduce template parameter for DIM in mapper
 */

#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG__

#ifndef DIM
	#if defined(UG_DIM_3) and not defined(UG_DIM_2) and not defined(UG_DIM_1)
		#define DIM 3
	#elif defined(UG_DIM_2) and not defined(UG_DIM_3) and not defined (UG_DIM_1)
		#define DIM 2
	#elif defined(UG_DIM_1) and not defined(UG_DIM_2) and not defined (UG_DIM_3)
		#define DIM 1
	#else
		#define DIM 3
		#warning "Assuming DIM=3 for now, as provided DIM is ambiguous."
	#endif
#endif

#include <vector>
#include <string>
#include <istream>
#include <ostream>
#include <cmath>

#include "common/common.h"

#include <ANN/ANN.h>
#include "common_typedefs.h"
#include "mvec.h"

// begin namespace ug
namespace ug {
	// begin namespace mpm
	namespace membrane_potential_mapping {

		// forward declarations
		template <class T> class Vm2uG;
		template <class T> class sPoint;
		template <class T> class uGPoint;

		template <class T> std::ostream& operator<< (std::ostream& output, const Vm2uG<T>& p);
		template <class T> std::ostream& operator<< (std::ostream& output, const sPoint<T>& p);
		template <class T> std::ostream& operator<< (std::ostream& output, const uGPoint<T>& p);

		/*!
		 * \brief Vm2uG represents a class to perform the k nearest neighbor search with a kd tree and can return the neighbors encapsulates in classes
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
			   * \param[in] dim dimension of the datapoints in the input file (default: 3)
			   * \param[in] maxPts maximum number of datapoints in the input files (default: 10000)
			   * \param[in] eps approximation factor for nearest neighbor search (default: 0)
			   * \param[in] k search for k nearest neighbors (default: 1)
			   */
			  Vm2uG(const std::string& dataFileBaseName, short int dim, int maxPts, number eps, short int k);

			  /*!
			   * \brief default destructor
			   */
			  ~Vm2uG();

			  /*!
			   * \brief interpolate linearly the membrane potentials
			   *
			   * \param[in] node the query point (uG)
			   * \param[in] timestep for which the nearest neighbor should be computed
			   * \param[in] cutoff maximum distance for nearest neighbor
			   * \param[in] k how many nearest neighbors should be considered for linear interpolation
			   * \return \c number the membrane potential (potentially linearly interpolated)
			   */
			 number interp_lin_vms(const T& timestep, number node[], number cutoff, int k);

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
			 number interp_bilin_vms(const T& timestep, number node[], number cutoff, int k);

			 /*!
			  * \brief return the associated membrane potential of the first nearest neighbor at cartesian coordinates x,y,z at timestep t: special case DIM=3.
			  *
			  * \param[in] x coordinate x
			  * \param[in] y coordinate y
			  * \param[in] z coordinate z
			  * \param[in] t timestep
			  *
			  * \return \c number the membrane potential
			  */
			  number get_potential(number x, number y, number z, const std::string& t) {
				 number node[3] = {x, y, z};
				 return vm_t(t, node).getVm();
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
			  * \return \c number the membrane potential (linearly interpolated)
			  */
			  number get_potential_lin(number x, number y, number z, const std::string& t, number cutoff, int k=2) {
				  number node[3] = {x,y,z};
				  return interp_lin_vms(t, node, cutoff, k);
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
			  * \return \c number the membrane potential (bilinearly interpolated)
			  */
			  number get_potential_bilin(number x, number y, number z, const std::string& t, number cutoff, int k=4) {
				  number node[3] = {x, y, z};
				  return interp_bilin_vms(t, node, cutoff, k);
			  }

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
			   * \return \c uGPoint<T> an grid point from UG with additional attributes, e. g. membrane potential, distance to query point, etc.
			   */
			  uGPoint<T> vm_t(const T& timestep, number node[]);

			  /*!
			   * \brief the same as vm_t except supplied is a list of query points; parallelized with OMP.
			   *
			   * \param[in] timestep the timestep for which a mmebrane potential shall be searched
			   * \param[in] nodes the query points (cartesian coordinates for grid points from UG)
			   *
			   * \return \c std::vector<uGPoint<T> > grid points from UG with additional attributes, e. g. membrane potential, distance to query points, etc.
			   *
			   * \see QDIM
			   * \see vm_t
			   * \see sPoint
			   * \see uGPoint
			   *
			   */
			  std::vector<uGPoint<T> > vm_t_many(const T& timestep, number nodes[][DIM]);



			  /*!
			   * \brief the same as vm_t but get the k nearest neighbors of one query points
			   *
			   * \param[in] timestep the timestep for which a membrane potential shall be searched
			   * \param[in] node the query point (cartesian coordinate of grid point from UG)
			   * \param[in] k how many nearest neighbors should be searched
			   *
			   * \return \c uGPoint<T> an grid point from UG with additional attributes, e. g. membrane potentials, distance to query points, etc.
			   * \see vm_t
			   * \see uGPoint
			   * \see sPoint
			  */
			  uGPoint<T> vm_t_k(const T& timestep, number node[], int k);

			  /*!
			   * \brief the same as vm_t_many but get the k nearest neighbors of many query points
			   *
			   * \param[in] timestep the timestep for which a membrane potential shall be searched
			   * \param[in] nodes the query points (cartesian coordinates for grid points from UG)
			   * \param[in] k how many nearest neighbors should be searched
			   *
			   * \return \c std::vector<uGPoint<T> > grid points from UG with additional attributes, e. g. membrane potentials, distance to query points, etc.
			   *
			   * \see QDIM
			   * \see vm_t_many
			   * \see sPoint
			   * \see uGPoint
			   */
			  std::vector<uGPoint<T> > vm_t_many_k(const T& timestep, number nodes[][DIM], int k);

			  // setters
			  void setK(short int k);
			  void setDim(short int dim);
			  void setTimestep(const T& timestep);
			  void setMaxPts(int maxPts);
			  void setEps(number eps);
			  void setPromise(bool promise);
			  void setdataFileBaseName(const std::string& dataFileBaseName);
			  void setdataFileExt(const std::string& dataFileExt);

			  // getters
			  bool treeBuild() { return isTreeBuild; }

		   protected:
			  bool promise;
			  bool isTreeBuild;

			  short int k;
			  short int dim;

			  int maxPts;

			  long timestep;

			  number eps;

			  std::string dataFileBaseName;
			  std::string dataFileExt;

		   private:
			  short int qPts;
			  int nPts;

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
			   * \brief reads a data or query point with its cartesian coordinates
			   *
			   * \param[in] node a point with its coordinates
			   *
			   * \return \c void
			   */
			  inline void readPt(number node[]);

			  /*!
			   * \brief rebuilds the kd-tree iff the timestep has changed
			   *
			   * \param[in] timestep which timestep is demanded
			   *
			   * \return \c void
			   */
			  void rebuildTree(const T& timestep);

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
		template <class T> class sPoint {

		   friend std::ostream& operator<< <> (std::ostream& output, const sPoint<T>& p);

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
			  sPoint(const std::vector<number>& coordinates, number Vm, number dist, int index, long timestep, const T& realfilename);

			  /*!
			   * \brief default constructor
			   */
			  sPoint();

			  /*!
			   * brief default destructor
			   */
			  ~sPoint();

			  // getters
			  inline number getVm() const;
			  inline number getDist() const;
			  inline number getIndex() const;
			  inline long getTimestep() const;
			  inline std::vector<number> getCoordinates() const;

		   protected:
			  int index;

			  long timestep;

			  number Vm;
			  number dist;

			  T realfilename;

			  std::vector<number> coordinates; /// coordinates extracted from a .hoc timestep file
		};

		/*!
		 * \brief uGPoint represents a class (composites sPoint) to store the k nearest neighbors (as sPoints) of a query point (UG grid point).
		 * \addtogroup mpm_plugin
		 */
		template <class T> class uGPoint {

		   friend std::ostream& operator<< <> (std::ostream& output, const uGPoint& p);

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
			  uGPoint(const std::vector<number>& coordinates, const std::vector<sPoint<T> >& nearestNeighbors);

			  /*!
			   * \brief default constructor
			   */
			  uGPoint();

			  /*!
			   * \brief default destructor
			   */
			  ~uGPoint();

			  // getters
			  inline std::vector<sPoint<T> > getNearestNeighbors() const;
			  inline std::vector<number> getCoordinates() const;
			  number getVm();
			  number getDist();

		   protected:
			  std::vector<number> coordinates; /// cartesian coordinates from an UG grid point
			  std::vector<sPoint<T> > nearestNeighbors;

		};
	// end namespace mpm
	}
// end namespace ug
}

// include vm2ug implementation
#include "vm2ug_impl.h"

#endif // __H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG__
