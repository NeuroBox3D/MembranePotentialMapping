/*!
 * \file vm2ug.h
 * \brief header for the membrane_potential_mapping
 * \see docs for additional information
 *
 * \author Stephan Grein
 * \date July, 2011
 *
 * Notes: Migrate to the UG internal kd-tree
 * TODO: Code optimization	and cleanup (cctor calls, typename vs. class for templates)
 */

#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG__

// macros (Please note: Symbol MPMDIM could not be resolved, e. g. in eclipse, is not an error.)
#ifndef DIM
	#if defined(UG_DIM_1) and not defined(UG_DIM_1) and not defined(UG_DIM_3)
		#define DIM 1
	#endif

	#if defined(UG_DIM_2) and not defined(UG_DIM_1) and not defined(UG_DIM_3)
		#define DIM 2
	#endif

	#if defined(UG_DIM_3) and not defined(UG_DIM_1) and not defined(UG_DIM_2)
		#define DIM 3
	#else
		#define DIM 3
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
			   * \param[in] promise specify if ordering of the datapoints in the input files can change (if true, tree needs to be rebuild less often)
			   */
			  Vm2uG(std::string dataFileBaseName=std::string("timesteps/"), std::string dataFileExt=std::string(".csv"), const bool promise_ = false);

			  /*!
			   * \brief enhanced constructor
			   *
			   * \param[in] dataFileBaseName basename of the input files (timesteps of the Neuron simulation)
			   * \param[in] dataFileExt extension of the input files (i.e.: .csv)
			   * \param[in] dim dimension of the datapoints in the input file (default: 3)
			   * \param[in] maxPts maximum number of datapoints in the input files (default: 10000)
			   * \param[in] eps approximation factor for nearest neighbor search (default: 0)
			   * \param[in] k search for k nearest neighbors (default: 1)
			   */
			  Vm2uG(std::string dataFileBaseName, const short int& dim, const int& maxPts, const double& eps, const short int& k);

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
			   * \return \c double the membrane potential (potentially linearly interpolated)
			   */
			 const double interp_lin_vms(const T& timestep, const double node[], const double cutoff, const int k);

			 /*!
			  * \brief interpolate bilinearly the membrane potentials
			  *
			  * \param[in] node the query point (uG)
			  * \param[in] timestep for which the nearest neighbor should be computed
			  * \param[in] cutoff maximum distance for nearest neighbor
			  * \param[in] k how many nearest neighbors should be considered for bilinear interpolation
			  *
			  * \return \c double the membrane potential (potentially bilinearly interpolated)
			  */
			 const double interp_bilin_vms(const T& timestep, const double node[], const double cutoff, const int k);

			 /*!
			  * \brief return the associated membrane potential of the first nearest neighbor at cartesian coordinates x,y,z at timestep t: special case DIM=3.
			  *
			  * \param[in] x coordinate x
			  * \param[in] y coordinate y
			  * \param[in] z coordinate z
			  * \param[in] t timestep
			  *
			  * \return \c double the membrane potential
			  */
			  double get_potential(double x, double y, double z, std::string t) {
				 double  node[3] = {x,y,z};
				 return vm_t(t, node).getVm();
			  }

			 /*!
			  * \brief return the associated membrane potential of the first nearest neighbor at cartesian coordinates x,y,z at timestep t: special case DIM=3.
			  *
			  * \param[in] x coordinate x
			  * \param[in] y coordinate y
			  * \param[in] z coordinate z
			  * \param[in] t timestep
			  *
			  * \return \c double the membrane potential (linearly interpolated)
			  */
			  double get_potential_lin(double x, double y, double z, std::string t, const double cutoff, const int k=2) {
				  double node[3] = {x,y,z};
				  return interp_lin_vms(t, node, cutoff, k);
			  }

			  /*!
			  * \brief return the associated membrane potential of the first nearest neighbor at cartesian coordinates x,y,z at timestep t: special case DIM=3.
			  *
			  * \param[in] x coordinate x
			  * \param[in] y coordinate y
			  * \param[in] z coordinate z
			  * \param[in] t timestep
			  *
			  * \return \c double the membrane potential (bilinearly interpolated)
			  */
			  double get_potential_bilin(double x, double y, double z, std::string t, const double cutoff, const int k=4) {
				  double node[3] = {x,y,z};
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
			  uGPoint<T> vm_t(const T& timestep, const double node[]);

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
			  std::vector<uGPoint<T> > vm_t_many(const T& timestep, const double nodes[][DIM]);



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
			  uGPoint<T> vm_t_k(const T& timestep, const double node[], const int& k);

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
			  std::vector<uGPoint<T> > vm_t_many_k(const T& timestep, const double nodes[][DIM], const int& k);

			  // setters
			  void setK(const short int& k);
			  void setDim(const short int& dim);
			  void setTimestep(const T& timestep);
			  void setMaxPts(const int& maxPts);
			  void setEps(const double& eps);
			  void setPromise(const bool& promise);
			  void setdataFileBaseName(std::string dataFileBaseName);
			  void setdataFileExt(std::string dataFileExt);

			  // getters
			  bool treeBuild() { return isTreeBuild; }

		   protected:
			  bool promise;
			  bool isTreeBuild;

			  short int k;
			  short int dim;

			  int maxPts;

			  long timestep;

			  double eps;

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
			  inline void printPt(std::ostream &out, const ANNpoint& p, const bool newline = true) const;

			  /*!
			   * \brief reads a data or query point with its cartesian coordinates
			   *
			   * \param[in] node a point with its coordinates
			   *
			   * \return \c void
			   */
			  inline void readPt(const double node[]);

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
			   * \return \c double the hashcode (runtime: log(n))
			   */
			  long genHash(const T& timestep);

			  /*!
			   * \brief check if two doubles can be considered equal (with accuracy < eps)
			   *
			   * \param[in] a first double
			   * \param[in] second double
			   *
			   * \return \c bool, true if considered equal
			   */
			  bool areSame(const double& a, const double& b);
		};

		/*!
		 * \brief sPoint represents a generalized n-dimensional point which has an associated membrane potential
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
			  sPoint(const std::vector<double>& coordinates, const double& Vm, const double& dist, const int& index, const long& timestep, const T realfilename);

			  /*!
			   * \brief default constructor
			   */
			  sPoint();

			  /*!
			   * brief default destructor
			   */
			  ~sPoint();

			  // getters
			  const inline double getVm() const;
			  const inline double getDist() const;
			  const inline double getIndex() const;
			  const inline long getTimestep() const;
			  const inline std::vector<double> getCoordinates() const;

		   protected:
			  int index;

			  long timestep;

			  double Vm;
			  double dist;

			  T realfilename;

			  std::vector<double> coordinates; /// coordinates extracted from a .hoc timestep file
		};

		/*!
		 * \brief uGPoint represents a class (composites sPoint) to store the k nearest neighbors (as sPoints) of a query point (UG grid point).
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
			  uGPoint(const std::vector<double>& coordinates, const std::vector<sPoint<T> >& nearestNeighbors);

			  /*!
			   * \brief default constructor
			   */
			  uGPoint();

			  /*!
			   * \brief default destructor
			   */
			  ~uGPoint();

			  // getters
			  const inline std::vector<sPoint<T> > getNearestNeighbors() const;
			  const inline std::vector<double> getCoordinates() const;
			  double getVm();
			  double getDist();

		   protected:
			  std::vector<double> coordinates; /// cartesian coordinates from an UG grid point
			  std::vector<sPoint<T> > nearestNeighbors;

		};
	// end namespace mpm
	}
// end namespace ug
}

// include vm2ug implementation
#include "vm2ug_impl.h"

#endif // __H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG__
