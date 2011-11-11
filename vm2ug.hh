/* Header {{{ */
/* Includes & Start guard {{{ */
#ifndef _VM2UG_HH_
#define _VM2UG_HH_
#include <ANN/ANN.h>
#include <vector>
#include <string>
#include <istream>

/* }}} */

/* Start namespace {{{ */
namespace vug {
/* }}} */

/* Macros {{{ */
#ifndef QDIM
#define QDIM 3
#define DOUBLE double
#endif
/* }}} */

/* Forward declarations {{{ */
template <class T> class Vm2uG;
template <class T> class sPoint;
template <class T> class uGPoint;

template <class T> std::ostream& operator<< (std::ostream& output, const Vm2uG<T>& p);
template <class T> std::ostream& operator<< (std::ostream& output, const sPoint<T>& p);
template <class T> std::ostream& operator<< (std::ostream& output, const uGPoint<T>& p);

template <class T = int> class Vm2uG;
template <class T = int> class sPoint;
template <class T = int> class uGPoint;
/* }}} */
/* }}} */

/* class Vm2uG {{{ */
/** Vm2uG represents a class to perform the k nearest neighbor search with a kd tree and can return the neighbors encapsulates in classes
  */
template <class T> class Vm2uG {

   friend std::ostream& operator<< <>(std::ostream& output, const Vm2uG<T>& p);

   public:
      /** main constructor 
      
      @param dataFileBaseName basename of the input files (timesteps of the Neuron simulation)
      @param dataFileExt extension of the input files (i.e. .csv)
      @param promise specify if ordering of the datapoints in the input files can change (if true, tree needs to be rebuild less often)
      
      @return Vm2uG
      */
      Vm2uG(std::string dataFileBaseName=std::string("timesteps/"), std::string dataFileExt=std::string(".csv"), const bool promise = false);

      /** enhanced constructor
   
      @param dataFileBaseName basename of the input files (timesteps of the Neuron simulation)
      @param dataFileExt extension of the input files (i.e.: .csv)
      @param dim dimension of the datapoints in the input file (default: 3)
      @param maxPts maximum number of datapoints in the input files (default: 10000)
      @param eps approximation factor for nearest neighbor search (default: 0)
      @param k search for k nearest neighbors (default: 1)
      
      @return Vm2uG
      */
      Vm2uG(std::string dataFileBaseName, const short int& dim, const int& maxPts, const double& eps, const short int& k);

      ~Vm2uG();
      

      double get_potential(double x, double y, double z, double t)	{return 0;}


      /** buildTree builds the tree with timestep specified in a class constructor
      
      @return void
      @see Vm2uG
      */
      void buildTree(const T& timestep);

      /** vm_t is a function for getting the nearest neighbor and returning its membrane potentials and distance to the query point
      
      @param node the query point (uG)
      @param timestep for which the nearest neighbor should be computed
      @return uGPoint with membrane potentials of the nearest neighbors (sPoint)

      @see uGPoint
      @see sPoint
      */
      uGPoint<T> vm_t(const T& timestep, const double node[]);
      
      /** vm_t_many is a function for getting the nearest neighbor of many query points w./w. o. parallelization (README)
      
      @param nodes the query points (uG)
      @param timestep for which the nearest neighbor should be computed
      @return std::vector of uGPoints with membrane potentials of the nearest neighbors (sPoint)

      @see QDIM
      @see vm_t

      @see uGPoint
      @see sPoint
      */
      
      std::vector<uGPoint<T> > vm_t_many(const T& timestep, const double nodes[][QDIM]); 

      /** vm_t_k is a function for getting the k nearest neighbors of many query points w./w. o. parallelization (README)
      
      @param node the query point
      @param timestep for which the nearest neighbors should be computed
      @return uGPoint with membrane potentials of the k nearest neighbors (sPoints) with decreasing distance

      @see uGPoint
      @see sPoint
      
      */
      uGPoint<T> vm_t_k(const T& timestep, const double node[], const int& k);

      /** same as vm_t_many but gets the k nearest neighbors of many query points
     
      @param nodes the query points
      @param timestep for which the nearest neighbors should be computed
      @param k how many nearest neighbors?

      @return std::vector of uGPoints with membrane potentials of the k nearest neighbors (sPoints) with decreasing distance for each point

      @SEE QDIM
      @see vm_t_many

      @see uGPoint
      @see sPoint

      */
      std::vector<uGPoint<T> > vm_t_many_k(const T& timestep, const double nodes[][QDIM], const int& k);
      
      void setK(const short int& k);
      
      void setDim(const short int& dim);

      void setTimestep(const T& timestep);
   
      void setMaxPts(const int& maxPts);

      void setEps(const double& eps);

      void setPromise(const bool& promise);

      void setdataFileBaseName(std::string dataFileBaseName);
      
      void setdataFileExt(std::string dataFileExt);
      
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
      
      /** prints a data or query point with its coordinates
      
      @param out where the output should go to (default: cout)
      @param p an ANNpoint p
      @param newline if true prints newline otherwise newline is omitted
      @return void
         
      */
      inline void printPt(std::ostream &out, const ANNpoint& p, const bool newline = true) const;

      /** reads a data or query point with its coordinates

      @param node a point with its coordinates
      @return void
      */
      inline void readPt(const double node[]);

      /** rebuilds the tree after timestep has changed
      
      @param timestep which timestep is demanded
      @return void
      */
      void rebuildTree(const T& timestep);
      
      /** generates a hash for the current timestep
      @param timestep which timestep is demanded
      @return the hashcode (logn)
      */
      long genHash(const T& timestep);
      
      /** check if two doubles can be considered equal (< eps)
      @param two doubles a and b
      @return logical value
      */
      bool areSame(const double& a, const double& b);

};
/* }}} */

/* class sPoint {{{ */
/** sPoint represents a generalized n-dimensional point which has an associated membrane potential Vm.
  */
template <class T> class sPoint {

   friend std::ostream& operator<< <> (std::ostream& output, const sPoint<T>& p);
   
   public:
      /** main constructor
      @param coordinates coordinates of a hoc point from Neuron
      @param Vm the associated membrane potential of the hoc point
      @param dist the distance to the nearest neighbor (nearest grid point) in uG
      @param index the index of the point w. r. t. the datafile in which the point is included
      @param timestep specifies the timestep in which the hoc point should be evaluated
      
      @see uGPoint
      @see Vm2uG
      */
      sPoint(const std::vector<double>& coordinates, const double& Vm, const double& dist, const int& index, const long& timestep, const T realfilename);
      sPoint();
      ~sPoint();

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

      std::vector<double> coordinates; /* {hoc} */

};
/* }}} */

/* class uGPoint {{{ */
/** uGPoint represents a class (composites sPoint) to store the k nearest neighbors (as sPoints) of a query point (uG grid point).
  */

template <class T> class uGPoint {
   
   friend std::ostream& operator<< <> (std::ostream& output, const uGPoint& p);

   public:

      /** main constructor
      @param coordinates coordinates of a query point (a grid point in uG)
      @param nearestNeighbors the k nearest neighbors (sPoints) of the query point

      @see sPoint
      @see Vm2uG
      */
      uGPoint(const std::vector<double>& coordinates, const std::vector<sPoint<T> >& nearestNeighbors);
      uGPoint();
      ~uGPoint();

      const inline std::vector<sPoint<T> > getNearestNeighbors() const;
      const inline std::vector<double> getCoordinates() const;

      /* returns uG DOUBLE */
      DOUBLE getVm();
      DOUBLE getDist();
      
   protected:
      std::vector<double> coordinates; /* {uG} */
      std::vector<sPoint<T> > nearestNeighbors; 

};
/* }}} */

/* Footer {{{ */
/* End namespace {{{ */
}
/* }}} */

/* End guard {{{ */
#endif
/* }}} */
/* }}} */
