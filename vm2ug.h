/* Header {{{ */
/* Includes & Start guard {{{ */
#ifndef _VM2UG_HH_
#define _VM2UG_HH_
#include <ANN/ANN.h>
#include <vector>
#include <string>
#include <istream>
#include <cmath>
#include "common_typedefs.h"
#include "mvec.h"

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
      Vm2uG(std::string dataFileBaseName=std::string("timesteps/"), std::string dataFileExt=std::string(".csv"), const bool promise_ = false); 

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
      
      const double interp_bilin_vms(const T& timestep, const double node[], const double cutoff, const int k) {

    	  if (k < 4)
    		  return vm_t(timestep, node).getVm();

    	  double Vm_intp = 0.0;
    	  /* TODO:
    	   * pseudo code:
    	   *
    	   * Get K nearest neighbors
    	   * -> for all pairwise different points (Q11, Q12, Q21, Q22) do
    	   * ---> create plane with points Q11, Q12, Q21, Q22 iff det(Q11-Q12, Q11-Q21, Q11-Q22) = 0! (use sarrus rule to calculate determinant)!
    	   * ---> project query point (node) onto plane (Lotfußpunkt), see mathenexus: AbstandPE_Hilfsger_Paraform.htm
    	   * ---> calculate distance, and set current_cutoff = distance
    	   * ---> get coordinates of Lotfußpunkt L (given by previous step!)
    	   * ---> project L onto line Q12-Q22 => yields R1 (see interp_lin_vms for that!)
    	   * ---> project L onto line Q11-Q21 => yields R2
    	   * ----> calculate f(R2) = dist(Q12,R2) / dist(Q12,Q22) * f(Q12) + dist(Q22, R2)/dist(Q12,Q22) * f(Q22)
    	   * ----> calculate f(R1) = analogue ...
    	   * ------> calculate dist(R1,R2) => yields R3
    	   * ------> calculate dist(L, R1) => yields L1
    	   * ------> calculate dist(L, R2) => yields L2
    	   * --------> calculate L1/R3 * f(R1) + L2/R3 * f(R2) --> yields Vm bilinearly interpolated at point L which is nearest to query node
    	   *
    	   * please note: due to the fact that the points (Q11,Q12,Q21,Q22) which span the plane lie in that plane,
    	   * 				we have a degenerated 2D case, in which we can use bilinear interpolation!
    	   */
    	  return Vm_intp;
      }

      const double interp_lin_vms(const T& timestep, const double node[], const double cutoff, const int k) {

    	  // need at least two nearest neighbors!
    	  if (k < 2)
    		  return vm_t(timestep, node).getVm();

    	  double Vm_intp = 0.0;

    	  uGPoint<T> nearest = vm_t(timestep, node);
    	  std::vector<sPoint<T> > nearestPoints = nearest.getNearestNeighbors();


    	  if (nearestPoints[0].getVm() > cutoff)
    		  nearestPoints = vm_t_k(timestep, node, k);
    	  else
    		  return nearestPoints[0].getVm();

       	  mvecd3 u;

    	  for (int i=0; i < 3; i++)
    		  u.push_back(node[i]);

    	  double current_cutoff = cutoff;
    	  double c_norm;
    	  double a_norm;
    	  double dist;
    	  double rhs;
    	  double lhs;
    	  double r_star;
    	  double d_minus;
    	  double d_plus;
    	  double d_sum;

    	  typedef typename std::vector<sPoint<T> >::const_iterator SPIT;
    	  /**
    	   * iterate over all pairwise different points (nearest neighbors) and calculate
    	   * closest point on line ("Lotpunkt") w. r. t. the query point (node).
    	   *
    	   * iterate over all to find best one.
    	   */
    	  for (SPIT it1 = nearestPoints.begin(); it1 < nearestPoints.end(); it1++) {
    		  for (SPIT it2 = nearestPoints.begin(); it2 < nearestPoints.end(); it2++) {
    			  /*
    			   * calculate distance to the "Lotpunkt" (nearest point on line G)
    			   * w. r. t. the query point (node)
    			   */
    			  if (it1 != it2) {
    				  sPoint<T> t1 = *it1;
    				  mvecd3 m1 = t1.getCoordinates();
    				  sPoint<T> t2 = *it2;
    				  mvecd3 m2 = t2.getCoordinates();
    				  mvecd3 a = m2 - m1;
    				  mvecd3 b = u - m1;
    				  mvecd3 c = a % b;

    				  c_norm = c.norm(EUCLIDEAN);
    				  a_norm = a.norm(EUCLIDEAN);
    				  dist = (c_norm/a_norm);

    				  /**
    				   * construct now the line through the "Ortsvektor" m1
    				   * and the "Richtungsvektor" a. g: m1 + R * a
    				   *
    				   * calculate the scalar lambda "R" which gives us the
    				   * the coordinates of the point ("Lotpunkt") on line G
    				   */
    				  if (areSame(dist, current_cutoff)) {

    					  std::cout << "Point found!" << std::endl;

    					  rhs = a * m1;
    					  lhs = a * a;
    					  r_star = -rhs / lhs;

    					  mvecd3 pointOnG;
    					  mvecd3 R = std::vector<double>(3, r_star);

    					  pointOnG = m1 + (R % a);

    					  mvecd3 d_minus_v = m1 - pointOnG;
    					  mvecd3 d_plus_v  = pointOnG - m2;
    					  mvecd3 d_sum_v = m1 - m2;

    					  d_minus = d_minus_v.norm(EUCLIDEAN);
    					  d_plus  = d_plus_v.norm(EUCLIDEAN);
    					  d_sum   = d_sum_v.norm(EUCLIDEAN);

    					  Vm_intp = 0.0;
    					  Vm_intp += (d_minus/d_sum) * (t1.getVm());
    					  Vm_intp += (d_plus/d_sum) * (t2.getVm());
    					  current_cutoff = d_sum;
    				  }
    			  }
    		  }
    	  }
    	  return Vm_intp;
      }

 /*     const double interp_lin_vms(const T& timestep, const double node[], const double cutoff, const int maxIter) {

    				std::vector<double> a = ((*it2) - (*it1));
    				std::vector<double> b = (u - (*it1));
     				std::vector<double> c = vector_product(a, b);
     				double c_norm = norm(c);
     				double a_norm = norm(a);
     				double dist = (c/a);
     				// if nearest NEURON point is more distant than the cutoff distance
     				if (dist < current_cutoff)
     					std::cout << "point found" << std::endl;
     					std::vector<double> temp = *it1;
     					std::vector<double> temp2 = *it2;
     					const double rhs = vector_dot(a, temp);
     					const double lhs = vector_dot(a, a);
     					const double r_star = -rhs / lhs;
     					std::vector<double> pointOnG = temp;
     					std::vector<double> R;
     					for (unsigned int i = 0; i < 3; i++)
     						R.push_back(r_star);

     					pointOnG += vector_product(R, a);

     					std::vector<double> d_minus_v = temp - pointOnG;
     					std::vector<double> d_plus_v = pointOnG - temp2;
     					std::vector<double> d_sum_v = temp - temp2;

     					const double d_minus = norm(d_minus_v);
     					const double d_plus = norm(d_plus_v);
     					const double d_sum = norm(d_sum_v);


     					 sPoint<T> t1 = (*it1);
     					 sPoint<T> t2 = (*it2);

     					Vm_intp = 0.0;
     					Vm_intp += (d_minus/d_sum) * (t1.getVm());
     					Vm_intp += (d_plus/d_sum) * (t2.getVm());
     					current_cutoff = d_sum;


    			 }

    		  }
    	  }
    	  return Vm_intp;
      } */

      double get_potential(double x, double y, double z, std::string t) {
         double  node[3] = {x,y,z};
         return vm_t(t, node).getVm();
      }

      std::string foo(std::string s) { return s; }
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
#endif /* _VM2UG_HH_ */
/* }}} */
/* }}} */
