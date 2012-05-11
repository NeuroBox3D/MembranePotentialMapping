
#include "vm2ug.h"

using namespace vug;

template<class T> const double Vm2uG<T>::interp_bilin_vms(const T& timestep,
		const double node[], const double cutoff, const int k) {
	// need at least four nearest neighbors for bilinear interpolation
	if (k < 4)
		return vm_t(timestep, node).getVm();

	double Vm_intp = 0.0;

	uGPoint<T> nearest = vm_t(timestep, node);
	std::vector<sPoint<T> > nearestPoints = nearest.getNearestNeighbors();

	if (nearestPoints[0].getVm() > cutoff)
		nearestPoints = vm_t_k(timestep, node, k).getNearestNeighbors();
		//nearestPoints = vm_t_k(timestep, node, k);
	else
		return nearestPoints[0].getVm();

	mvecd3 u; // query point

	for (int i = 0; i < 3; i++)
		u.push_back(node[i]);

	typedef typename std::vector<sPoint<T> >::const_iterator SPIT;

	double current_cutoff = cutoff;

	for (SPIT it1 = nearestPoints.begin(); it1 < nearestPoints.end(); it1++)
		for (SPIT it2 = nearestPoints.begin(); it2 < nearestPoints.end(); it2++)
			for (SPIT it3 = nearestPoints.begin(); it3 < nearestPoints.end();
					it3++)
				for (SPIT it4 = nearestPoints.begin();
						it4 < nearestPoints.end(); it4++) {
						std::vector<double> coords1 = (*it1).getCoordinates();
						std::vector<double> coords2 = (*it2).getCoordinates();
						std::vector<double> coords3 = (*it3).getCoordinates();
						std::vector<double> coords4 = (*it4).getCoordinates();
						mvecd3 Q11(coords1);
						mvecd3 Q12(coords2);
						mvecd3 Q21(coords3);
						mvecd3 Q22(coords4);

						mvecd3 c1 = Q11 - Q12;
						mvecd3 c2 = Q11 - Q21;
						mvecd3 c3 = Q11 - Q22;
						std::vector<mvecd3> mvecs;
						mvecs.push_back(c1);
						mvecs.push_back(c2);
						mvecs.push_back(c3);

						// iff all 4 points lay in one plane, then interpolate
						if (mvecd3::det(mvecs) != 0) {

							// construct a "Hilfsgerade"
							mvecd3 p = mvecs[0];
							mvecd3 q = p - mvecs[1];
							mvecd3 r = p - mvecs[2];

							mvecd3 n = r % q;
							mvecd3 o = n - p;


							double rhs = -(p * n);
							double lhs = o * p;

							double sigma = rhs / lhs;

							mvecd3 mysigma;
							for (size_t t = 0; t < 3; t++)
								mysigma.push_back(sigma * n[t]);

							mvecd3 pointOnPlane = u + mysigma; // Lotfußpunkt on the plane, orthogole Projektion

							double dist = (pointOnPlane - u).norm(EUCLIDEAN);

							if (dist < current_cutoff) {
								double local_rhs;
								double local_dist;
								double local_lhs;

								// project L onto line Q12-Q22 => yields R1
								// project L onto line Q11-Q21 => yields R2

								// point 1
								double c_norm = 0;
								double a_norm = 0;
								double r_star = 0;

								mvecd3 a = Q12 - Q22;
								mvecd3 b = pointOnPlane - Q22;
								mvecd3 c = a % b;

								c_norm = c.norm(EUCLIDEAN);
								a_norm = a.norm(EUCLIDEAN);
								local_dist = (c_norm / a_norm);
								local_rhs = a * Q22;
								local_lhs = a * a;
								r_star = -local_rhs / local_lhs;

								mvecd3 pointOnG;
								mvecd3 R = std::vector<double>(3, r_star);

								pointOnG = Q22 + (R % a);

								mvecd3 d_minus_v = Q22 - pointOnG;
								mvecd3 d_plus_v = pointOnG - Q12;
								mvecd3 d_sum_v = Q22 - Q12;

	/*							d_minus = d_minus_v.norm(EUCLIDEAN);
								d_plus = d_plus_v.norm(EUCLIDEAN);
								d_sum = d_sum_v.norm(EUCLIDEAN); */

								// point 2
								mvecd3 R1 = pointOnG;
								a = Q11 - Q21;
								b = pointOnPlane - Q21;
								c = a % b;

								c_norm = c.norm(EUCLIDEAN);
								a_norm = a.norm(EUCLIDEAN);
								local_dist = (c_norm / a_norm);
								local_rhs = a * Q21;
								local_lhs = a * a;
								r_star = -local_rhs / local_lhs;

								R = std::vector<double>(3, r_star);

								pointOnG = Q21 + (R % a);

								d_minus_v = Q21 - pointOnG;
								d_plus_v = pointOnG - Q11;
								d_sum_v = Q21 - Q11;

/*								d_minus = d_minus_v.norm(EUCLIDEAN);
								d_plus = d_plus_v.norm(EUCLIDEAN);
								d_sum = d_sum_v.norm(EUCLIDEAN); */

								mvecd3 R2 = pointOnG;

								// calculates bilinearly interpolated membrane potential
								double r1 = (Q12 - R2).norm(EUCLIDEAN)
										/ (Q12 - Q22).norm(EUCLIDEAN)
										* (*it1).getVm();
								r1 += (Q22 - R2).norm(EUCLIDEAN)
										/ (Q12 - Q22).norm(EUCLIDEAN)
										* (*it3).getVm();

								double r2 = (Q12 - R2).norm(EUCLIDEAN)
										/ (Q12 - Q22).norm(EUCLIDEAN)
										* (*it2).getVm();
								r2 += (Q22 - R2).norm(EUCLIDEAN)
										/ (Q12 - Q22).norm(EUCLIDEAN)
										* (*it4).getVm();

								double r3 = (R1 - R2).norm(EUCLIDEAN);
								double l1 = (pointOnPlane - R1).norm(EUCLIDEAN);
								double l2 = (pointOnPlane - R2).norm(EUCLIDEAN);

								Vm_intp = l1 / r3 * r1 + l2 / r3 * r2;
								current_cutoff = local_dist;
							}
						}
					}

	return Vm_intp;
}

template<class T> const double Vm2uG<T>::interp_lin_vms(const T& timestep,
		const double node[], const double cutoff, const int k) {

	// need at least two nearest neighbors!
	if (k < 2)
		return vm_t(timestep, node).getVm();

	double Vm_intp = 0.0;

	uGPoint<T> nearest = vm_t(timestep, node);
	std::vector<sPoint<T> > nearestPoints = nearest.getNearestNeighbors();

	if (nearestPoints[0].getVm() > cutoff)
		nearestPoints = vm_t_k(timestep, node, k).getNearestNeighbors();
	else
		return nearestPoints[0].getVm();

	mvecd3 u;

	for (int i = 0; i < 3; i++)
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
		for (SPIT it2 = nearestPoints.begin(); it2 < nearestPoints.end();
				it2++) {
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
				dist = (c_norm / a_norm);

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
					mvecd3 d_plus_v = pointOnG - m2;
					mvecd3 d_sum_v = m1 - m2;

					d_minus = d_minus_v.norm(EUCLIDEAN);
					d_plus = d_plus_v.norm(EUCLIDEAN);
					d_sum = d_sum_v.norm(EUCLIDEAN);

					Vm_intp = 0.0;
					Vm_intp += (d_minus / d_sum) * (t1.getVm());
					Vm_intp += (d_plus / d_sum) * (t2.getVm());
					current_cutoff = dist;
				}
			}
		}
	}
	return Vm_intp;
}

/* Vm2uG {{{ */
template <class T> Vm2uG<T>::Vm2uG(string dataFileBaseName_, string dataFileExt_, const bool promise_) {
   dim = 3;
   timestep = 0;
   maxPts = 100000;
   k = 1;
   eps = 0.0;
   isTreeBuild = false;
   this->dataFileBaseName = dataFileBaseName;
   dataFileExt = dataFileExt_;
   promise = promise_;
}

template <class T> Vm2uG<T>::Vm2uG (string dataFileBaseName, const short int& dim, const int& maxPts, const double& eps, const short int& k) {
   this->dataFileBaseName = dataFileBaseName;
   this->dim = dim;
   this->maxPts = maxPts;
   this->eps = eps;
}

template <class T> Vm2uG<T>::~Vm2uG () {
   if (this->isTreeBuild) { // if buildTree was called -> cleanup 
   delete [] nnIdx;
   delete [] dists;

   delete kdTree;
   delete queryPt; // needs to be delete in any case after destruction
   }
}

template <class T> void Vm2uG<T>::printPt(ostream &out, const ANNpoint& p, const bool newline) const {
   out << "(" << p[0];
   for (int i = 1; i < this->dim; i++) out << ", " << p[i];
   if (newline) out << ")" << endl;
}

template <class T> void Vm2uG<T>::readPt(const double node[]) {
   for (int i=0; i < dim; i++) queryPt[i] = node[i];
}

template <> void Vm2uG<double>::buildTree(const double& timestep) {
   if (this->isTreeBuild) {
      rebuildTree(timestep);
   }

   static ifstream dataStream;

   ostringstream s;
   
   s << timestep;

  // cout << timestep;

   dataStream.open((dataFileBaseName+s.str()+this->dataFileExt).c_str(), ios::in);
   if (!dataStream) {
      cerr << "Cannot open data file" << endl;
      exit(1);
   } else {
      dataIn =& dataStream;
   }

   nPts = 0;
   queryPt = annAllocPt(dim); // allocate initially queryPt (and deallocate it on destruction)
   dim++; // (n+1)th coordinate (Vm)
   dataPts = annAllocPts(maxPts, dim);
   nnIdx = new ANNidx[k];
   dists = new ANNdist[k];
   
   while (nPts < maxPts && !dataStream.eof()) {
      for (int i=0; i < dim; i++) *dataIn >> dataPts[nPts][i];
      nPts++;
   }

   nPts--; // remove last (empty) point
   
   kdTree = new ANNkd_tree(dataPts, nPts, dim-1);
   dim--;
   
   dataStream.close();

   this->timestep = genHash(timestep);
   this->isTreeBuild = true;

}

template <> void Vm2uG<string>::buildTree(const string& timestep) {
   //cout << timestep << endl;
   if (this->isTreeBuild) {
      rebuildTree(timestep);
   }

   static ifstream dataStream;

   dataStream.open((dataFileBaseName+timestep+this->dataFileExt).c_str(), ios::in);
   if (!dataStream) {
      cerr << "Cannot open data file" << dataFileBaseName+timestep+this->dataFileExt << endl;
      exit(1);
   } else {
      dataIn =& dataStream;
   }

   nPts = 0;
   queryPt = annAllocPt(dim); // allocate initially queryPt (and deallocate it on destruction)
   dim++; // (n+1)th coordinate (Vm)
   dataPts = annAllocPts(maxPts, dim);
   nnIdx = new ANNidx[k];
   dists = new ANNdist[k];
   
   while (nPts < maxPts && !dataStream.eof()) {
      for (int i=0; i < dim; i++) *dataIn >> dataPts[nPts][i];
      nPts++;
   }

   nPts--; // remove last (empty) point
   
   kdTree = new ANNkd_tree(dataPts, nPts, dim-1);
   dim--;
   
   dataStream.close();

   this->timestep = genHash(timestep);
   this->isTreeBuild = true;
}

template <class T> void Vm2uG<T>::buildTree(const T& timestep) {

   if (this->isTreeBuild) {
      rebuildTree(timestep);
   }

   static ifstream dataStream;
   ostringstream s;
   s << timestep; 

   //cout << timestep;
   dataStream.open((dataFileBaseName+s.str()+this->dataFileExt).c_str(), ios::in);
   if (!dataStream) {
      cerr << "Cannot open data file" << endl;
      exit(1);
   } else {
      dataIn =& dataStream;
   }
   
   nPts = 0;
   queryPt = annAllocPt(dim); // allocate initially queryPt (and deallocate it on destruction)
   dim++; // (n+1)th coordinate (Vm)
   dataPts = annAllocPts(maxPts, dim);
   nnIdx = new ANNidx[k];
   dists = new ANNdist[k];
   
   while (nPts < maxPts && !dataStream.eof()) {
      for (int i=0; i < dim; i++) *dataIn >> dataPts[nPts][i];
      nPts++;
   }

   nPts--; // remove last (empty) point
   
   kdTree = new ANNkd_tree(dataPts, nPts, dim-1);
   dim--;
   
   dataStream.close();

   this->timestep = genHash(timestep);
   this->isTreeBuild = true;
}



template <class T> std::vector<uGPoint<T> > Vm2uG<T>::vm_t_many(const T& timestep, const double nodes[][QDIM]) {

   annQueryPtsIdxArray.clear();
   annQueryPtsDistsArray.clear();
   
   const int size = sizeof(nodes[0])/sizeof(nodes[0][0]); // number of nodes 
   
   queryPts = annAllocPts(maxPts, dim);

   qPts = 0;

   for (int i=0; i < size; i++) annQueryPtsIdxArray.push_back(new ANNidx[k]);
   for (int i=0; i < size; i++) annQueryPtsDistsArray.push_back(new ANNdist[k]);

   while (qPts < size) {
      for (int i=0; i < dim; i++) queryPts[qPts][i] = nodes[qPts][i];
      qPts++;
   }

   if (timestep == this->timestep) {
         for (int i=0; i < size; i++) {
            this->kdTree->annkSearch(queryPts[i], this->k, annQueryPtsIdxArray[i], this->annQueryPtsDistsArray[i], this->eps);
         }  
   } else {
      this->timestep = timestep;
      rebuildTree(timestep);
         for (int i=0; i < size; i++) {
            this->kdTree->annkSearch(queryPts[i], this->k, annQueryPtsIdxArray[i], this->annQueryPtsDistsArray[i], this->eps);
      }
   }
   
   std::vector<uGPoint<T> > uGs;
   
   for (int i=0; i < size; i++) { 

      std::vector<double> dataP;
      std::vector<double> queryP;
      std::vector<sPoint<T> > sP;

      int j;
      int l;
      for (j=0; j<dim; j++) queryP.push_back(queryPts[i][j]);
      
      for (l=0; l<k;l++) {

         for (j=0; j<dim; j++) dataP.push_back(dataPts[annQueryPtsIdxArray[i][l]][j]);
      
         sP.push_back(sPoint<T>(dataP, sqrt(annQueryPtsIdxArray[i][l]), dataPts[annQueryPtsIdxArray[i][l]][dim], annQueryPtsIdxArray[i][l], genHash(timestep), timestep));

         dataP.clear();
      }

      uGs.push_back(uGPoint<T>(queryP,sP));
      
      queryP.clear();
      sP.clear();
   }
   
   delete[] queryPts; // cleanup queryPts (dont waste space)

   return uGs;
}

template <> uGPoint<string> Vm2uG<string>::vm_t(const string& timestep, const double node[]) {
   if (!this->isTreeBuild) buildTree(timestep);
   if (genHash(timestep) == this->timestep) { // timestep
   
   readPt(node);

   this->kdTree->annkSearch(queryPt, this->k, this->nnIdx, this->dists, this->eps);

   } else {
      this->timestep = genHash(timestep); //timestep
      rebuildTree(timestep);
      vm_t(timestep, node);
   }

   std::vector<double> dataP;
   std::vector<double> queryP;
   std::vector<sPoint<string> > sP;

   int i;
   int j;

   for (i=0; i<dim; i++) queryP.push_back(queryPt[i]);

   for (j=0; j < k; j++) {

      for (i=0; i<dim; i++) dataP.push_back(dataPts[nnIdx[j]][i]);

      sP.push_back(sPoint<string>(dataP, sqrt(dists[j]), dataPts[nnIdx[j]][dim], nnIdx[j], genHash(timestep), timestep)); //timestep und string
      
      dataP.clear();
   }

   return uGPoint<string>(queryP, sP);

}

template <class T> uGPoint<T> Vm2uG<T>::vm_t(const T& timestep, const double node[]) { // BOUNDARY NODE from uG
if (!this->isTreeBuild) buildTree(timestep);
   if (genHash(timestep) == this->timestep) { // timestep
   
   readPt(node);

   this->kdTree->annkSearch(queryPt, this->k, this->nnIdx, this->dists, this->eps);

   } else {
      this->timestep = genHash(timestep); //timestep
      rebuildTree(timestep);
      vm_t(timestep, node);
   }

   std::vector<double> dataP;
   std::vector<double> queryP;
   std::vector<sPoint<T> > sP;

   int i;
   int j;

   for (i=0; i<dim; i++) queryP.push_back(queryPt[i]);

   for (j=0; j < k; j++) {

      for (i=0; i<dim; i++) dataP.push_back(dataPts[nnIdx[j]][i]);

      sP.push_back(sPoint<T>(dataP, sqrt(dists[j]), dataPts[nnIdx[j]][dim], nnIdx[j], genHash(timestep), timestep)); //timestep
      
      dataP.clear();
   }

   return uGPoint<T>(queryP, sP);
}

template <class T> uGPoint<T> Vm2uG<T>::vm_t_k(const T& timestep, const double node[], const int& k) {
   if (this->k < k) {
      this->k = k;
      rebuildTree(timestep); // timestep, was: k
   }
   return vm_t(timestep, node);
}

template <class T> std::vector<uGPoint<T> > Vm2uG<T>::vm_t_many_k(const T& timestep, const double nodes[][QDIM], const int& k) {
   if (this->k < k) {
      this->k = k; 
      rebuildTree(timestep); // timestep, was: k
   }
   return vm_t_many(timestep, nodes);
}

template <class T> void Vm2uG<T>::setK(const short int& k) {
   this->k = k;
}
      
template <class T> void Vm2uG<T>::setDim(const short int& dim) {
   this->dim = dim;
   this->promise = false;
   rebuildTree(this->timestep);
}

template <class T> void Vm2uG<T>::setTimestep(const T& timestep) {
   this->timestep = timestep;
   rebuildTree(this->timestep);
}
  
template <class T> void Vm2uG<T>::setMaxPts(const int& maxPts) {
   this->maxPts = maxPts;
   rebuildTree(this->timestep);
}

template <class T> void Vm2uG<T>::setEps(const double& eps) {
   this->eps = eps;
}

template <class T> void Vm2uG<T>::setPromise(const bool& promise) {
   this->promise = promise;
}

template <class T> void Vm2uG<T>::setdataFileBaseName(string dataFile) {
   this->dataFileBaseName = dataFile;
   rebuildTree(this->timestep);
}

template <class T> void Vm2uG<T>::setdataFileExt(string dataFileExt) {
   this->dataFileExt = dataFileExt;
   rebuildTree(this->timestep);
}


template <class T> bool Vm2uG<T>::areSame(const double& a, const double& b) {
   return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}

template <class T> long Vm2uG<T>::genHash(const T& timestep) {
   locale loc;
   const collate<char>& coll = use_facet<collate<char> >(loc);
   ostringstream s;
   s << timestep; 
   string tmp = (this->dataFileBaseName+s.str()+this->dataFileExt).c_str();
   return coll.hash(tmp.data(),tmp.data()+tmp.length());
}

template <class T> void Vm2uG<T>::rebuildTree(const T& timestep) {
      cout << "tree needs rebuild!" << endl;
   // tree needs to be rebuild if indices and or datapoints vary during timestepping (default: promise=false, i.e. we need to rebuild all)
   if (!promise) {
      delete [] nnIdx;
      delete [] dists;
      delete kdTree; // note: deletes also dataPts memory (hopefully)
   }

   static ifstream dataStream;
   ostringstream s;
   s << timestep;
//   cout << timestep;

   dataStream.open((dataFileBaseName+s.str()+this->dataFileExt).c_str(), ios::in);
   if (!dataStream) {
      cerr << "Cannot open data file" << endl;
      exit(1);
   } else {
      dataIn =& dataStream;
   }
   
   nPts = 0;
   //queryPt = annAllocPt(dim); // let queryPt allocated, because queryPt dimension will not change
   dim++; // (n+1)th coordinate (Vm)
   if (!promise) {
   dataPts = annAllocPts(maxPts, dim);
   nnIdx = new ANNidx[k];
   dists = new ANNdist[k];
   }
   
   while (nPts < maxPts && !dataStream.eof()) {
      for (int i=0; i < dim; i++) *dataIn >> dataPts[nPts][i];
      nPts++;
   }
   nPts--; // remove last point
   
   // tree needs to be rebuild if indices and or datapoints vary during timestepping (default: promise=false, i.e. we need to rebuild)
   if (!promise) this->kdTree = new ANNkd_tree(dataPts, nPts, dim-1);
   
   dim--; // (n+1)th coordinate not needed anymore (was: Vm)

   dataStream.close(); 

   this->timestep=genHash(timestep);
}

template <class T> ostream& operator<<(ostream& output, const Vm2uG<T>& p) {
   output << "current instance of Vm2uG holds the properties:" << endl;
   output << "|" << "        dim: " << p.dim << " (dimension of query points)" << endl;
   output << "|" << "          k: " << p.k << " (k nearest neihbors)" << endl;
   output << "|" << "        eps: " << p.eps << " (approximation factor)" << endl;
   output << "|" << "   BaseName: " << p.dataFileBaseName << " (current basename of files which hold timestep data)" << endl;
   output << "|" << "  Extension: " << p.dataFileExt << " (current extension for those files)" << endl;
   output << "|" << "   timestep: " << p.timestep << " (initial timestep for building the tree)" << endl;
   output << "|" << "     maxPts: " << p.maxPts << " (maximum number of points the tree can hold)" << endl;
   output << "|" << "    promise: " << (p.promise ? "true" : "false") << " (indicates if indices between dataFileBaseNames with timestep data do notvary)" << endl << endl;
   return output;
}
/* }}} */

/* sPoint {{{ */

template <class T> sPoint<T>::sPoint(const std::vector<double>& coordinates, const double& dist, const double& Vm, const int& index, const long& timestep, const T realfilename) {
   if (coordinates.size() != 0) {
      try {
         this->coordinates = coordinates;
         this->dist = dist;
         this->Vm = Vm;
         this->index = index;
         this->timestep = timestep; // was: timestep
         this->realfilename = realfilename;
      } catch (const char* const str) {
      cout << "Initialization failed: " << str << endl;
      }
   }
   else {
      cout << "You created an EMPTY point (no coords). Check?!" << endl;
   }
}

template <class T> sPoint<T>::sPoint() { }

template <class T> sPoint<T>::~sPoint() { }

template <class T> const double sPoint<T>::getVm() const {
   return this->Vm;
}

template <class T> const double sPoint<T>::getDist() const {
   return this->dist;
}

template <class T> const double sPoint<T>::getIndex() const {
   return this->index;
}

template <class T> const long sPoint<T>::getTimestep() const {
   return this->timestep;
}

template <class T> const std::vector<double> sPoint<T>::getCoordinates() const {
   return this->coordinates;
}

template <class T> ostream& operator<<(ostream& output, const sPoint<T>& p) {
    output << "{hoc} (";
    unsigned int i;
    for (i=0; i < p.coordinates.size()-1; i++) output << p.coordinates[i] << ",";
    output << p.coordinates[i] << ")" << endl;
    output << "\tVm: " << p.Vm << " (mV)" << endl;
    output << "\tDist to query Point {uG}: " << p.dist << endl;
    output << "\tIndex {in timestep (Hashcode=" << p.timestep << ") file}: " << p.index << endl;
    output << "\tFilename (short): " << p.realfilename << endl;
    return output;  
}

/* }}} */

/* uGPoint {{{ */
template <class T> uGPoint<T>::uGPoint(const std::vector<double>& coordinates, const std::vector<sPoint<T> >& nearestNeighbors) {
   if (coordinates.size() != 0) {
      try {
         this->coordinates = coordinates;
         this->nearestNeighbors = nearestNeighbors;
      } catch (const char* const str) {
      cout << "Initialization failed: " << str << endl;
      }
   }
   else {
      cout << "You created an EMPTY point (no coords). Check?!" << endl;
   }
}

template <class T> uGPoint<T>::uGPoint() { }

template <class T> uGPoint<T>::~uGPoint() { }

template <class T> double uGPoint<T>::getVm() {
   return nearestNeighbors[0].getVm();
}

template <class T> double uGPoint<T>::getDist() {
   return nearestNeighbors[0].getDist();
}


template <class T> const std::vector<double> uGPoint<T>::getCoordinates() const {
   return this->coordinates;
}

template <class T> const std::vector<sPoint<T> > uGPoint<T>::getNearestNeighbors() const {
   return this->nearestNeighbors;
}

template <class T> ostream& operator<<(ostream& output, const uGPoint<T>& p) {
   output << "----------------------------------------------------------------------------------------" << endl;
   output << "query Point {uG} (";
   unsigned int i;
   for (i=0; i < p.coordinates.size()-1; i++) output << p.coordinates[i] << ", ";
   output << p.coordinates[i] << ")" << endl;
   output << "----------------------------------------------------------------------------------------" << endl;
   for (i=0; i < p.nearestNeighbors.size(); i++) output << "nearest neighbor #" << i << ": " << p.nearestNeighbors[i] << endl;
   output << "****************************************************************************************" << endl;
   return output;
   }
/* }}} */

