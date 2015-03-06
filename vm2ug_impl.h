/*!
 * \file vm2ug_impl.h
 * \brief implementation of vm2ug (\see vm2ug.h)
 *
 * \date Created on: Apr 27, 2012
 * \author Stephan Grein
 */

#ifndef __H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG_IMPL__
#define __H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG_IMPL__

// mpm includes
#include "vm2ug.h"
#include "mvec.h"
#include "common_typedefs.h"
#include "transformator.h"

namespace ug {
namespace membrane_potential_mapping {



template<class T>
template <size_t dim>
number Vm2uG<T>::interp_bilin_vms
(
	const T& timestep,
	const MathVector<dim>& coords,
	number cutoff,
	size_t k
)
{
	// need at least four nearest neighbors for bilinear interpolation
	if (k < 4)
		return vm_t(timestep, coords).getVm();

	number Vm_intp = 0.0;

	uGPoint nearest = vm_t(timestep, coords);
	std::vector<sPoint > nearestPoints = nearest.getNearestNeighbors();

	if (nearestPoints[0].getDist() > cutoff)
		nearestPoints = vm_t_k(timestep, coords, k).getNearestNeighbors();
		// nearestPoints = vm_t_k(timestep, node, k);
	else
		return nearestPoints[0].getVm();

	// TODO: This had better be reworked properly to work for dimensions 1 and 2 as well.
	// At the moment, for dim=1 and dim=2, we use zero values in the higher coordinates.
	// Does this produce correct results?
	mvecd3 u; // query point

	if (dim > 3)
		UG_THROW("Bilinear interpolation is only intended for dim=3 at most.")

	for (size_t i = 0; i < std::min((size_t) 3,dim); i++)
		u.push_back(coords[i]);

	if (dim < 3)
		for (size_t i = dim; i < 3; i++)
			u.push_back(0.0);

	typedef typename std::vector<sPoint>::const_iterator SPIT;

	number current_cutoff = cutoff;

	for (SPIT it1 = nearestPoints.begin(); it1 < nearestPoints.end(); it1++) {
		for (SPIT it2 = nearestPoints.begin(); it2 < nearestPoints.end(); it2++) {
			for (SPIT it3 = nearestPoints.begin(); it3 < nearestPoints.end(); it3++) {
				for (SPIT it4 = nearestPoints.begin(); it4 < nearestPoints.end(); it4++) {
					std::vector<number> coords1 = (*it1).getCoordinates();
					std::vector<number> coords2 = (*it2).getCoordinates();
					std::vector<number> coords3 = (*it3).getCoordinates();
					std::vector<number> coords4 = (*it4).getCoordinates();
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

						// constructs a "Hilfsgerade"
						mvecd3 p = mvecs[0];
						mvecd3 q = p - mvecs[1];
						mvecd3 r = p - mvecs[2];

						mvecd3 n = r % q;
						mvecd3 o = n - p;


						number rhs = -(p * n);
						number lhs = o * p;

						number sigma = rhs / lhs;

						mvecd3 mysigma;
						for (size_t t = 0; t < 3; t++)
							mysigma.push_back(sigma * n[t]);

						mvecd3 pointOnPlane = u + mysigma; // Lotfusspunkt on the plane, orthogonale Projektion

						number dist = (pointOnPlane - u).norm(EUCLIDEAN);

						if (dist < current_cutoff) {
							number local_rhs;
							number local_dist;
							number local_lhs;

							// project L onto line Q12-Q22 => yields R1
							// project L onto line Q11-Q21 => yields R2

							// point 1
							number c_norm = 0;
							number a_norm = 0;
							number r_star = 0;

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
							mvecd3 R = std::vector<number>(3, r_star);

							pointOnG = Q22 + (R % a);

							mvecd3 d_minus_v = Q22 - pointOnG;
							mvecd3 d_plus_v = pointOnG - Q12;
							mvecd3 d_sum_v = Q22 - Q12;

							/* d_minus = d_minus_v.norm(EUCLIDEAN);
							 *  d_plus = d_plus_v.norm(EUCLIDEAN);
							 *  d_sum = d_sum_v.norm(EUCLIDEAN);
							 */

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

							R = std::vector<number>(3, r_star);

							pointOnG = Q21 + (R % a);

							d_minus_v = Q21 - pointOnG;
							d_plus_v = pointOnG - Q11;
							d_sum_v = Q21 - Q11;

							/*  d_minus = d_minus_v.norm(EUCLIDEAN);
							 *  d_plus = d_plus_v.norm(EUCLIDEAN);
							 *  d_sum = d_sum_v.norm(EUCLIDEAN);
							 */

							mvecd3 R2 = pointOnG;

							// calculates bilinearly interpolated membrane potential
							number r1 = (Q12 - R2).norm(EUCLIDEAN)
									/ (Q12 - Q22).norm(EUCLIDEAN)
									* (*it1).getVm();
							r1 += (Q22 - R2).norm(EUCLIDEAN)
									/ (Q12 - Q22).norm(EUCLIDEAN)
									* (*it3).getVm();

							number r2 = (Q12 - R2).norm(EUCLIDEAN)
									/ (Q12 - Q22).norm(EUCLIDEAN)
									* (*it2).getVm();
							r2 += (Q22 - R2).norm(EUCLIDEAN)
									/ (Q12 - Q22).norm(EUCLIDEAN)
									* (*it4).getVm();

							number r3 = (R1 - R2).norm(EUCLIDEAN);
							number l1 = (pointOnPlane - R1).norm(EUCLIDEAN);
							number l2 = (pointOnPlane - R2).norm(EUCLIDEAN);

							Vm_intp = l1 / r3 * r1 + l2 / r3 * r2;
							current_cutoff = local_dist;
						}
					}
				}
			}
		}
	}

	return Vm_intp;
}



template<class T>
template <size_t dim>
number Vm2uG<T>::interp_lin_vms
(
	const T& timestep,
	const MathVector<dim>& coords,
	number cutoff,
	size_t k
)
{
	// need at least two nearest neighbors!
	if (k < 2)
		return vm_t(timestep, coords).getVm();

	number Vm_intp = 0.0;

	uGPoint nearest = vm_t(timestep, coords);
	std::vector<sPoint> nearestPoints = nearest.getNearestNeighbors();

	if (nearestPoints[0].getDist() > cutoff)
		nearestPoints = vm_t_k(timestep, coords, k).getNearestNeighbors();
	else
		return nearestPoints[0].getVm();

	// TODO: This had better be reworked properly to work for dimensions 1 and 2 as well.
	// At the moment, for dim=1 and dim=2, we use zero values in the higher coordinates.
	// Does this produce correct results?
	mvecd3 u; // query point

	if (dim > 3)
		UG_THROW("Linear interpolation is only intended for dim=3 at most.")

	for (size_t i = 0; i < std::min((size_t) 3,dim); i++)
		u.push_back(coords[i]);

	if (dim < 3)
		for (size_t i = dim; i < 3; i++)
			u.push_back(0.0);

	number current_cutoff = cutoff;
	number c_norm;
	number a_norm;
	number dist;
	number rhs;
	number lhs;
	number r_star;
	number d_minus;
	number d_plus;

	typedef typename std::vector<sPoint >::const_iterator SPIT;
	/*
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
				sPoint t1 = *it1;
				mvecd3 m1 = t1.getCoordinates();
				sPoint t2 = *it2;
				mvecd3 m2 = t2.getCoordinates();
				mvecd3 a = m2 - m1;
				mvecd3 b = u - m1;
				mvecd3 c = a % b;

				c_norm = c.norm(EUCLIDEAN);
				a_norm = a.norm(EUCLIDEAN);
				dist = (c_norm / a_norm);

				/*
				 * construct now the line through the "Ortsvektor" m1
				 * and the "Richtungsvektor" a. g: m1 + R * a
				 *
				 * calculate the scalar lambda "R" which gives us the
				 * the coordinates of the point ("Lotpunkt") on line G
				 */
				if (areSame(dist, current_cutoff)) {

					std::cout << "Point found!" << std::endl;

					rhs = a * b;
					lhs = a * a;
					r_star = -rhs / lhs;
					mvecd3 pointOnG = m1 + (a * r_star);

					mvecd3 d_minus_v = m1 - pointOnG;
					mvecd3 d_plus_v = pointOnG - m2;

					d_minus = d_minus_v.norm(EUCLIDEAN);
					d_plus = d_plus_v.norm(EUCLIDEAN);

					Vm_intp = 0.0;
					Vm_intp += (d_minus / a_norm) * (t1.getVm());
					Vm_intp += (d_plus / a_norm) * (t2.getVm());
					current_cutoff = dist;
				}
			}
		}
	}
	return Vm_intp;
}



template <class T> Vm2uG<T>::Vm2uG(const std::string& dataFileBaseName_, const std::string& dataFileExt_, bool promise_) {
   dim_file = 3;
   timestep = 0;
   maxPts = 100000;
   k = 1;
   eps = 0.0;
   isTreeBuilt = false;
   dataFileBaseName = dataFileBaseName_;
   dataFileExt = dataFileExt_;
   promise = promise_;
}


template <class T> Vm2uG<T>::Vm2uG (const std::string& dataFileBaseName, size_t dim_file, size_t maxPts, number eps, size_t k) {
   this->dataFileBaseName = dataFileBaseName;
   this->dim_file = dim_file;
   this->maxPts = maxPts;
   this->eps = eps;
}


template <class T> Vm2uG<T>::~Vm2uG () {
   if (this->isTreeBuilt) { /* if buildTree was called -> cleanup */
   delete [] nnIdx;
   delete [] dists;

   delete kdTree;
   delete queryPt; /* needs to be delete in any case after destruction */
   }
}


template <class T> void Vm2uG<T>::printPt(std::ostream &out, const ANNpoint& p, const bool newline) const {
   out << "(" << p[0];
   for (size_t i = 1; i < this->dim_file; i++) out << ", " << p[i];
   if (newline) out << ")" << std::endl;
}


template <class T>
template <size_t dim>
void Vm2uG<T>::readPt(const MathVector<dim>& node)
{
	size_t min_dim = std::min(dim, dim_file);
	for (size_t i = 0; i < min_dim; ++i) queryPt[i] = node[i];
	for (size_t i = dim; i < dim_file; ++i) queryPt[i] = 0.0;


#ifndef NDEBUG
	if (dim < dim_file)
	{
		UG_LOG("Warning: Query point has less coordinates (" << dim
			   << ") than specified in geometry file (" << dim_file << ").\n"
			   "Missing coordinates are assumed to be 0.");
	}
	else if (dim > dim_file)
	{
		UG_LOG("Warning: Query point has more coordinates (" << dim
			   << ") than specified in geometry file (" << dim_file << ").\n"
			   "Excess coordinates are ignored.");
	}
#endif
}


template <class T>
template <size_t dim>
void Vm2uG<T>::readPts(const std::vector<MathVector<dim> >& nodes)
{
	size_t min_dim = std::min(dim, dim_file);
	size_t sz = nodes.size();
	for (size_t j = 0; j < sz; ++j)
	{
		for (size_t i = 0; i < min_dim; ++i) queryPts[j][i] = nodes[j][i];
		for (size_t i = dim; i < dim_file; ++i) queryPts[j][i] = 0.0;
	}

#ifndef NDEBUG
	if (dim < dim_file)
	{
		UG_LOG("Warning: Query point has less coordinates (" << dim
		   << ") than specified in geometry file (" << dim_file << ").\n"
		   "Missing coordinates are assumed to be 0.");
	}
	else if (dim > dim_file)
	{
		UG_LOG("Warning: Query point has more coordinates (" << dim
			   << ") than specified in geometry file (" << dim_file << ").\n"
			   "Excess coordinates are ignored.");
	}
#endif
}


#ifdef MPMNEURON
template <class T> void Vm2uG<T>::buildTree() {
	if (this->isTreeBuilt) {
		 rebuildTree();
	} else {
		nPts = 0;
		queryPt = annAllocPt(dim_file);
		dim_file++;
		dataPts = annAllocPts(maxPts, dim_file);
		nnIdx = new ANNidx[k];
		dists = new ANNdist[k];

		typedef std::vector<std::pair<std::vector<number>, number> > VMS;
		VMS vms = (m_transformator->get_vms()).front();
		for (VMS::const_iterator cit = vms.begin(); cit != vms.end(); cit++) {
			int i = 0;
			for (std::vector<number>::const_iterator cit2 = (cit->first).begin(); cit2 != (cit->first).end(); cit2++) {
				dataPts[nPts][i] = *cit2;
				i++;
			}
			dataPts[nPts][i] = cit->second; // was cit->second
			nPts++;
		}

		nPts--;
		kdTree = new ANNkd_tree(dataPts, nPts, dim_file-1);
		dim_file--;
		this->isTreeBuilt = true;
	}
}
#endif


template <class T> void Vm2uG<T>::buildTree(const T& timestep) {

   if (this->isTreeBuilt) {
	  rebuildTree(timestep);
   }

   static std::ifstream dataStream;
   std::ostringstream s;
   s << timestep;

   dataStream.open((dataFileBaseName+s.str()+this->dataFileExt).c_str(), std::ios::in);
   if (!dataStream) {UG_THROW("Can not open data file \"" << dataFileBaseName+s.str()+this->dataFileExt << "\".");}
   else dataIn =& dataStream;

   nPts = 0;
   queryPt = annAllocPt(dim_file); // allocate initially queryPt (and deallocate it on destruction)
   dim_file++; // (n+1)th coordinate (Vm)
   dataPts = annAllocPts(maxPts, dim_file);
   nnIdx = new ANNidx[k];
   dists = new ANNdist[k];

   while (nPts < maxPts && !dataStream.eof()) {
	  for (size_t i = 0; i < dim_file; i++) *dataIn >> dataPts[nPts][i];
	  nPts++;
   }

   nPts--; // remove last (empty) point

   kdTree = new ANNkd_tree(dataPts, nPts, dim_file-1);
   dim_file--;

   dataStream.close();

   this->timestep = genHash(timestep);
   ts_oType = timestep;
   this->isTreeBuilt = true;
}


template <class T>
template <size_t dim>
std::vector<uGPoint> Vm2uG<T>::vm_t_many(const T& timestep, const std::vector<MathVector<dim> >& nodes) {
	return vm_t_many_k(timestep, nodes, k);
}


template <class T>
template <size_t dim>
uGPoint Vm2uG<T>::vm_t(const T& timestep, const MathVector<dim>& node) {
	return vm_t_k(timestep, node, k);
}


template <class T>
template <size_t dim>
uGPoint Vm2uG<T>::vm_t_k(const T& timestep, const MathVector<dim>& node, size_t k) {
	if (!this->isTreeBuilt) buildTree(timestep);
		if (genHash(timestep) == this->timestep) {
			readPt(node);
			this->kdTree->annkSearch(queryPt, k, this->nnIdx, this->dists, this->eps);
		} else {
			rebuildTree(timestep);
			vm_t_k(timestep, node, k);
		}

		std::vector<number> dataP;
		std::vector<number> queryP;
		std::vector<sPoint> sP;

		for (size_t i = 0; i < dim_file; i++) queryP.push_back(queryPt[i]);

		for (size_t j = 0; j < k; j++) {

		  for (size_t i = 0; i < dim_file; i++) dataP.push_back(dataPts[nnIdx[j]][i]);

		  sP.push_back(sPoint(dataP, sqrt(dists[j]), dataPts[nnIdx[j]][dim_file], nnIdx[j], genHash(timestep)));

		  dataP.clear();
		}

		return uGPoint(queryP, sP);
}


template <class T>
template <size_t dim>
std::vector<uGPoint> Vm2uG<T>::vm_t_many_k(const T& timestep, const std::vector<MathVector<dim> >& nodes, size_t k) {
	if (!this->isTreeBuilt) buildTree(timestep);

	annQueryPtsIdxArray.clear();
	annQueryPtsDistsArray.clear();

	qPts = nodes.size();

	queryPts = annAllocPts(qPts, dim_file);

	for (int i=0; i < qPts; i++) annQueryPtsIdxArray.push_back(new ANNidx[k]);
	for (int i=0; i < qPts; i++) annQueryPtsDistsArray.push_back(new ANNdist[k]);

	readPts(nodes);

	if (genHash(timestep) == this->timestep) {
		#ifdef _OPENMP
		#pragma omp parallel
		{
		#pragma omp for private(i)
		#endif
		for (int i=0; i < qPts; i++) {
			this->kdTree->annkSearch(queryPts[i], k, annQueryPtsIdxArray[i], this->annQueryPtsDistsArray[i], this->eps);
		}
		#ifdef _OPENMP
		}
		#endif
	} else {
		rebuildTree(timestep);
		#ifdef _OPENMP
		#pragma omp parallel
		{
		#pragma omp for private(i)
		#endif
		for (int i=0; i < qPts; i++) {
			this->kdTree->annkSearch(queryPts[i], k, annQueryPtsIdxArray[i], this->annQueryPtsDistsArray[i], this->eps);
		}
		#ifdef _OPENMP
		}
		#endif
	}

	std::vector<uGPoint > uGs;

	for (size_t i = 0; i < qPts; i++) {
		std::vector<number> dataP;
		std::vector<number> queryP;
		std::vector<sPoint> sP;

		for (size_t j = 0; j < dim_file; j++)
			queryP.push_back(queryPts[i][j]);

		for (size_t l = 0; l < k; l++) {
			for (size_t j = 0; j < dim_file; j++) dataP.push_back(dataPts[annQueryPtsIdxArray[i][l]][j]);
			sP.push_back(sPoint(dataP, sqrt(annQueryPtsDistsArray[i][l]), dataPts[annQueryPtsIdxArray[i][l]][dim_file], annQueryPtsIdxArray[i][l], genHash(timestep)));
			dataP.clear();
		}

		uGs.push_back(uGPoint(queryP,sP));

		queryP.clear();
		sP.clear();
	}

   delete[] queryPts; // cleanup queryPts (do not waste space/introduce mem leak)

   return uGs;
}


#ifdef MPMNEURON
template <class T> number Vm2uG<T>::vm_t(number node[]) {
	if (!this->isTreeBuilt) {
		buildTree();
	}

	readPt(node);

	this->kdTree->annkSearch(queryPt, this->k, this->nnIdx, this->dists, this->eps);

   return dataPts[nnIdx[0]][dim_file];
}
#endif


template <class T> void Vm2uG<T>::setK(size_t k) {
   this->k = k;
}

template <class T> void Vm2uG<T>::setDim(size_t dim_file) {
   this->dim_file = dim_file;
   this->promise = false;
   rebuildTree(ts_oType);
}

template <class T> void Vm2uG<T>::setTimestep(const T& timestep) {
   ts_oType = timestep;
   rebuildTree(timestep);
}

template <class T> void Vm2uG<T>::setMaxPts(size_t maxPts) {
   this->maxPts = maxPts;
   rebuildTree(ts_oType);
}

template <class T> void Vm2uG<T>::setEps(number eps) {
   this->eps = eps;
}

template <class T> void Vm2uG<T>::setPromise(bool promise) {
   this->promise = promise;
}

template <class T> void Vm2uG<T>::setdataFileBaseName(const std::string& dataFile) {
   this->dataFileBaseName = dataFile;
   rebuildTree(ts_oType);
}

template <class T> void Vm2uG<T>::setdataFileExt(const std::string& dataFileExt) {
   this->dataFileExt = dataFileExt;
   rebuildTree(ts_oType);
}


template <class T> bool Vm2uG<T>::areSame(number a, number b) {
   return std::fabs(a - b) < std::numeric_limits<number>::epsilon();
}

template <class T> long Vm2uG<T>::genHash(const T& timestep) {
	std::locale loc;
	const std::collate<char>& coll = std::use_facet<std::collate<char> >(loc);
	std::ostringstream s;
	s << timestep;
	std::string tmp = (this->dataFileBaseName+s.str()+this->dataFileExt).c_str();
	return coll.hash(tmp.data(),tmp.data()+tmp.length());
}

template <class T> void Vm2uG<T>::rebuildTree(const T& timestep) {
   //std::cout << "tree needs rebuild!" << std::endl;
   // tree needs to be rebuild if indices and or datapoints vary during timestepping (default: promise=false, i.e. we need to rebuild all)
   if (!promise) {
	  delete [] nnIdx;
	  delete [] dists;
	  delete kdTree; // note: deletes also dataPts memory (hopefully)
   }

   static std::ifstream dataStream;
   std::ostringstream s;
   s << timestep;
//   cout << timestep;

   dataStream.open((dataFileBaseName+s.str()+this->dataFileExt).c_str(), std::ios::in);
   if (!dataStream) {UG_THROW("Can not open data file \"" << dataFileBaseName+s.str()+this->dataFileExt << "\".");}
   else dataIn =& dataStream;


   nPts = 0;
   //queryPt = annAllocPt(dim_file); // let queryPt allocated, because queryPt dimension will not change
   dim_file++; // (n+1)th coordinate (Vm)
   if (!promise) {
	   dataPts = annAllocPts(maxPts, dim_file);
	   nnIdx = new ANNidx[k];
	   dists = new ANNdist[k];
   }

   while (nPts < maxPts && !dataStream.eof()) {
	  for (size_t i = 0; i < dim_file; i++) *dataIn >> dataPts[nPts][i];
	  nPts++;
   }
   nPts--; // remove last point

   // tree needs to be rebuild if indices and or datapoints vary during timestepping (default: promise=false, i.e. we need to rebuild)
   if (!promise) this->kdTree = new ANNkd_tree(dataPts, nPts, dim_file-1);

   dim_file--; // (n+1)th coordinate not needed anymore (was: Vm)

   dataStream.close();

   this->timestep=genHash(timestep);
   ts_oType = timestep;
}

#ifdef MPMNEURON
template <class T> void Vm2uG<T>::rebuildTree() {
	if (!promise) {
		delete [] nnIdx;
		delete [] dists;
		delete kdTree;
	}

	nPts = 0;
	dim_file++;
	if (!promise) {
		dataPts = annAllocPts(maxPts, dim_file);
		nnIdx = new ANNidx[k];
		dists = new ANNdist[k];
	}

	typedef std::vector<std::pair<std::vector<number>, number> > VMS;
	VMS vms = (m_transformator->get_vms()).front();
	for (VMS::const_iterator cit = vms.begin(); cit != vms.end(); cit++) {
		int i = 0;
		for (std::vector<number>::const_iterator cit2 = (cit->first).begin(); cit2 != (cit->first).end(); cit2++) {
			dataPts[nPts][i] = *cit2;
			i++;
		}

	dataPts[nPts][i+1] = cit->second;
	nPts++;
	}

	nPts--;
	if (!promise) this->kdTree = new ANNkd_tree(dataPts, nPts, dim_file-1);
	dim_file--;
}
#endif

template <class T> std::ostream& operator<<(std::ostream& output, const Vm2uG<T>& p) {
   output << "current instance of Vm2uG holds the properties:" << std::endl;
   output << "|" << "   dim_file: " << p.dim_file << " (dimension of points in file)" << std::endl;
   output << "|" << "          k: " << p.k << " (k nearest neihbors)" << std::endl;
   output << "|" << "        eps: " << p.eps << " (approximation factor)" << std::endl;
   output << "|" << "   BaseName: " << p.dataFileBaseName << " (current basename of files which hold timestep data)" << std::endl;
   output << "|" << "  Extension: " << p.dataFileExt << " (current extension for those files)" << std::endl;
   output << "|" << "   timestep: " << p.ts_oType << " (initial timestep for building the tree)" << std::endl;
   output << "|" << "     maxPts: " << p.maxPts << " (maximum number of points the tree can hold)" << std::endl;
   output << "|" << "    promise: " << (p.promise ? "true" : "false") << " (indicates if indices between dataFileBaseNames with timestep data do notvary)" << std::endl << std::endl;
   return output;
}


} // end namespace mpm
} // end namespace ug


#endif //__H__UG__MEMBRANE_POTENTIAL_MAPPING__VM2UG_IMPL__
