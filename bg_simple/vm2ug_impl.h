/* Vm2uG {{{ */
template <class T> Vm2uG<T>::Vm2uG(string dataFileBaseName_, string dataFileExt_, const bool promise_) {
   dim = 3;
   timestep = 0;
   maxPts = 100000;
   k = 1;
   eps = 0.0;
   isTreeBuild = false;
   dataFileBaseName = dataFileBaseName;
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

   cout << timestep;

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
   cout << timestep << endl;
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

   cout << timestep;
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
   cout << timestep;

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

template <class T> DOUBLE uGPoint<T>::getVm() {
   return (DOUBLE)nearestNeighbors[0].getVm();
}

template <class T> DOUBLE uGPoint<T>::getDist() {
   return (DOUBLE)nearestNeighbors[0].getDist();
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
