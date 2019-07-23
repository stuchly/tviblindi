#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <queue>
#include <limits>
#include <math.h>
#include <unordered_map>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdarg>
#include <numeric>
#include <cfloat>


#include "RcppArmadillo.h"
#include "vptree.h"


using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
RcppExport SEXP createKNNgraph_cos(SEXP X,SEXP D, SEXP K, SEXP N){
  Rcpp::NumericVector _X(X);

  int _D = Rcpp::as<int >(D);
  int _K = Rcpp::as<int >(K);
  int _N = Rcpp::as<int >(N);
  Rcpp::NumericMatrix _IND(_N,_K);
  Rcpp::NumericMatrix _DIST(_N,_K);



  VpTree<DataPoint, cosine_distance>* tree = new VpTree<DataPoint, cosine_distance>();

  std::vector<DataPoint> obj_X(_N, DataPoint(_D, -1, _X.begin()));

  for(int n = 0; n < _N; n++) obj_X[n] = DataPoint(_D, n, _X.begin() + n * _D);
  tree->create(obj_X);

  printf("Building tree...\n");
  std::vector<DataPoint> indices;
  std::vector<double> distances;

  for(int n = 0; n < _N; n++) {

    if(n % 10000 == 0) printf(" - point %d of %d\n", n, _N);

    // Find nearest neighbors
    indices.clear();
    distances.clear();
    tree->search(obj_X[n], _K + 1, &indices, &distances);
    for (int i=0; i<_K; i++){
      _IND(n,i)=indices[i].index();
      _DIST(n,i)=distances[i];
    }
  }


  return Rcpp::List::create(Named("IND",   _IND),Named("DIST", _DIST));

}

// [[Rcpp::export]]
RcppExport SEXP createKNNgraph_eu(SEXP X,SEXP D, SEXP K, SEXP N){
  Rcpp::NumericVector _X(X);

  int _D = Rcpp::as<int >(D);
  int _K = Rcpp::as<int >(K);
  int _N = Rcpp::as<int >(N);
  Rcpp::NumericMatrix _IND(_N,_K);
  Rcpp::NumericMatrix _DIST(_N,_K);


  VpTree<DataPoint, euclidean_distance>* tree = new VpTree<DataPoint, euclidean_distance>();

  std::vector<DataPoint> obj_X(_N, DataPoint(_D, -1, _X.begin()));

  for(int n = 0; n < _N; n++) obj_X[n] = DataPoint(_D, n, _X.begin() + n * _D);
  tree->create(obj_X);

  printf("Building tree...\n");
  std::vector<DataPoint> indices;
  std::vector<double> distances;

  for(int n = 0; n < _N; n++) {

    if(n % 10000 == 0) printf(" - point %d of %d\n", n, _N);

    // Find nearest neighbors
    indices.clear();
    distances.clear();
    tree->search(obj_X[n], _K + 1, &indices, &distances);
    for (int i=0; i<_K; i++){
      _IND(n,i)=indices[i].index();
      _DIST(n,i)=distances[i];
    }
  }


  return Rcpp::List::create(Named("IND",   _IND),Named("DIST", _DIST));

}



// ######################################################################################
// assign values to 0- to n-simplices

std::string
  hash_func( // name for simplex
    std::vector<int> _in
  ) { 
    std::string _out = std::to_string(_in[0]);
    if (_in.size() == 1) return _out;
    for (int i = 1; i < _in.size(); ++i)
    {
      _out = _out + " ";
      _out = _out + std::to_string(_in[i]);
    }
    return _out;
  }

std::vector< std::vector<int> >
  subsets( // all (n-1)-element subsets from n-element set
    std::vector<int> simplex)
  {
    std::vector<int> tmp;
    std::vector< std::vector<int> > _out;
    for (int i = 0; i < simplex.size(); ++i)
    {
      tmp = simplex;
      tmp.erase(tmp.begin() + i);
      _out.push_back(tmp);
    }
    return _out;
  }

std::vector< std::vector<int> >
  boundaries( // for each (n)-simplex, which (n+1)-simplices is it a boundary for?
    List f,
    bool sort_input = true
  ) {
    // init output
    std::vector< std::vector<int> > _boundaries;
    std::vector<int> def;
    def.push_back(-1);
    for (int i = 0; i < f.size(); ++i)
    {
      _boundaries.push_back(def);
    }
    
    if (sort_input)
    {
      for (int i = 0; i < f.size(); ++i)
      {
        NumericVector v = f[i];
        std::sort(v.begin(), v.end());
        f[i] = v;
      }
    }
    
    std::unordered_map< std::string, std::vector<int> > subset_map;
    // keys are (n+1)-simplex subsets, values are indices of (n+1)-simplices
    for (int i = 0; i < f.size(); ++i)
    {
      std::vector<int> simplex = f[i];
      if (simplex.size() > 2)
      {
        std::vector<std::vector<int> > ss = subsets(simplex);
        for (int j = 0; j < ss.size(); ++j)
        {
          std::string key = hash_func(ss[j]);
          std::unordered_map< std::string, std::vector<int> >::iterator it = subset_map.find(key);
          if (it == subset_map.end())
          {
            std::vector<int> val;
            val.push_back(i);
            subset_map.insert(make_pair(key, val));
          }
          else
          {
            std::vector<int> val = it->second;
            val.push_back(i);
            it->second = val;
          }
        }
      }
    }
    
    std::unordered_map< int, std::vector<int> > assoc;
    // keys are n-simplex indices, values are correspondent (n+1)-simplex indices
    for (int i = 0; i < f.size(); ++i)
    {
      std::vector<int> simplex = f[i];
      if (simplex.size() > 1)
      {
        std::string key = hash_func(simplex);
        std::unordered_map< std::string, std::vector<int> >::iterator it = subset_map.find(key);
        if (it != subset_map.end())
        {
          assoc.insert(make_pair (i, it->second));
        }
      }
    }
    
    std::unordered_map<int, std::vector<int> >::iterator it;
    
    for (it = assoc.begin(); it != assoc.end(); ++it)
    {
      int index = it->first;
      std::vector<int> val = it->second;
      _boundaries[index] = val;
    }
    for (int i = 0; i < _boundaries.size(); ++i)
    {
      std::vector<int> v;
      v = _boundaries[i];
      std::sort(v.begin(), v.end());
      _boundaries[i] = v;
    }
    
    return _boundaries;
  }

// [[Rcpp::export]]
std::vector< std::vector<int> >
  faces( // for each (n+1)-simplex, which (n)-simplices are its faces?
    List f,
    bool sort_input = true
  ) {
    // init output
    std::vector< std::vector<int> > _faces;
    
    std::vector<int> def; // default: no faces ~ -1
    def.push_back(-1);
    for (int i = 0; i < f.size(); ++i)
      _faces.push_back(def);
    
    if (sort_input)
    { // sort input simplex indices
      for (int i = 0; i < f.size(); ++i)
      {
        NumericVector v = f[i];
        std::sort(v.begin(), v.end());
        f[i] = v;
      }
    }
    
    // create a map of simplex names (keys) to indices (values)
    std::unordered_map<std::string, int> F;
    for (int i = 0; i < f.size(); ++i)
      F.insert(make_pair(hash_func(f[i]), i));
    
    // create a map of simplex indices (keys) to indices of their faces (values)
    std::unordered_map<int, std::vector<int> > C;
    for (int i = 0; i < f.size(); ++i)
    {
      std::vector<int> simplex = f[i];
      if (simplex.size() > 2) // (only if n > 1)
      {
        std::vector<int> vals;
        
        std::vector< std::vector<int> > subs = subsets(simplex);
        for (int j = 0; j < subs.size(); ++j)
        {
          std::string key = hash_func(subs[j]);
          int new_index = F.at(key);
          vals.push_back(new_index);
        }
        
        C.insert(make_pair(i, vals));
      }
    }
    
    // transfer values in map into a vector
    for (auto it = C.begin(); it != C.end(); ++it)
    {
      int index = it->first;
      std::vector<int> vals = it->second;
      _faces[index] = vals;
    }
    
    // sort values within entries of the vector
    for (int i = 0; i < _faces.size(); ++i)
    {
      std::vector<int> v;
      v = _faces[i];
      sort(v.begin(), v.end());
      _faces[i] = v;
    }
    
    return _faces;
  }

class
  fVertex // a vertex object
  {
    int ind_;
    int dim_;
    std::vector<double> loc_;
    
  public:
    fVertex(int ind, NumericMatrix coord)
    { // construct object given entire coords matrix and a row index
      ind_ = ind;
      dim_ = coord.ncol();
      for (int i = 0; i < dim_; ++i)
      {
        double val = coord(ind-1, i); // re-indexing
        loc_.push_back(val);
      }
    }
    fVertex(int ind, std::vector<double> points)
    { // construct object given index and an STL coordinate vector
      ind_ = ind;
      dim_ = points.size();
      for (int i = 0; i < dim_; ++i)
        loc_.push_back(points[i]);
    }
    fVertex(int ind, arma::vec points)
    { // construct object given index and an Armadillo coordinate vector
      ind_ = ind;
      dim_ = points.size();
      for (int i = 0; i < dim_; ++i)
        loc_.push_back(points[i]);
    }
    std::vector<double> location()
    { // get STL coordinate vector
      std::vector<double> _out;
      for (int i = 0; i < dim_; ++i)
        _out.push_back(loc_[i]);
      return _out;
    }
    int dimension() const { return dim_; } // get vertex dimension
    int index() const { return ind_; } // get vertex index
    double location(int d) const { return loc_[d]; } // get single coordinate
  };


double
  eucl_dist( // get Euclidean distance between two vertices
    const fVertex & v1,
    const fVertex & v2
  ) {
    double dist = .0;
    for (int i = 0; i < v1.dimension(); ++i)
      dist += (v1.location(i)-v2.location(i))*(v1.location(i)-v2.location(i));
    return sqrt(dist);
  }

class
  fSimplex // a simplex object
  {
    std::string name_;
    int npt_; // number of vertices
    std::vector<fVertex> vert_; // vector of vertices which make up the simplex
  public:
    fSimplex(std::vector<fVertex> vert)
    { // construct object given a list of vertices
      std::vector<int> v;
      for (int i = 0; i < vert.size(); ++i)
        v.push_back(vert[i].index());
      name_ = hash_func(v);
      npt_ = vert.size();
      vert_ = vert;
    }
    fSimplex(std::vector<int> inds, NumericMatrix coord)
    { // construct object given entire coords matrix and row indices
      npt_ = inds.size();
      for (int i = 0; i < npt_; ++i)
        vert_.emplace_back(inds[i], coord);
      name_ = hash_func(inds);
    }
    fVertex circumcenter()
    { // calculate coordinates of the circumcenter of a simplex
      
      const int _dim = vert_[0].dimension(); // dimension of simplex
      arma::vec _x(_dim); // _x is the answer vector of a system of linear equations
      
      if (npt_ > 2) // if simplex is triangle, tetrahedron, etc...
      { // ... compute coordinates of circumcenter by solving a linear system Ax = b
        std::vector< std::vector<double> > L; // list of vertex coordinates
        for (int i = 0; i < npt_; ++i)
        {
          std::vector<double> l = vert_[i].location();
          L.push_back (l);
        }
        
        arma::mat A (npt_-1,_dim); // coefficient matrix
        for (int pt = 0; pt < npt_-1; ++pt)
        {
          for (int dim = 0; dim < _dim; ++dim)
          {
            A (pt, dim) = -2*L[0][dim] + 2*L[pt+1][dim];
          }
        }
        
        arma::vec b (npt_-1); // answer vector
        for (int pt = 0; pt < npt_-1; ++pt)
        {
          b (pt) = .0;
          for (int dim = 0; dim < _dim; ++dim)
          {
            b (pt) += L[pt+1][dim]*L[pt+1][dim] - L[0][dim]*L[0][dim];
          }
        }
        _x = arma::solve(A,b);
      } else {
        if (npt_ == 2)
        {
          for (int i = 0; i < _dim; ++i)
          {
            _x(i) = (vert_[0].location(i)/2)+(vert_[1].location(i)/2);
          }
        }
      }
      return fVertex(-1, _x);
    }
    double gabriel_value(fVertex center)
    { // calculate alpha-complex filtration value of simplex, given that it is Gabriel
      
      double circumradius = eucl_dist(center, vert_[0]);
      return circumradius*circumradius;
    }
    
    double diameter(fVertex center)
    { // calculate diameter ~ Rips filtration value
      
      double circumradius = eucl_dist(center, vert_[0]);
      return circumradius*2;
    }
    
    std::vector<int> indices()
    {
      std::vector<int> _i;
      for (int i = 0; i < npt_; ++i)
      {
        _i.push_back(vert_[i].index());
      }
      return _i;
    }
    
    int n() const { return npt_-1; }
    int npt() const { return npt_; }
    std::string name() const { return name_; }
    int index(int nn) const { return vert_[nn].index(); }
    fVertex vertex(int i) const { return vert_[i]; }
  };

// [[Rcpp::export]]
SEXP
  rips_complex_filtration_values( // get filtration values for given filtration$cmplx
    List f,
    NumericMatrix coord,
    double threshold
  ) { // follows algorithm outlined at http://gudhi.gforge.inria.fr/doc/latest/group__rips__complex.html
    std::cout << "Computing filtration values for Rips complex." << std::endl;
    
    std::vector<double> _vals;
    
    int max_pt = 0;
    for (int i = 0; i < f.size(); ++i)
    {
      _vals.push_back(0); // initialise filtration values as 0
      std::vector<int> simplex = f[i];
      int s = simplex.size();
      if (s > max_pt)
        max_pt = s; // find max number of vertices in simplex
    }
    
    std::cout << "Max dimension is " << max_pt-1 << "." << std::endl;
    
    // store fVertex and fSimplex objects in vectors
    std::vector<fSimplex> S;
    
    for (int i = 0; i < f.size(); ++i)
    {
      std::vector<int> s = f[i]; // single simplex
      S.emplace_back(s, coord);
    }
    
    std::cout << "Complex contains " << S.size() << " simplices." << std::endl;
    
    std::vector< std::vector<int> > c = faces(f);
    std::unordered_map<int,double> ref; // keys are simplex indices, values are filtration values
    
    int npt = 2;
    for (int i = 0; i < S.size(); ++i)
    {
      fSimplex s = S[i];
      if (s.npt() == npt)
      {
        fVertex center = s.circumcenter();
        double val = s.diameter(center);
        if (val <= threshold)
        {
          std::pair<int,double> ins = std::make_pair(i, val);
          ref.insert(ins);
        } else {
          std::pair<int, double> ins = std::make_pair(i, 0);
          ref.insert(ins);
        }
      }
    }
    
    for (npt = 3; npt <= max_pt; ++npt)
    {
      for (int i = 0; i < S.size(); ++i)
      {
        fSimplex s = S[i];
        if (s.npt() == npt)
        {
          std::vector<int> ff = c[i];
          if (ff[0] != -1)
          {
            double val = ref.at(ff[0]);
            for (int j = 1; j < ff.size(); ++j)
            {
              double tmp = ref.at(ff[j]);
              std::cout << "tmp: " << tmp << std::endl;
              if (tmp > val && val != 0)
                val = tmp;
              if (tmp == 0)
                val = 0;
            }
            std::pair<int, double> ins = std::make_pair(i, val);
            ref.insert(ins);
          }
        }
      }
    }
    
    std::unordered_map<int,double>::iterator it;
    for (it = ref.begin(); it != ref.end(); ++it)
    {
      int ind = it->first;
      double val = it->second;
      _vals[ind] = val;
    }
    
    return wrap(_vals);
  }


// [[Rcpp::export]]
SEXP
  alpha_complex_filtration_values( // get filtration values for given filtration$cmplx
    List f,
    NumericMatrix coord
  ) { // follows algorithm outlined at http://gudhi.gforge.inria.fr/doc/latest/group__alpha__complex.html
    std::cout << "Computing filtration values for alpha-complex." << std::endl;
    
    std::vector<double> _vals;
    
    int max_pt = 0;
    for (int i = 0; i < f.size(); ++i)
    {
      _vals.push_back(0); // initialise filtration values as 0
      std::vector<int> simplex = f[i];
      int s = simplex.size();
      if (s > max_pt)
        max_pt = s; // find max number of vertices in simplex
    }
    
    std::cout << "Max dimension is " << max_pt-1 << "." << std::endl;
    
    // store fVertex and fSimplex objects in vectors
    std::vector<fVertex> V;
    std::vector<fSimplex> S;
    
    for (int i = 0; i < f.size(); ++i)
    {
      std::vector<int> s = f[i]; // single simplex
      std::vector<fVertex> v; // vertices of single simplex
      for (int j = 0; j < s.size(); ++j)
      {
        V.emplace_back(s[j], coord);
        v.emplace_back(s[j], coord);
      }
      S.emplace_back(s, coord);
    }
    
    std::cout << "Complex contains " << S.size() << " simplices." << std::endl;
    
    std::vector< std::vector<int> > b = boundaries(f);
    
    // vantage-point tree for nearest neighbour look-ups:
    VpTree<fVertex, eucl_dist> tree;
    tree.create(V);
    
    std::unordered_map<int,double> ref; // keys are simplex indices, values are filtration values
    
    for (int i = S.size()-1; i > -1; --i)
    {
      fSimplex s = S[i];
      if (s.npt() == max_pt)
      {
        fVertex center = s.circumcenter();
        std::pair<int,double> ins = std::make_pair(i, s.gabriel_value (center));
        ref.insert(ins);
      }
    }
    
    for (int pt = max_pt-1; pt > 1; --pt)
    {
      for (int i = S.size()-1; i > -1; --i)
      {
        fSimplex s = S[i];
        if (s.npt() == pt)
        {
          fVertex center = s.circumcenter();
          double val = s.gabriel_value(center);
          
          // check if simplex is Gabriel...
          std::vector<fVertex> results;
          std::vector<double> distances;
          tree.search(center, 1, &results, &distances);
          
          std::vector<int> simplex_vert = s.indices();
          
          if (find(simplex_vert.begin(), simplex_vert.end(), results[0].index()) == simplex_vert.end()) // if not Gabriel
          {
            std::vector<int> coface_ind = b[i];
            
            if (coface_ind.size() == 1 && coface_ind[0] == -1)
            {} else
            {
              val = ref.at(coface_ind[0]);
              for (int j = 1; j < coface_ind.size(); ++j)
              {
                double alt = ref.at(coface_ind[j]);
                if (alt < val)
                {
                  val = alt;
                }
              }
            }
          }
          ref.insert (std::make_pair (i, val));
        }
      }
    }
    
    std::unordered_map<int,double>::iterator it;
    for (it = ref.begin(); it != ref.end(); ++it)
    {
      int ind = it->first;
      double val = it->second;
      _vals[ind] = val;
    }
    
    return wrap(_vals);
  }

// ######################################################################################
