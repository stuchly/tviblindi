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

// [[Rcpp::export]]
std::vector< std::vector<int> >
  boundaries( // for each (n > 0)-simplex, which (n+1)-simplices is it a boundary for?
    List f, // list of complexes (as vectors of vertex indices)
    bool sort_input = true
  ) {
    // init output
    std::vector< std::vector<int> > _boundaries;
    std::vector<int> def; // default: no codimension-1 cofaces; indicated by -1
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
          assoc.insert(make_pair(i, it->second));
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
    List f, // list of complexes (as vectors of vertex indices)
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
  vec_dot( // dot product of vectors
    std::vector<double> a,
    std::vector<double> b
  ) {
    double _res = .0;
    int len = a.size();
    for (int i = 0; i < len; ++i)
    {
      _res += a[i] * b[i];
    }
    return _res;
  }

std::vector<double>
  vec_diff( // vector difference
    std::vector<double> a,
    std::vector<double> b
  ) {
    std::vector<double> _res;
    int len = a.size();
    for (int i = 0; i < len; ++i)
    {
      _res.push_back(a[i] - b[i]);
    }
    return _res;
  }

double
  eucl_dist( // get Euclidean distance between two vertices
    const fVertex & v1,
    const fVertex & v2
  ) {
    double dist = 0;
    for (int i = 0; i < v1.dimension(); ++i)
      dist += (v1.location(i)-v2.location(i))*(v1.location(i)-v2.location(i));
    return sqrt(dist);
  }

std::vector<double> simplex_circumcenter(
    std::vector< std::vector<double> > vert
)
{
  const int N = vert.size();
  const int dim = vert[0].size();
  
  arma::vec x(dim);
  
  std::vector< std::vector<double> > E; // E is a list of vectors from one vertex to all the others
  std::vector<double> origin = vert[N - 1];
  for (int i = 0; i < N - 1; ++i)
  {
    std::vector<double> target = vert[i];
    E.push_back(vec_diff(origin, target));
  }
  
  arma::mat A(N - 1, N - 1);
  for (int i = 0; i < N - 1; ++i)
  {
    for (int j = 0; j < N - 1; ++j)
    {
      A(i, j) = 2 * vec_dot(E[i], E[j]);
    }
  }
  
  arma::vec b(N - 1);
  for (int i = 0; i < N - 1; ++i)
  {
    b(i) = vec_dot(E[i], E[i]);
  }
  
  x = arma::solve(A, b); // we solve for circumcenter as barycentric coordinates...
                         // ...but the last coordinate is missing and must be filled in
  
  std::vector<double> w;
  double last = 1;
  for (int i = 0; i < x.size(); ++i)
  {
    w.push_back(x[i]);
    last -= x[i];
  }
  w.push_back(last);  // w is now the complete vector of coords of the circumcenter (in barycentric coordinate system)
  
  std::vector<double> res(dim);
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < dim; ++j)
    {
      res[j] += vert[i][j] * w[i]; // barycentric coordinates are used as weights to get coordinates in the original system
    }
  }
  return res;
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
    { // get circumcenter of simplex

      const int _dim = vert_[0].dimension();
      std::vector<double> _x(_dim);

      if (npt_ > 2)
      {
        std::vector< std::vector<double> > pts;
        for (int i = 0; i < npt_; ++i)
        {
          std::vector<double> pt;
          for (int j = 0; j < _dim; ++j)
          {
            pt.push_back(vert_[i].location(j));
          }
          pts.push_back(pt);
        }
        _x = simplex_circumcenter(pts);
      } else {
        if (npt_ == 2)
        {
          for (int i = 0; i < _dim; ++i)
          {
            _x[i] = (vert_[0].location(i)/2)+(vert_[1].location(i)/2);
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
  alpha_complex_filtration_values_C( // get filtration values for given filtration$cmplx
    List f,
    NumericMatrix coord,
    std::vector<int> f_u // unique vertex indices that appear in filtration
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
    
    // store fVertex and fSimplex objects from filtration in vectors S and V...
    std::vector<fVertex> V;
    std::vector<fSimplex> S;
    for (int i = 0; i < f.size(); ++i)
    {
      std::vector<int> s = f[i];
      S.emplace_back(s, coord);
    }
    
    V.emplace_back(1, coord); // filler
    for (int i = 0; i < f_u.size(); ++i)
    {
      V.emplace_back(f_u[i], coord);
    }
    
    std::cout << "Complex contains " << S.size() << " simplices." << std::endl;
    
    std::vector< std::vector<int> > b = boundaries(f);
    
    // create a vantage-point tree for nearest neighbour look-ups:
    VpTree<fVertex, eucl_dist> tree;
    tree.create(V);
    
    std::unordered_map<int,double> ref; // keys are simplex indices, values are filtration values
    
    for (int i = S.size()-1; i > -1; --i)
    {
      fSimplex s = S[i];
      if (s.npt() == max_pt) // if simplex is of max size, use its Gabriel value (no codimension-1 cofaces can exist, therefore it is Gabriel)
      {
        fVertex center = s.circumcenter();
        double val = s.gabriel_value(center);
        // std::cout << "Simplex " << i+1 << " filt val: " << val << std::endl;
        std::pair<int,double> ins = std::make_pair(i, val);
        ref.insert(ins);
      }
    }
    
    for (int pt = max_pt-1; pt > 1; --pt) // iterate over over simplices of given given dimension (dimension decreasing with each iteration)
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
          
          int res_idx = results[0].index();
          
          if (find(simplex_vert.begin(), simplex_vert.end(), res_idx) == simplex_vert.end())
          { // ...if s is not Gabriel
            
            std::vector<int> coface_ind = b[i];

            if (coface_ind.size() == 1 && coface_ind[0] == -1)
            {} else
                {
              bool found = false;

              double new_val = std::numeric_limits<double>::max();
              for (int j = 0; j < coface_ind.size(); ++j)
              {
                double alt = ref.at(coface_ind[j]);

                if (alt < new_val)
                {
                  found = true;
                  new_val = alt; // pick minimum
                }
              }
              if (found)
              {
                val = new_val;
              }
            }
          } else {
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
