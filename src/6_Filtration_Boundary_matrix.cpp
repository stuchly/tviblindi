#include <Rcpp.h>
#include <vector>
#include <map>
#include <algorithm>
using namespace Rcpp;

std::vector< std::vector<int> > subs (   // returns all relevant subsets of vector, gained by removing one vertex
    std::vector<int>  simplex) {
  std::vector<int> tmp;
  std::vector< std::vector<int> > output;
  for(int ii = 0; ii < simplex.size(); ++ii) {
    tmp = simplex;
    tmp.erase(tmp.begin()+ii);
    output.push_back(tmp);
  }
  return output;
}

// [[Rcpp::export]]
RcppExport SEXP build_boundary_matrix (   // provided a list of simplices, makes list of previous simplices which are boundaries
    Rcpp::List      filtration,
    bool            sort_input = true) {
  
  std::vector <std::vector<int> > ff;  
  
  if (sort_input) { // sort input
    for (int ii = 0; ii < filtration.size(); ++ii) {
      NumericVector vv;
      vv = filtration[ii];
      std::sort (vv.begin(), vv.end());
      filtration[ii] = vv;
    }
  }
  
  std::map< std::vector<int>, int > assoc;
  
  for (int ii = 0; ii < filtration.size(); ++ii)
    assoc.insert(std::pair< std::vector<int>, int > (filtration[ii], ii));
  
  for (int ii = 0; ii < filtration.size(); ++ii) {
    std::vector< std::vector<int> > ss = subs (filtration[ii]);
    std::vector<int> h1;
    
    for (int jj = 0; jj < ss.size(); ++jj) {
      std::map< std::vector<int>, int >::iterator it = assoc.find(ss[jj]);
      if (it != assoc.end())
        h1.push_back(it->second+1);
    }
    
    // if (h1.size() == 0) h1.push_back(-1);
    
    ff.push_back(h1);
  }
  
  for (int ii = 0; ii < ff.size(); ++ii) { // sort output
    std::vector<int> vv;
    vv = ff[ii];
    std::sort (vv.begin(), vv.end());
    ff[ii] = vv;
  }
  
  return wrap(ff);
}