#include <Rcpp.h>
#include <algorithm>
#include <string>
#include <unordered_map>
using namespace Rcpp;

std::string
onesimplex_hash(
  int a,
  int b,
  bool sort = true
) {
  if (sort)
    if (a > b)
    {
      int swap = a;
      a = b;
      b = swap;
    }
  
  std::string _out = std::to_string(a);
  _out = _out + " ";
  _out = _out + std::to_string(b);
  return _out;
}

// [[Rcpp::export]]
std::vector< std::vector<size_t> >
  walk_1simplex_indices(
    std::vector< std::vector<int> > walks,
    std::vector< std::vector<int> > cmplx
  ) {
    std::vector< std::vector<size_t> > _out;
    
    std::unordered_map<std::string, size_t> assoc;
    
    std::vector<int>      v = walks[0];
    std::vector<int> starts = walks[1];
    
    for (size_t i = 0; i < cmplx.size(); ++i)
      if (cmplx[i].size() == 2) // add 1-simplices to a dictionary
    {
      std::string hash = onesimplex_hash(cmplx[i][0], cmplx[i][1]);
      if (assoc.find(hash) == assoc.end())
        assoc.insert(std::make_pair(hash, i));
    }
    
    std::vector<int>::iterator first, last;
    
    for (size_t i = 0; i < starts.size(); ++i)
    {
      first = v.begin() + starts[i] - 1;
      last = (i == starts.size()-1) ? v.begin() + v.size() : v.begin() + starts[i+1] - 1;
      
      std::vector<int> path(first, last);
      
      std::vector<size_t> repre;
      
      std::unordered_map<std::string, size_t>::iterator it;
      for (size_t j = 0; j < path.size()-1; ++j)
      {
        std::string hash = onesimplex_hash(path[j], path[j+1]);
        it = assoc.find(hash);
        repre.push_back((it == assoc.end()) ? -1 : it->second+1);
      }
      _out.push_back(repre);
    }
    
    return _out;
  }