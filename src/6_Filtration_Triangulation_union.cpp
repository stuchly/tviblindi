#include <Rcpp.h>
#include <unordered_map>
using namespace Rcpp;

std::string
get_hash( // get name of simplex
          std::vector<int> simplex
) {
  std::string _hash = std::to_string(simplex[0]);
  if (simplex.size() == 1) return _hash;
  for (size_t i = 1; i < simplex.size(); ++i)
  {
    _hash = _hash + " ";
    _hash = _hash + std::to_string(simplex[i]);
  }
  return(_hash);
}

// fix: sort by original indices

// [[Rcpp::export]]
std::vector<int>
  unique_simplex_idcs(
    std::vector<std::vector<int> > cmplx, // list of simplices
    bool sort_input
  ) {
    std::vector<int> _unif; // output
    
    if (sort_input)
    {
      for (size_t i = 0; i < cmplx.size(); ++i)
      {
        std::vector<int> v = cmplx[i];
        std::sort(v.begin(), v.end());
        cmplx[i] = v;
      }
    }
    
    // 'hashes': associative map where key = simplex name, value = index
    typedef std::unordered_map<std::string, size_t> hash_table;
    typedef std::unordered_map<std::string, size_t>::iterator hash_table_iterator;
    
    hash_table hashes;
    for (size_t i = 0; i < cmplx.size(); ++i)
    {
      std::string hash = get_hash(cmplx[i]);
      hash_table_iterator it = hashes.find(hash);
      if (it == hashes.end())
      { // if given simplex not present in the map already, add its information
        hashes.insert(std::make_pair(hash, i));
      }
    }
    
    for (hash_table_iterator it = hashes.begin(); it != hashes.end(); ++it)
      _unif.push_back(it->second+1); // re-indexing
    
    std::sort(_unif.begin(), _unif.end());
    
    return _unif;
  }