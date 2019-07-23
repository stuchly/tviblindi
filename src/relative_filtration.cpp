#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::IntegerVector
  relative_filtration(
    std::vector<int> indices,
    std::vector< std::vector<int> > filtration,
    bool sort_output
  ) {
    std::vector<size_t> _out;
    
    typedef std::unordered_set<int>           index_map;
    typedef std::unordered_set<int>::iterator index_map_it;
    
    index_map    im;
    index_map_it it;
    
    for (size_t i = 0; i < indices.size(); ++i)
      im.insert(indices[i]);
    
    for (size_t i = 0; i < filtration.size(); ++i)
    {
      size_t len = filtration[i].size();
      for (size_t j = 0; j < len; ++j)
      {
        it = im.find(filtration[i][j]);
        if (it == im.end())
          break;
        if (j == len-1)
          _out.push_back(i+1); // re-indexing
      }
    }
    
    if (sort_output)
      std::sort(_out.begin(), _out.end());
    
    return(wrap(_out));
  }
