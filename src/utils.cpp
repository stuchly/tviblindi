/*
  reindex_cmplx:
      replace all vertex indices in a complex according to a map of new indices
  lowerbound_rb:
      create a modified reduced boundary matrix: set all entries in rows >= threshold to 0; adjust nonzero_col, low and dim accordingly
  path_1simplex_indices:
      get representation of walks in terms of 1-simplex indices in cmplx (replacement for function index_path)
*/

#include <Rcpp.h>
#include <algorithm>
#include <string>
#include <unordered_map>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector< std::vector<int> >
  reindex_cmplx(
    std::vector< std::vector<int> > cmplx,
    std::vector<int> inds // map of new indices
  ) {
    for (size_t i = 0; i < cmplx.size(); ++i)
      for (size_t j = 0; j < cmplx[i].size(); ++j)
       cmplx[i][j] = inds[cmplx[i][j]-1];
    
    return cmplx;
  }

// [[Rcpp::export]]
List
  lowerbound_rb( // set all entries in rows # >= threshold to 0
    List reduced_boundary_matrix,
    int threshold
  ) {
    if (threshold < 1)
      return wrap("lowerbound_rb: invalid threshold");
    
    std::vector< std::vector<int> > boundary = reduced_boundary_matrix[0];
    std::vector<int> nzc = reduced_boundary_matrix[1];
    std::vector<int> low = reduced_boundary_matrix[2];
    std::vector<int> dim = reduced_boundary_matrix[3];
    
    std::vector<size_t> zero_columns;
    
    for (size_t i = 0; i < boundary.size(); ++i)
    {
      for (size_t j = 0; j < boundary[i].size(); ++j)
        if (boundary[i][j] <= threshold)
          boundary[i][j] = -1;
      boundary[i].erase(std::remove(boundary[i].begin(), boundary[i].end(), -1), boundary[i].end());
      
      if (boundary[i].size() == 0)
        zero_columns.push_back(i);
      else
        low[i] = *std::max_element(boundary[i].begin(), boundary[i].end());
    }
    
    for (int i = zero_columns.size()-1; i >= 0; --i)
    {
      size_t this_ind = zero_columns[i];
      
      boundary.erase(boundary.begin()+this_ind);
      nzc.erase(nzc.begin()+this_ind);
      low.erase(low.begin()+this_ind);
      dim.erase(dim.begin()+this_ind);
    }
      
      return Rcpp::List::create(Named("boundary", boundary), Named("nonzero_col", nzc), Named("low", low), Named("dim", dim));
  }

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
  path_1simplex_indices(
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
