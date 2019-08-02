#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

namespace straight {

  void print_vector(std::vector<int> v){
    int N=v.size();
    for (int i=0;i<N;i++) std::cout << v[i]<<", ";
    std::cout <<"\n";
    
  }

  std::vector<int> sym_diff(std::vector<int> aa, std::vector<int> bb) {
    std::vector<int>::iterator it;
    std::vector<int> output(aa.size()+bb.size());
    it = std::set_symmetric_difference(aa.begin(), aa.end(), bb.begin(), bb.end(), output.begin());
    output.resize(it-output.begin());
    return output;
  }

  std::vector<int> get_rep_int(std::vector<int> cycle, Rcpp::List boundary, std::vector<int> Ilow) {
    IntegerVector low(boundary[2]);
    IntegerVector nzc(boundary[1]);

    Rcpp::List boundaryV = boundary[0];

    std::vector<int> repre;

    std::sort(cycle.begin(), cycle.end());

    int ll = cycle[cycle.size()-1];           // ll <- max(cycle)

    int rr;
 
    while(cycle.size() > 0) {                 // while(ll > -Inf) { ...while cycle is not empty
      rr = Ilow[ll];                            // rr <- Ilow[ll]
    
      // if(rr < 0)
      //   return wrap("not a cycle?");    // if(rr == 0)  stop("not a cycle?")
  
      cycle = sym_diff(boundaryV[rr], cycle);    // cycle <- symdiff(boundary[[rr]],cycle)
  
      repre.push_back(nzc[rr]);            // repre <- c(repre, nzc[rr]) ---bacha nzc[rr]!!!
      ll = cycle[cycle.size()-1];                   // ll <- max(cycle)
    }                                           // }
    return repre;                         // repre
  }
}

// // [[Rcpp::export]]
// SEXP get_rep_heurestic(std::vector<int> cycle, Rcpp::List R,Rcpp::List B, bool update=false) {
 
//   IntegerVector low(R[2]);
//   IntegerVector nzc(R[1]);
//   Rcpp::List boundaryV = R[0];
//   Rcpp::List BV=B[0];
  

//   std::vector<int> Ilow(BV.length(),-1); //1-simplex last (complex size safe)

//   for(int ii = 0; ii < low.length(); ii++) {
//     Ilow[low[ii]] = ii;                 // Ilow[i] ~ which non-zero column has entry 1 in given row
//   }
 

//   std::vector<int > repre;
//   std::vector<int > Vl;

//   std::sort(cycle.begin(), cycle.end());

//   int ll = cycle[cycle.size()-1];           // ll <- max(cycle)

//   int rr;

//   while(cycle.size() > 0) {                 // while(ll > -Inf) { ...while cycle is not empty
    
//     R_CheckUserInterrupt();
//     rr = Ilow[ll];                            // rr <- Ilow[ll]
   
//     if(rr < 0){
//       std::vector<int> Bv=BV[ll-1];
   
//       std::vector<int > V_ll=straight::get_rep_int(Bv,R,Ilow);
     
//       if (V_ll.back()!=ll) V_ll.push_back(ll); // is it safe assertion?

//       std::sort(V_ll.begin(), V_ll.end());
//       cycle = straight::sym_diff(V_ll, cycle);
//       if (Vl.size()>0){
//         Vl=straight::sym_diff(Vl,V_ll);
//           } else{
//         Vl=V_ll;
//       }
//       repre.push_back(-ll);
//       ll = cycle[cycle.size()-1];
//     }
//     else
//       {
//         cycle = straight::sym_diff(boundaryV[rr], cycle);
//         repre.push_back(nzc[rr]);            // repre <- c(repre, nzc[rr]) !!!nzc[rr]
//         ll = cycle[cycle.size()-1];                   // ll <- max(cycle)
//       }
//   }

//   if (Vl.size()>0){
//     //std::sort(Vl.begin(),Vl.end());
//     boundaryV.push_back(Vl);
//     rr=boundaryV.length()-1;
//     ll=Vl.back();
//     Ilow[ll]=rr;
//     nzc.push_back(-ll); 
//     low.push_back(ll);
   
//     if (update) { //is this safe?? memory leak?
//       R[0]=boundaryV;
//       R[2]=low;
//       R[1]=nzc;
//     }
//   }
  
//   return wrap(repre);                         // repre
// }

// [[Rcpp::export]]
SEXP get_rep_straight(std::vector<int> cycle, Rcpp::List R,Rcpp::List B, bool update=false) {
  
  IntegerVector low(R[2]);
  std::vector<int > nzc=R[1];
  IntegerVector dims(R[3]);
  std::vector<double > vals=R[4];
  Rcpp::List boundaryV = R[0];
  Rcpp::List BV=B[0];
  
  
  std::vector<int> Ilow(BV.length(),-1); //1-simplex last (complex size safe)
  // std::cout<<"hu1"<<std::endl;
  for(int ii = 0; ii < low.length(); ii++) {
    Ilow[low[ii]] = ii;                 // Ilow[i] ~ which non-zero column has entry 1 in given row
  }
  // std::cout<<"hu2"<<std::endl;
  
  std::vector<int > repre;
  std::vector<int > Vl;
  
  std::sort(cycle.begin(), cycle.end());
  // std::cout<<"hu2.1"<<std::endl;
  int ll = cycle.back();       // ll <- max(cycle)
  // std::cout<<"hu2.2"<<std::endl;
  int rr;
  
  while(cycle.size() > 0) {                 // while(ll > -Inf) { ...while cycle is not empty
    R_CheckUserInterrupt();
    // std::cout<<"hu3.1"<<std::endl;
    rr = Ilow[ll];                            // rr <- Ilow[ll]
    //std::cout<<rr<<std::endl;
    if(rr < 0){
      // std::cout<<"hurr"<<std::endl;
      cycle.pop_back();
      Vl.push_back(ll);
      repre.push_back(-ll);
      ll = cycle.back();
    }
    else
      {
        cycle = straight::sym_diff(boundaryV[rr], cycle);
        repre.push_back(nzc[rr]);            // repre <- c(repre, nzc[rr]) !!!nzc[rr]
        ll = cycle.back();                   // ll <- max(cycle)
      }
  }
  // std::cout<<"hu3"<<std::endl;
  
  if (Vl.size()>0){
    std::reverse(Vl.begin(),Vl.end());
    //std::sort(Vl.begin(),Vl.end());
    boundaryV.push_back(Vl);
    rr=boundaryV.length()-1;
    ll=Vl.back();
    Ilow[ll]=rr;
    if (vals.back()<std::numeric_limits<double>::infinity()) nzc.push_back(BV.length()+1); else  nzc.push_back(nzc.back()+1); 
    low.push_back(ll);
    dims.push_back(1);
    vals.push_back(std::numeric_limits<double>::infinity());
    
    if (update) { //is this safe?? memory leak?
      R[0]=boundaryV;
      R[2]=low;
      R[1]=nzc;
      R[3]=dims;
      R[4]=vals;
    }
  }
  
  return wrap(repre);                         // repre
}

// [[Rcpp::export]]
SEXP
get_rep_straight_modified(
                          std::vector<int> cycle,
                          Rcpp::List R,
                          bool update=false
                          ) {
  IntegerVector low(R[2]);
  IntegerVector nonzero_col(R[1]);
  Rcpp::List rb = R[0];

  std::vector<int> which_low(*std::max_element(low.begin(), low.end())+1, -1);
  for(int i = 0; i < low.length(); i++)
    which_low[low[i]] = i;

  std::vector<int> repre;
  std::vector<int> Vl;

  std::sort(cycle.begin(), cycle.end());

  int cycle_max = cycle.back();

  int cycle_which_max;

  while(cycle.size() > 0)
    {
      R_CheckUserInterrupt();
      cycle_which_max = which_low[cycle_max];
   
      if(cycle_which_max < 0)
        {
          cycle.pop_back();
          Vl.push_back(cycle_max);
          repre.push_back(-cycle_max);
          cycle_max = cycle.back();
        }
      else
        {
          cycle = straight::sym_diff(rb[cycle_which_max], cycle);
          repre.push_back(nonzero_col[cycle_which_max]);
          cycle_max = cycle.back();
        }
    }

  if (Vl.size() > 0)
    {
      std::reverse(Vl.begin(),Vl.end());
      rb.push_back(Vl);
      cycle_which_max = rb.length()-1;
      cycle_max=Vl.back();
      which_low[cycle_max] = cycle_which_max;
      nonzero_col.push_back(-cycle_max); 
      low.push_back(cycle_max);
   
      if (update)
        {
          R[0]=rb;
          R[2]=low;
          R[1]=nonzero_col;
        }
    }
  
  return wrap(repre);
}

// [[Rcpp::export]]
std::vector< std::vector<int> >
get_reps_straight_modified(
                           std::vector< std::vector<int> > cycles,
                           Rcpp::List                      R, // reduced boundary matrix
                           bool                            update = false
                           ) {
  IntegerVector low(R[2]);
  IntegerVector nonzero_col(R[1]);
  Rcpp::List    rb = R[0];
    
  std::vector<int> which_low(*std::max_element(low.begin(), low.end())+1, -1);
  for(int i = 0; i < low.length(); i++)
    which_low[low[i]] = i;
    
  std::vector< std::vector<int> > representations;
    
  for (size_t i = 0; i < cycles.size(); ++i)
    {
      std::vector<int> cycle = cycles[i];
      std::vector<int> repre;
      std::vector<int> Vl;
      
      std::sort(cycle.begin(), cycle.end());
      
      int cycle_max = cycle.back();
      
      int cycle_which_max;
      
      while(cycle.size() > 0)
        {
          R_CheckUserInterrupt();
          cycle_which_max = which_low[cycle_max];
        
          if(cycle_which_max < 0)
            {
              cycle.pop_back();
              Vl.push_back(cycle_max);
              repre.push_back(-cycle_max);
              cycle_max = cycle.back();
            }
          else
            {
              cycle = straight::sym_diff(rb[cycle_which_max], cycle);
              repre.push_back(nonzero_col[cycle_which_max]);
              cycle_max = cycle.back();
            }
          if (Vl.size() > 0)
            {
              std::reverse(Vl.begin(),Vl.end());
              rb.push_back(Vl);
              cycle_which_max = rb.length()-1;
              cycle_max=Vl.back();
              which_low[cycle_max] = cycle_which_max;
              nonzero_col.push_back(-cycle_max); 
              low.push_back(cycle_max);
            }
        }
      
      representations.push_back(repre);
      std::cout << ".";
    }
  std::cout << std::endl;
    
  if (update)
    {
      R[0]=rb;
      R[2]=low;
      R[1]=nonzero_col;
    }
  return representations;
}
