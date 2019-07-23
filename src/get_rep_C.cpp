#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;


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



// [[Rcpp::export]]
SEXP get_rep_C(std::vector<int> cycle, Rcpp::List boundary) {
  IntegerVector low(boundary[2]);
  IntegerVector nzc(boundary[1]);

  IntegerVector ss(nzc.length());
  for(int ii = 0; ii < nzc.length(); ii++)
    ss[ii] = ii;

  IntegerVector Ilow(*std::max_element(low.begin(), low.end())+1,-1);
  for(int ii = 0; ii < low.length(); ii++) {
    Ilow[low[ii]] = ss[ii];                 // Ilow[i] ~ which non-zero column has entry 1 in given row
  }

  Rcpp::List boundaryV = boundary[0];

  std::vector<int> repre;

  std::sort(cycle.begin(), cycle.end());

  int ll = cycle[cycle.size()-1];           // ll <- max(cycle)

  int rr;
  while(cycle.size() > 0) {                 // while(ll > -Inf) { ...while cycle is not empty
    rr = Ilow[ll];                            // rr <- Ilow[ll]
    // std::cout << rr<<"\n";
    if(rr < 0)
      return wrap("not a cycle?");    // if(rr == 0)  stop("not a cycle?")
    
    cycle = sym_diff(boundaryV[rr], cycle);    // cycle <- symdiff(boundary[[rr]],cycle)

    repre.push_back(nzc[rr]);            // repre <- c(repre, nzc[rr])
    ll = cycle[cycle.size()-1];                   // ll <- max(cycle)
  }                                           // }
  return wrap(repre);                         // repre
}


std::vector<int> get_rep_int(std::vector<int> cycle, Rcpp::List boundary, std::vector<int> Ilow) {
  IntegerVector low(boundary[2]);
  IntegerVector nzc(boundary[1]);

  // IntegerVector ss(nzc.length());
  // for(int ii = 0; ii < nzc.length(); ii++)
  //   ss[ii] = ii;

  // IntegerVector Ilow(*std::max_element(low.begin(), low.end())+1,-1);
  // for(int ii = 0; ii < low.length(); ii++) {
  //   Ilow[low[ii]] = ss[ii];                 // Ilow[i] ~ which non-zero column has entry 1 in given row
  // }

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

// [[Rcpp::export]]
SEXP get_rep_lazy(std::vector<int> cycle, Rcpp::List R,Rcpp::List B, bool update=false) {
  // Rcpp::Environment base = Environment("package:base");
  // Rcpp::Function readline = base["readline"];
  // Rcpp::Function as_numeric = base["as.numeric"];

  IntegerVector low(R[2]);
  IntegerVector nzc(R[1]);
  IntegerVector dim(R[3]);
  Rcpp::List boundaryV = R[0];
  Rcpp::List BV=B[0];
  
  IntegerVector ss(nzc.length());
  for(int ii = 0; ii < nzc.length(); ii++)
    ss[ii] = ii;

  // std::vector<int> Ilow(*std::max_element(cycle.begin(), cycle.end())+1,-1);
  std::vector<int> Ilow(BV.length(),-1); //1-simplex last (complex size safe)

  // std::cout<<"pred low \n";
  for(int ii = 0; ii < low.length(); ii++) {
    Ilow[low[ii]] = ss[ii];                 // Ilow[i] ~ which non-zero column has entry 1 in given row
  }
  // std::cout<<"bum\n";

  std::vector<int> repre;

  std::sort(cycle.begin(), cycle.end());

  int ll = cycle[cycle.size()-1];           // ll <- max(cycle)

  int rr;

  while(cycle.size() > 0) {                 // while(ll > -Inf) { ...while cycle is not empty
   
    // std::cout<<cycle.size()<<"\n";
    // int drink = Rcpp::as<int>(as_numeric(readline("> ")));
    // std::cout<<"uhehe: "<<cycle.back()<<"\n";
    R_CheckUserInterrupt();
    rr = Ilow[ll];                            // rr <- Ilow[ll]
   
    if(rr < 0){
      std::vector<int> Bv=BV[ll-1];
   
      std::vector<int > V_ll=get_rep_int(Bv,R,Ilow);
     
      if (V_ll.back()!=ll) V_ll.push_back(ll); // is it safe assertion?
      
      std::sort(V_ll.begin(), V_ll.end());
      cycle = sym_diff(V_ll, cycle);
      repre.push_back(-ll);
      // std::cout<<ll<<"\n";
      // print_vector(V_ll);
      
      //if (update){//this should be safe 
      boundaryV.push_back(V_ll);
      rr=boundaryV.length()-1;
      Ilow[ll]=rr;
      nzc.push_back(-ll); 
      low.push_back(ll);
      //}
          
      ll = cycle[cycle.size()-1];
         
    }
    else
      {
        //std::vector<int> L(1);
      
        cycle = sym_diff(boundaryV[rr], cycle);
        repre.push_back(nzc[rr]);            // repre <- c(repre, nzc[rr]) !!!nzc[rr]
        ll = cycle[cycle.size()-1];                   // ll <- max(cycle)
      }
  }                                           // }

  // std::cout<<"pred upd \n";
  if (update) { //is this safe?? memory leak?
    R[0]=boundaryV;
    R[2]=low;
    R[1]=nzc;
  }
  // std::cout<<"end \n";

 
  return wrap(repre);                         // repre
}


std::vector<int> get_rep_int_2(std::vector<int> cycle, Rcpp::List boundary, std::vector<int> Ilow) {
  IntegerVector low(boundary[2]);
  IntegerVector nzc(boundary[1]);


  Rcpp::List boundaryV = boundary[0];

  std::vector<int> repre;

  std::sort(cycle.begin(), cycle.end());

  int ll = cycle[cycle.size()-1];           // ll <- max(cycle)

  int rr;
 
  while(cycle.size() > 0) {                 // while(ll > -Inf) { ...while cycle is not empty
    rr = Ilow[ll];                            // rr <- Ilow[ll]
  
    cycle = sym_diff(boundaryV[rr], cycle);    // cycle <- symdiff(boundary[[rr]],cycle)
  
    repre.push_back(-nzc[rr]);            // repre <- c(repre, nzc[rr]) ---bacha nzc[rr]!!!
    ll = cycle[cycle.size()-1];                   // ll <- max(cycle)
  }                                           // }
  return repre;                         // repre
}


// [[Rcpp::export]]
SEXP get_rep_lazy_2(std::vector<int> cycle, Rcpp::List R,Rcpp::List B, bool update=false) {
  // Rcpp::Environment base = Environment("package:base");
  // Rcpp::Function readline = base["readline"];
  // Rcpp::Function as_numeric = base["as.numeric"];

  std::vector<std::vector<int > > VLL, CYCLE;
  std::vector<int > LL,II;
  
  IntegerVector low(R[2]);
  IntegerVector nzc(R[1]);
  IntegerVector dim(R[3]);
  Rcpp::List boundaryV = R[0];
  Rcpp::List BV=B[0];
  
  IntegerVector ss(nzc.length());
  for(int ii = 0; ii < nzc.length(); ii++)
    ss[ii] = ii;

  // std::vector<int> Ilow(*std::max_element(cycle.begin(), cycle.end())+1,-1);
  std::vector<int> Ilow(BV.length(),-1); //1-simplex last (complex size safe)

  std::cout<<"pred low \n";
  for(int ii = 0; ii < low.length(); ii++) {
    Ilow[low[ii]] = ss[ii];                 // Ilow[i] ~ which non-zero column has entry 1 in given row
  }
  std::cout<<"bum\n";

  std::vector<int> repre;

  std::sort(cycle.begin(), cycle.end());

  int ll = cycle[cycle.size()-1];           // ll <- max(cycle)

  int rr;

  while(cycle.size() > 0) {                 // while(ll > -Inf) { ...while cycle is not empty
    CYCLE.push_back(cycle);
    // std::cout<<cycle.size()<<"\n";
    // int drink = Rcpp::as<int>(as_numeric(readline("> ")));
    // std::cout<<"uhehe: "<<cycle.back()<<"\n";
    R_CheckUserInterrupt();
    rr = Ilow[ll];                            // rr <- Ilow[ll]
    LL.push_back(ll);
    if(rr < 0){
      II.push_back(-1);
      // std::cout<<ll<<"\n";
      // std::cout<<cycle.back()<<"\n";
      std::vector<int> Bv=BV[ll-1];
   
      std::vector<int > V_ll=get_rep_int_2(Bv,R,Ilow);
      // std::vector<int > V_lll;
      // for (int i=0;i<V_ll.size();i++){
      //   std::vector<int > Bvv=boundaryV[V_ll[i]];
      //   for (int j=0;j<Bvv.size();j++) V_lll.push_back(Bvv[j]);
          
      // }
     
      if (V_ll.back()!=ll) V_ll.push_back(ll); // is it safe assertion?
      
      // if (V_ll.back()==ll){
      //   std::cout<<"-------------\n";
      //   print_vector(cycle);
      //   print_vector(Bv);
      //   print_vector(V_ll);
        
      //   return(Rcpp::List::create(Rcpp::Named("V_ll",Rcpp::wrap(V_ll)),Rcpp::Named("cycle",Rcpp::wrap(cycle)),Rcpp::Named("Bv",Rcpp::wrap(Bv))));
      //   std::cout<<"au\n";
        
      // }
      
      std::sort(V_ll.begin(), V_ll.end());
      
      VLL.push_back(V_ll);
      cycle = sym_diff(V_ll, cycle);
      
      std::sort(cycle.begin(),cycle.end());
      if (cycle.back()<0){
        std::cout<<"humpf\n";
        boundaryV.push_back(cycle);
        break;
      }
      repre.push_back(-ll);
      
      boundaryV.push_back(V_ll);
      rr=boundaryV.length()-1;
      Ilow[ll]=rr;
      nzc.push_back(-ll); 
      low.push_back(ll);
       
          
      ll = cycle[cycle.size()-1];
         
    }
    else
      {
        //std::vector<int> L(1);
        II.push_back(1);
        VLL.push_back(boundaryV[rr]);
        cycle = sym_diff(boundaryV[rr], cycle);
        repre.push_back(nzc[rr]);            // repre <- c(repre, nzc[rr]) !!!nzc[rr]
        ll = cycle[cycle.size()-1];                   // ll <- max(cycle)
      }
  }                                           // }

  // std::cout<<"pred upd \n";
  if (update) { //is this safe?? memory leak?
    R[0]=boundaryV;
    R[2]=low;
    R[1]=nzc;
  }
  // std::cout<<"end \n";

  return(Rcpp::List::create(Rcpp::Named("VLL",Rcpp::wrap(VLL)),Rcpp::Named("cycle",Rcpp::wrap(CYCLE)),Rcpp::Named("ll",Rcpp::wrap(LL)),Rcpp::Named("II",Rcpp::wrap(II))));
  return wrap(repre);                         // repre
}


/*
  GET CYCLE REPRESENTATION
   
  Find representations of cycles using reduced boundary matrix. If reduced boundary matrix does not contain
  complexes sufficient to kill the cycle (produce an empty set by sequence of addition over Z2), nonreduced boundary
  matrix can be provided to calculate additional representation of such 'immortal' subcycles.
   
  • get_cycle_representation:
  get representation of a cycle (which is defined by a vector of simplex indices).
  Input:
  cycle ..................... integer vector of simplex indices (1-based);
  reduced_boundary_matrix ... reduced boundary matrix object (as produced by reduce_boundaryC);
  nonreduced_boundaries ..... non-reduced boundary matrix object, as produced by build_boundary. If provided,
  the algorithm will create new complexes in order to cover subcycles (if needed);
  update_rb ................. Boolean, if nonreduced_boundaries provided, should the provided reduced boundary
  matrix be amended to include the newly constructed complexes? This modifies the
  existing object!
  Output:
  resulting integer vector of column indices of reduced boundary matrix.
*/

// [[Rcpp::export]]
std::vector<int>
get_cycle_representation(
                         std::vector<int> cycle,
                         List             reduced_boundary_matrix,
                         const List       nonreduced_boundaries = 0,
                         bool             update_rb             = false
                         ) {
  bool new_complexes = nonreduced_boundaries.containsElementNamed("cmplx");
  // if nonreduced boundary matrix provided, subcycles can be resolved
    
  update_rb = new_complexes && update_rb;
    
  List b;
  if (new_complexes)
    b = nonreduced_boundaries["cmplx"];
    
  List rb                   = reduced_boundary_matrix["boundary"];
  IntegerVector nonzero_col = reduced_boundary_matrix["nonzero_col"];
  IntegerVector low         = reduced_boundary_matrix["low"];
    
  IntegerVector which_low(*std::max_element(low.begin(), low.end())+1, -1);
  for (size_t i = 0; i < low.length(); ++i)
    which_low[low[i]] = i;
  // which_low[i] ~ which non-zero column has entry '1' in row i
    
  std::vector<int> repre; // cycle representation
  std::sort(cycle.begin(), cycle.end()); // sort: increasing
    
  while (cycle.size() > 0)
    {
      R_CheckUserInterrupt();
      
      int cycle_max       = cycle.back();
      int cycle_which_max = which_low[cycle_max];
      // cycle_which_max ~ which non-zero column will cancel out the max element of cycle by addition
      
      if (cycle_which_max < 0) // if no such element present, subcycle representation needed...
        {
          if (new_complexes) // ...generate it, if we are allowed to
            {
              // std::cout << ".";
          
              std::vector<int> subrepre; // subcycle representation
              std::vector<int> subcycle = b[cycle_max-1]; // -1: re-indexing to 0-based
              std::sort(subcycle.begin(), subcycle.end()); // sort: increasing
          
              while (subcycle.size() > 0)
                {
                  R_CheckUserInterrupt();
            
                  int subcycle_max       = subcycle.back();
                  int subcycle_which_max = which_low[subcycle_max];
            
                  if (subcycle_which_max < 0)
                    {
                      warning("Failed with nonreduced boundaries provided. Not a cycle?");
                      std::vector<int> v(1, -1); return v;
                    }
            
                  subcycle = sym_diff(rb[subcycle_which_max], subcycle);
                  subrepre.push_back(nonzero_col[subcycle_which_max]);
                }
          
              if (subrepre[subrepre.size()-1] != cycle_max)
                subrepre.push_back(cycle_max);
          
              cycle = sym_diff(subrepre, cycle);
              repre.push_back(-cycle_max);
          
              // update reduced boundary matrix (internally)
              rb.push_back(subrepre); // add subcycle representation to boundary vectors
          
              cycle_which_max = rb.length()-1;
              which_low[cycle_max] = cycle_which_max;
          
              nonzero_col.push_back(-cycle_max);
              low.push_back(cycle_max);
          
            } else {
            warning("Failed without nonreduced boundaries provided. Not a cycle?");
            std::vector<int> v(1, -1); return v;
          }
        } else {
        cycle = sym_diff(rb[cycle_which_max], cycle);
        repre.push_back(nonzero_col[cycle_which_max]);
        cycle_max = cycle[0];
      }
    }
    
  if (update_rb)
    {
      reduced_boundary_matrix["boundary"]    = rb;
      reduced_boundary_matrix["nonzero_col"] = nonzero_col;
      reduced_boundary_matrix["low"]         = low;
    }
    
  return repre;
}

