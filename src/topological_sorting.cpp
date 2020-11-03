#include <Rcpp.h>
#include <graph.hpp>

using namespace Rcpp;

/**
 * function to build a directed acyclic graph 
 * @param otu_list: vector with relative abundances for one sample/patient
 * @param cutoff: integer with cutoff to decide direction of edge
 * 
 * @return Graph object with nodes and edges as adjacency list  * 
 */

Graph construct_graph(std::vector<std::size_t> otu_list, double cutoff){
  //generate object Graph with as many nodes as otus
  Graph otu_graph = Graph(otu_list.size());
  
  for(std::size_t source_otu{0}; source_otu < otu_list.size() -1; source_otu++){
    double source_abund = otu_list.at(source_otu);
    if(source_abund == 0){
      otu_graph.get_otu_with_level_0_().emplace_back(source_otu);
      continue;
    }
    
    for(std::size_t target_otu{source_otu+1}; target_otu < otu_list.size(); target_otu++){
      double target_abund = otu_list.at(target_otu);
      if(target_abund == 0) {continue;}
      
      double logratio = std::log(source_abund/target_abund);
      if(logratio < -cutoff){
        otu_graph.add_edge(source_otu, target_otu);
      }
      else if(logratio > cutoff){
        otu_graph.add_edge(target_otu, source_otu);
      }
    }
  }
  
  //check last element of otu-list if it has abundance of 0
  if(otu_list.back() == 0){
    otu_graph.get_otu_with_level_0_().emplace_back(otu_list.size()-1);
  }
  
  /*for(int i = 0; i < otu_graph.get_otu_with_level_0_().size(); i++){
    Rcpp::Rcout<<otu_graph.get_otu_with_level_0_().at(i)<<"\n";
  }*/
    
  return(otu_graph);
  
}


/**
 * function to calculate topological sorting of a given vector (patient/sample);
 * builds directed acyclic graph and computes levels for each node given a cutoff
 * 
 * @param sample: vector with otu abundance values 
 * @param cutoff: double
 * 
 */
// [[Rcpp::export]]
std::vector<double> calculate_topological_sorting(std::vector<std::size_t> sample,double cutoff){
  Graph g = construct_graph(sample,cutoff);
  std::vector<std::size_t> levels;
  double max_level = (double)g.compute_levels(levels);

  // normalize level values between 0 and 1
  // replace levels with 0, which are saved in level_0 vector
  std::vector<double> normalized_levels;
  std::size_t otu_level_0_iterator{0};
  for(std::size_t i{0}; i < levels.size(); i++){
    if(i == g.get_otu_with_level_0_().at(otu_level_0_iterator)){
      normalized_levels.emplace_back(0);
      otu_level_0_iterator++;
      continue;
    }
    normalized_levels.emplace_back((double)levels.at(i)/max_level);
  }
  
  return normalized_levels;
}




/***R
out<-calculate_topological_sorting(t,0.1)
*/
