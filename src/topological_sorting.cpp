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

  for(std::size_t node_level_idx{0}; node_level_idx < levels.size(); node_level_idx++){
    //test if this node was saved as level_0 node earlier
    if(node_level_idx == g.get_otu_with_level_0_().at(otu_level_0_iterator)){
      normalized_levels.emplace_back(0);
      otu_level_0_iterator++;
      //handle case if iterator arrives at last element of level_0 list -> upper if would crash
      if(otu_level_0_iterator==g.get_otu_with_level_0_().size()){
        otu_level_0_iterator=g.get_otu_with_level_0_().size()-1;
      }
      continue;
    }
    normalized_levels.emplace_back((double)levels.at(node_level_idx)/max_level);
  }
  
  return normalized_levels;
}




/***R
otu_names<-rownames(otu)
out<-lapply(otu,function(x){
  return(calculate_topological_sorting(x,0.1))
})
out<-data.frame(out,row.names = otu_names)
#calculate_topological_sorting(otu$X18.SPF.CD,0.1)
*/
