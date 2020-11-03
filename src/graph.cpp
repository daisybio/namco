#include "graph.hpp"


// construct Graph object with num_nodes many adjacency lists
// initialize the list with otus of lvl 0
Graph::Graph(std::size_t num_nodes):
  adj_lists_(num_nodes){
}

std::size_t Graph::num_nodes() const{
  return adj_lists_.size();
}

std::vector<std::size_t> & Graph::get_otu_with_level_0_(){ 
  return Graph::otu_with_level_0_;
}

void Graph::add_edge(std::size_t source, std::size_t target){
  
  try{
    if(source >= num_nodes()){
      throw std::range_error("Source node" + std::to_string(source) + " does not exist.");
    }
    if(target >= num_nodes()){
      throw std::range_error("Target node" + std::to_string(target) + " does not exist.");
    }
    
    adj_lists_.at(source).emplace_back(target);
  }catch(std::exception &ex){
    forward_exception_to_r(ex);
  }catch(...){
    ::Rf_error("C++ exception in add_edge(). (unknown reason).");
  }
  
}

std::size_t Graph::compute_levels(std::vector<std::size_t> & levels) const{
  std::vector<std::size_t> in_degrees(num_nodes(), 0);
  for(const std::list<std::size_t> & adj_list : adj_lists_){
    for(std::size_t target : adj_list) {
      in_degrees.at(target)++;
    }
  }
  
  levels.clear();
  for (std::size_t node{0}; node < num_nodes(); node++) {
    levels.emplace_back(0);
  }
  
  std::list<std::size_t> minima;
  std::list<std::size_t> new_minima;
  for (std::size_t node{0}; node < num_nodes(); node++) {
    if (in_degrees.at(node) == 0) {
      minima.emplace_back(node);
    }
  }
  
  std::size_t current_level{1};
  std::size_t num_nodes_with_levels{0};
  
  while (minima.size() > 0) {
    new_minima.clear();
    for (std::size_t node : minima) {
      levels.at(node) = current_level;
      num_nodes_with_levels++;
      for (std::size_t target : adj_lists_.at(node)) {
        in_degrees.at(target)--;
        if (in_degrees.at(target) == 0) {
          new_minima.emplace_back(target);
        }
      }
    }
    minima = new_minima;
    current_level++;
  }

  return (current_level);
}


void Graph::print_graph(){
  for(int v=0;v < num_nodes(); v++){
    Rcpp::Rcout << "\n Adjacency list of vertex "<< v << "\n head ";
    for(auto x : adj_lists_[v]){
      Rcpp::Rcout << "-> " << x;
    }
    Rcpp::Rcout << "\n";
  }
}

