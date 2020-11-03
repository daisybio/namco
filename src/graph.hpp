#ifndef INCLUDE_GRAPH_HPP_
#define INCLUDE_GRAPH_HPP_


class Graph {

public:
  Graph(std::size_t num_nodes);
  
  //Graph(const Graph & graph);
  
  void add_edge(std::size_t source, std::size_t target);
  
  std::size_t compute_levels(std::vector<std::size_t> & levels) const;
  
  std::size_t num_nodes() const;
  
  void print_graph();
  
  std::vector<std::size_t> & get_otu_with_level_0_();
  
private:
  
  std::vector<std::list<std::size_t>> adj_lists_;
  
  std::vector<std::size_t> otu_with_level_0_;
  
};


#include "graph.cpp"
#endif /* INCLUDE_GRAPH_HPP_ */