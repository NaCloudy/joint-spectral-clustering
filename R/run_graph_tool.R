### Make multilayer adjacency matrices into single multilayer graphml
# filename = "../temp/Example1b.graphml"
adjacencies_to_graphml <- function(Adj_list, filename) {
  
  require(igraph)
  require(Matrix)
  
  m <- length(Adj_list)
  n <- ncol(Adj_list[[1]])
  gs1 <- lapply(1:m, function(i) {
    A <- Adj_list[[i]]
    diag(A) <- 1
    g <- graph_from_adjacency_matrix(i*A, mode = "undirected", weighted = "weight")
    get.data.frame(g)
  })
  
  df1 <- Reduce(rbind, gs1)
  
  g_all1 <- graph_from_edgelist(as.matrix(df1[,1:2]), directed = F)
  E(g_all1)$weight <- df1[,3]
  g_all1 <- delete_edges(g_all1, E(g_all1)[which_loop(g_all1)])
  igraph::write.graph(g_all1,filename, format = "graphml")
  
}


run_graph_tool <- function(Adj_list, K, id = NULL) {
  if(is.null(id)) {
    id = round(runif(1)*10000000)
  }
  require(reticulate)
  # use_python("/usr/local/bin/python3") # did not work
  #graph_file <-  paste("temp/multilayer", id, ".graphml", sep="")
  #result_file <- paste("temp/result", id,".csv", sep = "")
  #adjacencies_to_graphml(Adj_list, graph_file)
  #command_gt <- paste("/usr/local/bin/python3 Python/graphtool-script.py", K, id)
  
  python_path <- "/root/miniconda3/envs/gt/bin/python3"
  graph_file <-  paste("temp/multilayer", id, ".graphml", sep="")
  result_file <- paste("temp/result", id,".csv", sep = "")
  adjacencies_to_graphml(Adj_list, graph_file)
  command_gt <- paste(python_path, "Python/graphtool-script.py", K, id)
  
  system(command_gt)
  result = read.csv(result_file)
  system(paste("rm", graph_file))
  system(paste("rm", result_file))
  as.numeric(factor(result$X0))
}
