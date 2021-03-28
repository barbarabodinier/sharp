Graph=function(adjacency, node_label=NULL, node_colour=NULL, node_shape=NULL,
               weighted=NULL, satellites=FALSE){
  # Checking input values (weighted)
  if (!is.null(weighted)){
    if (!weighted%in%c(TRUE,FALSE)){
      stop("Invalid input for argument 'weighted'. Possible values are: NULL, TRUE or FALSE.")
    }
  }

  # Extracting the adjacency matrix from the output of CalibrateNetwork()
  if (!is.matrix(adjacency)){
    adjacency=Adjacency(stability=adjacency)
  }

  # Setting row and column names if none
  if (is.null(rownames(adjacency))&(is.null(colnames(adjacency)))){
    rownames(adjacency)=colnames(adjacency)=paste0("var",1:ncol(adjacency))
  } else {
    if (is.null(rownames(adjacency))){
      rownames(adjacency)=colnames(adjacency)
    }
    if (is.null(colnames(adjacency))){
      colnames(adjacency)=rownames(adjacency)
    }
  }

  # Checking input values (node label)
  if (!is.null(node_label)){
    if (length(node_label)!=ncol(adjacency)){
      stop(paste0("Invalid input for argument 'node_label'. It must be a vector of length ",
                  ncol(adjacency)," (the same as the number of nodes in the adjacency matrix)."))
    }
  }

  # Checking input values (node colour)
  if (!is.null(node_colour)){
    if (length(node_colour)==1){
      node_colour=rep(node_colour, ncol(adjacency))
    } else {
      if (length(node_colour)!=ncol(adjacency)){
        stop(paste0("Invalid input for argument 'node_colour'. It must be a vector of length ",
                    ncol(adjacency)," (the same as the number of nodes in the adjacency matrix)."))
      }
    }
  }

  # Checking input values (node shape)
  if (!is.null(node_shape)){
    if (length(node_shape)==1){
      node_shape=rep(node_shape, ncol(adjacency))
    } else {
      if (length(node_shape)!=ncol(adjacency)){
        stop(paste0("Invalid input for argument 'node_shape'. It must be a vector of length ",
                    ncol(adjacency)," (the same as the number of nodes in the adjacency matrix)."))
      }
    }

    if (!any(node_shape%in%c("circle", "square", "triangle", "star"))){
      stop(paste0("Invalid input for argument 'node_shape'. Possible values for the entries of the vector are: circle, square, triangle or star."))
    }
  }

  # Adding shapes if required
  if (!is.null(node_shape)){
    if (any(node_shape=="star")){
      igraph::add_shape("star", clip=igraph::shape_noclip,
                        plot=mystar, parameters=list(vertex.norays=5))
    }

    if (any(node_shape=="triangle")){
      igraph::add_shape("triangle", clip=igraph::shape_noclip,
                        plot=mytriangle)
    }
  }

  # Default node colours
  if (is.null(node_colour)){
    node_colour=rep("skyblue", ncol(adjacency))
  }

  # Default node shapes
  if (is.null(node_shape)){
    node_shape=rep("circle", ncol(adjacency))
  }

  # Default node labels
  if (is.null(node_label)){
    node_label=colnames(adjacency)
  }

  # Formatting node characteristics
  names(node_colour)=colnames(adjacency)
  names(node_label)=colnames(adjacency)
  names(node_shape)=colnames(adjacency)

  # Formatting adjacency matrix
  if (!is.null(weighted)){
    if (!weighted){
      adjacency=ifelse(adjacency!=0,yes=1,no=0)
      weighted=NULL
    }
  }

  # Estimating igraph object
  mygraph=igraph::graph_from_adjacency_matrix(adjacency, mode="undirected", weighted=weighted)
  mydegrees=igraph::degree(mygraph)

  # Including/excluding satellites (nodes with no edges)
  if (!satellites){
    mygraph=igraph::delete.vertices(mygraph, v=names(mydegrees)[mydegrees==0])
  }

  # Formatting vertices
  mydegrees=igraph::degree(mygraph)
  igraph::V(mygraph)$size=as.numeric(as.character(cut(mydegrees, breaks=4, labels=c(3, 4, 5, 6))))
  igraph::V(mygraph)$label=node_label[igraph::V(mygraph)$name]
  igraph::V(mygraph)$color=node_colour[igraph::V(mygraph)$name]
  igraph::V(mygraph)$shape=node_shape[igraph::V(mygraph)$name]
  igraph::V(mygraph)$frame.color=igraph::V(mygraph)$color
  igraph::V(mygraph)$label.family="sans"
  igraph::V(mygraph)$label.cex=as.numeric(as.character(cut(mydegrees, breaks=4, labels=c(0.4, 0.45, 0.5, 0.55))))
  igraph::V(mygraph)$label.color="grey20"

  # Formatting edges
  igraph::E(mygraph)$color="grey60"
  if (is.null(weighted)){
    igraph::E(mygraph)$width=0.5
  } else {
    igraph::E(mygraph)$width=igraph::E(mygraph)$weight
  }

  return(mygraph)
}


mystar=function(coords, v=NULL, params) {
  vertex.color=params("vertex", "color")
  if ((length(vertex.color)!=1)&!is.null(v)) {
    vertex.color=vertex.color[v]
  }
  vertex.size=1/200*params("vertex", "size")
  if ((length(vertex.size)!=1)&!is.null(v)) {
    vertex.size=vertex.size[v]
  }
  norays=params("vertex", "norays")
  if ((length(norays)!=1)&!is.null(v)) {
    norays=norays[v]
  }

  mapply(coords[,1], coords[,2], vertex.color, vertex.size, norays,
         FUN=function(x, y, bg, size, nor) {
           graphics::symbols(x=x, y=y, fg=bg, bg=bg,
                             stars=matrix(c(size,size/2), nrow=1, ncol=nor*2),
                             add=TRUE, inches=FALSE)})
}


mytriangle=function(coords, v=NULL, params) {
  vertex.color=params("vertex", "color")
  if ((length(vertex.color)!=1)&!is.null(v)) {
    vertex.color=vertex.color[v]
  }
  vertex.size=1/200*params("vertex", "size")
  if ((length(vertex.size)!=1)&!is.null(v)) {
    vertex.size=vertex.size[v]
  }
  norays=3

  mapply(coords[,1], coords[,2], vertex.color, vertex.size, norays,
         FUN=function(x, y, bg, size, nor) {
           graphics::symbols(x=x, y=y, fg=bg, bg=bg,
                             stars=matrix(c(size,size/2), nrow=1, ncol=nor*2),
                             add=TRUE, inches=FALSE)})
}

