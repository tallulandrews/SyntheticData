# Generation of Synthetic Spatial Transcriptomics Data (VISIUM)
# Written by: Adam Southworth
# Under the Supervision of: Dr. Tallulah Andrews
# Date: 2023-08-24



#' create_coords_VISIUM_hexes
#' Creates a base grid of spots with no data; the default parameters represent an average spatial transcriptomics slice size. 
#' @param nrow The number of hexes along the x-axis of the grid, numeric
#' @param ncol The number of hexes along the y-axis of the grid, numeric
#' @param pixels_x The number of pixels along the x-axis of a hypothetical image associated with the grid, numeric
#' @param pixels_y The number of pixels along the x-axis of a hypothetical image associated with the grid, numeric
#'
#' @return A dataframe with barcoded spots and their associated X and Y coordinates both on the grid (hexagonal spot) and image(pixel). This return, and any modifications made to it by other helper functions, is referred to as a coordinate dataframe in other docstrings.
create_coords_VISIUM_hexes <- function(nrow=77, ncol=127, pixels_x=c(1980,24014), pixels_y=c(2780, 23665)) { 
  row_coords1 <- seq(1, nrow, by =2)
  row_coords2 <- seq(2, nrow, by =2)
  col_coords1 <- seq(1, ncol, by =2)
  col_coords2 <- seq(2, ncol, by =2)
  
  all_coords <- cbind(
    c(rep(row_coords1, times=length(col_coords1)), rep(row_coords2, times=length(col_coords2))), 
    c(rep(col_coords1, each=length(row_coords1)), rep(col_coords2, each=length(row_coords2)))
  )
  
  rownames(all_coords) <- paste("Spot", 1:nrow(all_coords), sep="")
  
  pixels_x <- round(seq(pixels_x[1], pixels_x[2], length = max(c(row_coords1, row_coords2))))
  
  pixels_y <- round(seq(pixels_y[1], pixels_y[2], length = max(c(col_coords1, col_coords2))))
  
  all_coords_pixels <- cbind(pixels_x[all_coords[,1]], pixels_y[all_coords[,2]])
  
  return(cbind(rep(1, nrow(all_coords)), all_coords, all_coords_pixels))
}

#' add_uniform
#' Applies a uniform label across all hexes. 
#' @param coords A coordinate dataframe
#' @param id The label to apply, string
#' @param plot Whether or not to plot the grid after applying the label. Usually used for testing, boolean
#'
#' @return A modififed coordinate dataframe
add_uniform <- function(coords, id = "ashared1", plot = FALSE) {
  coords <- data.frame(coords)
  
  coords$data <- id
  
  if (plot) {
    plot_hexes(coords)
  } 
  
  return(coords)
}

#' add_single_stripe 
#' Divides the grid into a number of vertical (or horizontal) stripes and changes the label of one of those stripes. 
#' @param coords A coordinate dataframe
#' @param vertical Creates a vertical stripe if TRUE, horizontal otherwise, boolean
#' @param nstripes The number of stripes the grid should be divided into. Effectively this is the width of the stripe, integer
#' @param target The target stripe for which the label is changed. Effectively this is the location of the stripe on the grid, integer in a range from 1 to nstripes
#' @param id The label to apply, string
#' @param plot Whether or not to plot the grid after applying the labels. Usually used for testing, boolean
#'
#' @return A modified coordinate dataframe
add_single_stripe = function(coords, vertical=TRUE, nstripes=3,target = "1", id = "ashared1", plot=FALSE) { 
  if (vertical) {
    out <- split(rownames(coords), cut(coords[,2], nstripes)); 
    # cut converts the X (column) IDs into factors at the number of points specified by the number of stripes, then returns a list of the hexes split by this
  } else {
    out <- split(rownames(coords), cut(coords[,3], nstripes));
    # cut converts the Y (row) IDs into factors at the specified number of points then returns a list of vectors (inside the vectors are hexes)
  }
  
  names(out) <- paste("stripe", if (vertical){"V"}else{"H"}, 1:nstripes, sep="_")
  
  for (i in 1:length(out)) {
    if (grepl(target, names(out)[i])) {  # if the target is in the stripe we want to mark those ones only. Others are NAs
      names(out)[i] <- id
    }
    else {
      names(out)[i] <- NA
    }
  }
  
  coords <- data.frame(coords)
  
  for (i in 1:length(out)) {
    if (!is.na(names(out)[i])) {
      coords$data[ rownames(coords) %in% out[[i]] ] <- names(out)[i]
    }
  }
  
  if (plot) {
    plot_hexes(coords)
  } 
  return(coords)
}


#' add_rectangle
#' Splits up the grid into smaller boxes and applies a label to one of them.
#' @param coords A coordinate dataframe
#' @param xfrac The number of columns the grid will be split into, integer
#' @param targx Which column of those split by xfrac will be selected, integer
#' @param yfrac The number of rows the grid will be split into, integer
#' @param targy Which column of those split by yfrac will be selected, integer
#' @param id The label to apply, string
#' @param plot Whether or not to plot the grid after applying the labels. Usually used for testing, boolean
#'
#' @return A modified coordinate dataframe
add_rectangle = function(coords, xfrac = 3, targx = 1, yfrac = 3, targy = 1, id = "ashared1", plot=FALSE) { 
  xdivided <- split(rownames(coords), cut(coords[,2], xfrac)) # splits up the x values of the grid by the fraction we specified
  ydivided <- split(rownames(coords), cut(coords[,3], yfrac)) # same but for y
  
  out <- vector(mode = "list", length = 0)
  
  lenx <- length(xdivided[[targx]])
  leny <- length(ydivided[[targy]])
  
  if (leny > lenx)  { # checks the region of overlap between the two segments if the segment on the y-axis is bigger
    for (i in 1:length(xdivided[[targx]])) { 
      if (is.element(xdivided[[targx]][i], ydivided[[targy]])) {
        out <- c(out, xdivided[[targx]][i])
      }
    }
  }
  else { # checks the region of overlap between the two segments if the segment on the x-axis is bigger or if they are the same size
    for (i in 1:length(ydivided[[targy]])) { 
      if (is.element(ydivided[[targy]][i], xdivided[[targx]])) {
        out <- c(out, ydivided[[targy]][i])
      }
    }
  }
  
  coords <- data.frame(coords)
  
  coords$data[is.element(rownames(coords), out)] <- id  
  
  if (plot) {
    plot_hexes(coords)
  }
  return(coords)
}

#' add_random_points
#' Labels a totally random assortment of spots on the grid. 
#' @param coords A coordinate dataframe
#' @param nspots The number of spots to be labeled, numeric
#' @param nsets The number of times to select from all spots, integer
#' @param id The label to apply, string
#' @param plot Whether or not to plot the grid after applying the labels. Usually used for testing, boolean
#' @return A modified coordinate dataframe
add_random_points <- function(coords, nspots=0.1*nrow(coords), nsets=1, id="ashared1", plot=FALSE) {  # by default selects 10% of the spots once
  coords <- data.frame(coords)
  
  for (i in 1:nsets) {
    spots <- sample(1:nrow(coords), nspots)
    coords$data[spots] <- id
  }
  if (plot) {
    plot_hexes(coords)
  }
  return(coords)
}

#' add_stripe_or_rect
#' Using helper functions, adds either a stripe or a rectangle (can be random which shape, or specified) with random size and location.
#' @param slice A coordinate dataframe
#' @param id The label to apply, string
#' @param shape The shape to be made, may either be 1 to specify stripe or 2 to specify rectangle. Otherwise the shape chosen will be random, integer
#'
#' @return A modified coordinate dataframe
add_stripe_or_rect <- function(slice, id="ashared1", shape=0) {
  
  if ((sample(c(TRUE, FALSE), 1) && shape == 0) || shape == 1) {  # adds stripe 
    randdirection <- sample(c(TRUE, FALSE), 1)  # randomly determines if regions split horizontally or vertically
    randstripes <- sample(c(6:8), 1)    # randomly determines the number of splits
    randtarget <- sample(c(1:randstripes), 1)    # picks which part to target
    
    slice <- add_single_stripe(slice, randdirection, randstripes, randtarget, id=id)
    
  }
  else {  # adds rectangle
    randxfrac <- sample(c(3:6), 1) # randomly determine parameters. The limits for randxfrac and randyfrac can be changed for more variable stripe sizes
    randyfrac <- sample(c(3:6), 1)
    randtargx <- sample(c(1:randxfrac), 1)
    randtargy <- sample(c(1:randyfrac), 1)
    
    slice <- add_rectangle(slice, randxfrac, randtargx, randyfrac, randtargy, id=id)
    
  }
  return(slice)
}

#' add_neighbour
#' Adds a cluster with a given label in a neighbouring (or potentially overlapping) position relative to clusters of other labels.
#' @param slice A coordinate dataframe
#' @param id_originals The labels in the grid that the new cluster will neighbour, vector of strings
#' @param id_neighbour The label of the new cluster, string
#' @param stripe The type of shape being added to the grid, stripe if TRUE, rectangle otherwise, boolean
#'
#' @return A modified coordinate dataframe
add_neighbour <- function(slice, id_originals, id_neighbour, stripe) {
  neighbours_original <- c()
  is_neighbour <- 0
  counter <- 0
  
  for (i in 1:length(id_originals))
    neighbours_original <- c(neighbours_original, get_cluster_neighbours(slice, id_originals[i]))
  
  if (stripe == TRUE) {
    while (is_neighbour == 0) {
      counter = counter + 1
      temp <- slice  # to store the original slice in case we need to switch back to it (i.e. if the addition of the neighbour fails)
  
      slice <- add_stripe_or_rect(slice, id=id_neighbour, shape=1)
      new_neighbour <- rownames(slice[slice$data == id_neighbour, ])
      
      if (length(intersect(new_neighbour, neighbours_original)) > 0 & nlevels(as.factor(slice$data)) > nlevels(as.factor(temp$data))) {  # do we intersect with the neighbours and are also not totally overlapping?
        is_neighbour <- 1
      }
      else {
        slice <- temp # reset the slice and try again!
      }
    }
  }
  else if (stripe == FALSE) {  # similar to before but add a rectangle
    while (is_neighbour == 0)  {
      counter = counter + 1
      temp <- slice 
      
      slice <- add_stripe_or_rect(slice, id=id_neighbour, shape=2)
      new_neighbour <- rownames(slice[slice$data == id_neighbour, ])
      
      if (length(intersect(new_neighbour, neighbours_original)) > 0 & nlevels(as.factor(slice$data)) > nlevels(as.factor(temp$data))) {
        is_neighbour <- 1
      }
      else {
        slice <- temp  
      }
    }
  }
  return(slice)
}

#' get_cluster_neighbours
#' Finds all spots that neighbour all regions of a particular label.
#' @param slice A coordinate dataframe
#' @param id_cluster The label to try and find neighbours to, string
#'
#' @return A vector of neighbouring spot barcodes
get_cluster_neighbours <- function(slice, id_cluster) {
  cluster <- slice[slice$data == id_cluster, ]
  neighbours <- c()
  
  cluster_max_x <- max(cluster$X2)  # finds maximum and minimum values of each cluster, we will increment by 1 to find neighbours
  cluster_min_x <- min(cluster$X2)
  cluster_max_y <- max(cluster$X3)
  cluster_min_y <- min(cluster$X3)
  
  if (cluster_max_x + 1 %in% cluster$X2) {  # will the area just beyond the maximum or minimum values actually exist?
    neighbours <- c(neighbours, rownames(slice[slice$X2 == (cluster_max_x + 1) & (slice$X3 >= cluster_min_y) & (slice$X3 <= cluster_max_y), ]))
  }
  
  if (cluster_min_x - 1 %in% cluster$X2) {  
    neighbours <- c(neighbours, rownames(slice[slice$X2 == (cluster_min_x - 1) & (slice$X3 >= cluster_min_y) & (slice$X3 <= cluster_max_y), ]))
  }
  
  if (cluster_max_y + 1 %in% cluster$X3) { 
    neighbours <- c(neighbours, rownames(slice[slice$X3 == (cluster_max_y + 1) & (slice$X2 >= cluster_min_x) & (slice$X2 <= cluster_max_x), ]))
  }
  
  if (cluster_min_y - 1 %in% cluster$X3) {  
    neighbours <- c(neighbours, rownames(slice[slice$X3 == (cluster_min_y - 1) & (slice$X2 >= cluster_min_x) & (slice$X2 <= cluster_max_x), ]))
  }
  
  return(neighbours)
}

#' degree_of_similarity
#' Compares two slices to see how similar they are with respect to a particular label or labels
#' @param slice1 A coordinate dataframe
#' @param slice2 A coordinate dataframe
#' @param id1 The label on slice 1 to be used in the comparison, string
#' @param id2 The label on slice 2 to be used in the comparison, string
#'
#' @return A list contaning the exact number of similar hexes and the fraction of each original slice that is spatially similar to the other with respect to the given label or labels
degree_of_similarity <- function(slice1, slice2, id1, id2) {
  
  unique_spots_slice1 <- rownames(slice1)[slice1$data == id1]
  unique_spots_slice2 <- rownames(slice2)[slice2$data == id2]
  
  similar_spatials <- length(intersect(unique_spots_slice1, unique_spots_slice2)) # the number of hexes of the given id (id may differ depending on slice) who are found in the same spot in both slices
  
  similar_per_slice1 <- similar_spatials / length(unique_spots_slice1)
  similar_per_slice2 <- similar_spatials / length(unique_spots_slice2)
  
  output <- list(similar_spatials, similar_per_slice1, similar_per_slice2)
  names(output) <- c("Sim_Count", "SimFrac_Slice1", "SimFrac_Slice2")
  
  return(output)
}

#' plot_hexes
#' A function that plots a single slice's hexagonal grid based on the data in each spot; note that colouring is inaccurate in some instances. 
#' @param coords A coordinate dataframe
#'
#' @return NULL
plot_hexes <- function(coords) {
  
  plot(coords$X2[!is.na(coords$data)], coords$X3[!is.na(coords$data)], 
       pch=16, col=as.factor(coords$data),
       xlim=c(0, max(coords$X2)), ylim=c(0, max(coords$X3)), 
       xlab="", ylab="")
}

#' plot_side_by_side
#' Plots two slices' hexagonal grids side by side based on the data in each spot; note that colouring is inaccurate in some instances. 
#' @param doublecoords A list of two coordinate dataframes
#'
#' @return NULL
plot_side_by_side <- function(doublecoords) {  # plots two hex grids side by side based on the data in each cell. Takes in a list of two hex grids
  coords1 <- doublecoords[[1]]
  coords2 <- doublecoords[[2]]
  
  par(mfrow=c(1,2))
  
  plot(coords1$X2[!is.na(coords1$data)], coords1$X3[!is.na(coords1$data)], 
       pch=16, col=as.factor(coords1$data), 
       xlim=c(0, max(coords1$X2)), ylim=c(0, max(coords1$X3)), 
       xlab="", ylab="") 
  
  plot(coords2$X2[!is.na(coords2$data)], coords2$X3[!is.na(coords2$data)], 
       pch=16, col=as.factor(coords2$data), 
       xlim=c(0, max(coords2$X2)), ylim=c(0, max(coords2$X3)), 
       xlab="", ylab="") 
  
}

#' make_idents_ex1
#' Creates two slices (coordinate dataframes) with identical and uniform data in both
#' @return A list of two coordinate dataframes
make_idents_ex1 <- function() {  
  slice1 <- create_coords_VISIUM_hexes()  # create two blank hex sets to add data to
  slice2 <- create_coords_VISIUM_hexes()
  
  slice1 <- add_uniform(slice1) # apply identical data
  slice2 <- add_uniform(slice2)
  
  slices = list(slice1, slice2)
  names(slices) <- c("Slice 1", "Slice 2")
  
  return(slices)
}

#' make_opposites_ex2
#' Creates two slices with opposite and uniform data in both
#' @return A list of two coordinate dataframes
make_opposites_ex2 <- function() { 
  slice1 <- create_coords_VISIUM_hexes()  # create two blank hex sets to add data to
  slice2 <- create_coords_VISIUM_hexes()
  
  slice1 <- add_uniform(slice1, id="buni1")  # apply opposite data uniformly
  slice2 <- add_uniform(slice2, id="buni2")
  
  slices = list(slice1, slice2)
  names(slices) <- c("Slice 1", "Slice 2")
  
  return(slices)
}

#' make_simdata_simspat_ex3
#' Creates two slices with lots of spots matching each other; all matching spots are spatially identical
#' @return A list of two coordinate dataframes
make_simdata_simspat_ex3 <- function() {
  slice1 <- create_coords_VISIUM_hexes()  # create two blank hex sets to add data to
  slice2 <- create_coords_VISIUM_hexes()
  
  slice1 <- add_uniform(slice1)  # apply identical data
  slice2 <- add_uniform(slice2)
  
  randdirection <- sample(c(TRUE, FALSE), 1)  # # randomly determines if regions split horizontally or vertically. The random parameters are used on both slices
  randstripes <- sample(c(3:6), 1)    # randomly determines the number of splits
  randtarget <- sample(c(1:randstripes), 1)    # picks which part to target
  
  slice1 <- add_single_stripe(slice1, nstripes=randstripes, target=randtarget, vertical=randdirection, id="buni1")
  slice2 <- add_single_stripe(slice2, nstripes=randstripes, target=randtarget, vertical=randdirection, id="buni2") 
  
  slices = list(slice1, slice2)
  names(slices) <- c("Slice 1", "Slice 2")
  
  return(slices)
}

#' make_diffdata_simspat_ex4
#' Creates two slices with very little cells matching each other; all matching cells are spatially identical
#' @return A list of two coordinate dataframes
make_diffdata_simspat_ex4 <- function() { 
  slice1 <- create_coords_VISIUM_hexes()  # create two blank hex sets to add data to
  slice2 <- create_coords_VISIUM_hexes()
  
  slice1 <- add_uniform(slice1, id="buni1")  # apply opposite data uniformly
  slice2 <- add_uniform(slice2, id="buni2")
  
  randdirection <- sample(c(TRUE, FALSE), 1)  # randomly determines if regions split horizontally or vertically
  randstripes <- sample(c(3:6), 1)    # randomly determines the number of splits
  randtarget <- sample(c(1:randstripes), 1)    # picks which part to target
  
  slice1 <- add_single_stripe(slice1, nstripes=randstripes, target=randtarget, vertical=randdirection)
  slice2 <- add_single_stripe(slice2, nstripes=randstripes, target=randtarget, vertical=randdirection)
  
  slices = list(slice1, slice2)
  names(slices) <- c("Slice 1", "Slice 2")
  
  return(slices)
}

#' make_simdata_diffspat_ex5
#' Creates two slices with lots of cells matching each other; all matching cells are not spatially identical
#' @return A list of two coordinate dataframes
make_simdata_diffspat_ex5 <- function() {
  slice1 <- create_coords_VISIUM_hexes()  # create two blank hex sets to add data to
  slice2 <- create_coords_VISIUM_hexes()
  
  slice1 <- add_uniform(slice1)  # apply identical uniform data
  slice2 <- add_uniform(slice2)

  slice1 <- add_stripe_or_rect(slice1, id="buni1", shape=2)
  
  is_diffspat <- 0  # value to maintain a while loop to keep generating the second slice until it's spatially different from slice 1
  
  while (is_diffspat == 0) {

    slice2 <- add_stripe_or_rect(slice2, id="buni2", shape=2)
    
    if (degree_of_similarity(slice1, slice2, "buni1", "buni2")[[1]] == 0) {
      is_diffspat <- 1
    } else {
      slice2 <- add_uniform(slice2)  # resets the slice to try and add a new rectangle to it
    }
  }
  
  slices = list(slice1, slice2)
  names(slices) <- c("Slice 1", "Slice 2")
  
  return(slices)
}

#' make_diffdata_diffspat_ex6
#' Creates two slices with very little cells matching each other; all matching cells are not spatially identical
#' @return A list of two coordinate dataframes
make_diffdata_diffspat_ex6 <- function() {  # returns two slices with mostly different data and different spatial orientation
  slice1 <- create_coords_VISIUM_hexes()  # create two blank hex sets to add data to
  slice2 <- create_coords_VISIUM_hexes()
  
  slice1 <- add_uniform(slice1, id = "buni1")  # apply identical uniform data
  slice2 <- add_uniform(slice2, id = "buni2")
  
  slice1 <- add_stripe_or_rect(slice1, id="ashared1", shape=2)
  
  is_diffspat <- 0  # value to maintain a while loop to keep generating the second slice until it's spatially different from slice 1
  
  while (is_diffspat == 0) {

    slice2 <- add_stripe_or_rect(slice2, id="ashared1", shape=2)
    
    
    if (degree_of_similarity(slice1, slice2, id1="ashared1", id2="ashared1")[[1]] == 0) {
      is_diffspat <- 1
    } else {
      slice2 <- add_uniform(slice2, id="buni2")  # resets the slice to try and add a new rectangle to it
    }
  }
  
  slices = list(slice1, slice2)
  names(slices) <- c("Slice 1", "Slice 2")
  
  return(slices)
}

#' make_randoms_ex7
#' Creates two slices with randomly similar regions and no significant spatial correlation between them
#' @return A list of two coordinate dataframes
make_randoms_ex7 <- function() {  # returns two slices with random similarity and spatially random
  slice1 <- create_coords_VISIUM_hexes()  # create two blank hex sets to add data to
  slice2 <- create_coords_VISIUM_hexes()
  
  slice1 <- add_uniform(slice1, id="buni1")  # apply opposite data uniformly
  slice2 <- add_uniform(slice2, id="buni2")
  
  slice1 <- add_random_points(slice1)  # add the random shared points
  slice2 <- add_random_points(slice2)
  
  slices = list(slice1, slice2)
  names(slices) <- c("Slice 1", "Slice 2")
  
  return(slices)
}

#' make_somedata_diffloc_diffsize_ex8
#' Creates two slices with some similar regions; similar regions have different sizes and little spatial similarity
#' @return A list of two coordinate dataframes
make_somedata_diffloc_diffsize_ex8 <- function() {  # returns two slices with random similarity and spatially random
  slice1 <- create_coords_VISIUM_hexes()  # create two blank hex sets to add data to
  slice2 <- create_coords_VISIUM_hexes()
  
  slice1 <- add_uniform(slice1, id="buni1")
  slice2 <- add_uniform(slice2)
  
  slice1 <- add_stripe_or_rect(slice1, id="ashared1", shape=2)
  
  is_simspat <- 0
  
  while (is_simspat == 0) {
    
    slice2 <- add_stripe_or_rect(slice2, id="buni2", shape=2)
    
    if (degree_of_similarity(slice1, slice2, id1="ashared1", id2="ashared1")[[2]] == 1) { 
      is_simspat <- 1
    }
    else {
      slice2 <- add_uniform(slice2, id="ashared1")  # resets the slice to try and add a new rectangle to it
    }
  }
  
  slices = list(slice1, slice2)
  names(slices) <- c("Slice 1", "Slice 2")
  
  return(slices)
}

#' make_somesim_spatunique_ex9
#' Creates two slices with multiple groups of similar regions that are in different locations in proximity to each other
#' @return A list of two coordinate dataframes
make_somesim_spatunique_ex9 <- function() {
  slice1 <- create_coords_VISIUM_hexes()  # create two blank hex sets to add data to
  slice2 <- create_coords_VISIUM_hexes()
  
  slice1 <- add_uniform(slice1, id="buni1")
  slice2 <- add_uniform(slice2, id="buni2")
  
  slice1 <- add_stripe_or_rect(slice1, id="ashared1", shape=2) # adding the first set of idential but spatially unique data. The others will be in proximity to it (but in different spatial orientations)
  
  slice1 <- add_neighbour(slice = slice1, id_originals = c("ashared1"), id_neighbour = "cshared2", stripe = sample(c(TRUE, FALSE), 1)) # adds either a rectangle or a stripe first
  
  slice1 <- add_neighbour(slice = slice1, id_originals = c("ashared1", "cshared2"), id_neighbour = "dshared3", stripe = sample(c(TRUE, FALSE), 1)) # does the same thing with the next cluster
  
  # Now we will create slice 2 with a similar process
  
  slice2 <- add_stripe_or_rect(slice2, id="ashared1", shape=2) # adding the first set of identical but spatially unique data. The others will be in proximity to it (but in different spatial orientations)
  
  slice2 <- add_neighbour(slice = slice2, id_originals = c("ashared1"), id_neighbour = "cshared2", stripe = sample(c(TRUE, FALSE), 1)) # adds either a rectangle or stripe first
  
  slice2 <- add_neighbour(slice = slice2, id_originals = c("ashared1","cshared2"), id_neighbour = "dshared3", stripe = sample(c(TRUE, FALSE), 1)) # does the same thing with the next cluster
  
  slices = list(slice1, slice2)
  names(slices) <- c("Slice 1", "Slice 2")
  
  return(slices)
}

#' make_somesim_spatsim_ex10
#' Creates two slices with multiple groups of similar regions that are in the same locations in proximity to each other
#' @return A list of two coordinate dataframes
make_somesim_spatsim_ex10 <- function() {
  slice1 <- create_coords_VISIUM_hexes()  # create two blank hex sets to add data to
  slice2 <- create_coords_VISIUM_hexes()
  
  slice1 <- add_uniform(slice1, id="buni1")
  
  slice1 <- add_stripe_or_rect(slice1, id="ashared1", shape=2)
  
  slice1 <- add_neighbour(slice = slice1, id_originals = c("ashared1"), id_neighbour = "cshared2", stripe = sample(c(TRUE, FALSE), 1))
  
  slice1 <- add_neighbour(slice = slice1, id_originals = c("ashared1", "cshared2"), id_neighbour = "dshared3", stripe = sample(c(TRUE, FALSE), 1))
  
  slice2 <- slice1  # slice 2 has the same spatial orientation of clusters but we still have to change some data
  
  slice2$data[slice2$data == "buni1"] <- "buni2"
  
  slices = list(slice1, slice2)
  names(slices) <- c("Slice 1", "Slice 2")
  
  return(slices)
}

#' make_somesim_somespat_ex11
#' Creates two slices with multiple groups of similar regions who have some similarities and differences spatially
#' @return A list of two coordinate dataframes
make_somesim_somespat_ex11 <- function() {
  slice1 <- create_coords_VISIUM_hexes()  # create two blank hex sets to add data to
  slice2 <- create_coords_VISIUM_hexes()
  
  slice1 <- add_uniform(slice1, id="buni1")
  slice2 <- add_uniform(slice2, id="buni2")
  
  temp1 <- slice1
  temp2 <- slice2
  
  # A critical difference when making slice 1 is that all the new clusters must neighbour the original cluster. This is because we are building slice 2 based on slice 1 (random similarity) and forcing everything to neighbour the original means that there will be no isolated clusters in slice 2 (those that do not neighbour any other cluster)
  
  # For slice 2 some clusters will be guaranteed the same position as the corresponding cluster in slice 1, some will be different. One of these clusters will be randomly similar or different.
  
  while(nlevels(as.factor(slice1$data)) != 5 | nlevels(as.factor(slice2$data)) != 5) {
    slice1 <- temp1
    slice2 <- temp2
    
    slice1 <- add_stripe_or_rect(slice1, id="ashared1", shape=2)
    slice2$data[slice1$data == "ashared1"] <- "ashared1" # adds the same original cluster as in slice 1
    
    
    slice1 <- add_neighbour(slice = slice1, id_originals = c("ashared1"), id_neighbour = "cshared2", stripe = sample(c(TRUE, FALSE), 1)) # adds either a rectangle or a stripe first
    slice2 <- add_neighbour(slice = slice2, id_originals = c("ashared1"), id_neighbour = "cshared2", stripe = sample(c(TRUE, FALSE), 1)) # add neighbour with no spatial consideration, to avoid making both slices ashared1 
    
    
    slice1 <- add_neighbour(slice = slice1, id_originals = c("ashared1"), id_neighbour = "dshared3", stripe = sample(c(TRUE, FALSE), 1)) # does the same thing with the next cluster
    if (sample(c(TRUE, FALSE), 1)) {
      slice2 <- add_neighbour(slice = slice2, id_originals = c("ashared1", "cshared2"), id_neighbour = "dshared3", stripe = sample(c(TRUE, FALSE), 1)) # add neighbour with no spatial consideration
    }
    else {  # add the same cluster that was added to the first 
      slice2$data[slice1$data=="dshared3"] <- "dshared3"
    }
    
    slice1 <- add_neighbour(slice = slice1, id_originals = c("ashared1"), id_neighbour = "eshared4", stripe = sample(c(TRUE, FALSE), 1)) # addition of third cluster
    slice2$data[slice1$data=="eshared4"] <- "eshared4"   # add the same cluster that was added to the first. Otherwise some results might be too different
  }
  
  slices = list(slice1, slice2)
  names(slices) <- c("Slice 1", "Slice 2")
  
  return(slices)
}

#' make_manydiff_somesim_ex12
#' Creates two slices with lots of similar regions that all vary spatially
#' @return A list of two coordinate dataframes
make_manydiff_somesim_ex12 <- function() {
  slice1 <- create_coords_VISIUM_hexes()  # create two blank hex sets to add data to
  slice2 <- create_coords_VISIUM_hexes()
  
  slice1 <- add_uniform(slice1, id="buni1")
  slice2 <- add_uniform(slice2, id="buni2")
  
  for (i in 1:30) {
    slice1 <- add_stripe_or_rect(slice1, id = sample(c("cshared2", "dshared3", "eshared4"), 1), shape=0)
  }
  
  for (i in 1:30) {
    slice2 <- add_stripe_or_rect(slice2, id = sample(c("cshared2", "dshared3", "eshared4"), 1), shape=0)
  }
  
  slice1 <- add_stripe_or_rect(slice1, id="ashared1", shape=2)
  
  slice2 <- add_stripe_or_rect(slice2, id="ashared1", shape=2) 
  
  slices = list(slice1, slice2)
  names(slices) <- c("Slice 1", "Slice 2")
  
  return(slices)
}

#' generate_cluster_types
#' Generates a set of "clusters" by using metadata categories from an sc-RNA seq Seurat object
#' @param sc Seurat single-cell RNA experiment object from which categories will be drawn
#' @param cluster_count The number of clusters to be created, integer
#' @param metadata_id The metadata category in "sc" from which the items will be placed in each cluster, string
#' @param interactive Dataframe containing user-specified items to be placed in the categories. "Clusters" are found on the rows, and any item in the dataframe must exist as a factor level in the metadata category specified. More details about this dataframe are found in create_seurat_pair.
#'
#' @return A list of lists containing the names of the items in each cluster
generate_cluster_types <- function(sc, cluster_count = 3, metadata_id, interactive) {
  
  cell_types <- levels(FetchData(sc, vars = metadata_id)[[1]])
  type_count <- length(cell_types)
  
  clusters <- vector(mode = "list", length = cluster_count)
  
  if (is.data.frame(interactive)) {
    for (i in 1:cluster_count) {
      clusters[[i]] <- as.character(interactive[i, ])
    }
  }
  else {
    for (i in 1:cluster_count)  {
      cluster_to_add <- c(sample(cell_types, min(type_count/cluster_count, 5))) # split up the metadata levels and take some of them for each cluster
      
      cell_types <- cell_types[!(cell_types %in% cluster_to_add)]
      
      clusters[[i]] <- cluster_to_add
    }
  }
  return(clusters)
}

#' create_matrices_for_clusters
#' Subsets the raw count matrix of an sc-RNA Seurat object based on a given metadata category, and organises the subsetted matrices into clusters
#' @param sc Seurat single-cell RNA experiment object
#' @param cluster_count The number of clusters to be created, integer
#' @param metadata_id The metadata category in "sc" from which the clusters will be defined, string
#' @param interactive User-specified cluster data, dataframe
#'
#' @return A list of lists of sparse subsetted count matrices
create_matrices_for_clusters <- function(sc, cluster_count=3, metadata_id, interactive)  {
  Idents(sc) <- metadata_id
  
  clusters <- generate_cluster_types(sc, cluster_count, metadata_id, interactive) 
  
  matrices <- vector(mode = "list", length = cluster_count)
  
  for (i in 1:cluster_count) {
    seurat_matrix <- c()
    names(seurat_matrix) <- c()
    
    for (j in 1:length(clusters[[i]])) {
      seurat_matrix <- c(seurat_matrix, Matrix(GetAssayData(sc, slot="counts")[, WhichCells(sc, ident=clusters[[i]][j])], sparse=TRUE))
    }  
    names(seurat_matrix) <- clusters[[i]]
    matrices[[i]] <- seurat_matrix
  }
  
  names(matrices) <- c(1:length(matrices))
  
  return(matrices)
}

#' add_data_to_sim
#' Takes a simulated spatial pattern and a Seurat object to apply generated clusters of a raw count matrix to the spatial pattern based on the labels of each spot
#' @param spatial A coordinates dataframe
#' @param sc Seurat single-cell RNA experiment object
#' @param metadata_id The metadata category in "sc" from which the clusters will be defined, string
#' @param interactive User-specified cluster data, dataframe
#'
#' @return A raw count matrix barcoded with the same spot names as the spatial data, sparse matrix
add_data_to_sim <- function(spatial, sc, metadata_id, interactive) { 
  counts <- Matrix(data=0, nrow=nrow(sc), ncol=0, sparse=TRUE)
  
  cluster_levels <- levels(as.factor(spatial$data))
  names(cluster_levels) <- c(1:length(cluster_levels))
  
  matrices_for_cell_types <- create_matrices_for_clusters(sc, nlevels(as.factor(spatial$data)), metadata_id, interactive)
  
  for (i in 1:length(names(cluster_levels))) {
    cluster_to_make <- Matrix(data=0, nrow = nrow(sc), ncol=table(spatial$data)[i])  # make a matrix for all cells defined in one cluster type
    for (j in 1:length(matrices_for_cell_types[[i]])) {
      one_cell_matrix <- matrices_for_cell_types[[i]][[j]][, sample(1:ncol(matrices_for_cell_types[[i]][[j]]),table(spatial$data)[i],replace = TRUE)]  # select a random assortment of cells from one of the target cell types who belongs in the cluster to put into all the spots
      one_cell_matrix@Dimnames[[2]] <- rownames(spatial[spatial$data == cluster_levels[i], ])  # track which spots are in which cluster
      cluster_to_make <- cluster_to_make + one_cell_matrix  # gradually add in each cell type's worth of cells to the overall cluster
    }
    if (i == 1) {
      counts <- cluster_to_make
    }
    else {
      counts <- cbind(counts, cluster_to_make)
    }
  }
  return(counts) 
} 

#' create_synth_seurat
#' Takes a simulated spatial pattern and a Seurat object and creates a new Seurat object whose data is based on the labels in each spot of the pattern
#' @param spatial A coordinates dataframe
#' @param sc Seurat single-cell RNA experiment object
#' @param metadata_id The metadata category in "sc" from which the clusters will be defined, string
#' @param interactive User-specified cluster data, dataframe
#'
#' @return A spatial transcriptomics Seurat object
create_synth_seurat <- function(spatial, sc, metadata_id, interactive = 0) {  
  
  counts_synth <- add_data_to_sim(spatial, sc, metadata_id, interactive)  # creates a new count matrix
  
  colnames(spatial) <- c("tissue", "row", "col", "imagerow", "imagecol", "data")
  
  synth_seurat <- CreateSeuratObject(counts=counts_synth, assay="Spatial")  
  
  synth_seurat@images$image =  new(  # adds the data from the spatial grid
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = spatial[,1:5]
  )
  
  return(synth_seurat)
}

#' create_seurat_pair
#' Creates a pair of Seurat spatial transcriptomics objects based on a simulated spatial pattern and single-cell RNA object. 
#' @param pattern_id The type of synthetic spatial pattern to generate, integer from 1-12 inclusive
#' @param sc_data Seurat single-cell RNA experiment object
#' @param metadata_id The metadata category in "sc" from which the clusters will be defined, string
#' @param clusterdata User-specified cluster data, dataframe. The clusters will be stored in the rows. The first row represents the first shared cluster in both slices. The second row represents the unique cluster in slice 1. The third row represents the unique cluster in slice 2. The next three rows represent the next three shared clusters in both slices. Patterns 1 and 2 have ONE cluster per slice (one row in the dataframe). Patterns 3 through 8 inclusive feature TWO clusters per slice (two rows in df). Patterns 9 and 10 feature FOUR clusters per slice (four rows in df). Patterns 11 and 12 feature FIVE clusters per slice (five rows in df)
#'
#' @return A list of two spatial transcriptomics Seurat objects
create_seurat_pair <- function(pattern_id, sc_data, metadata_id, clusterdata = 0) { 
  
  pattern_names <- c("make_idents_ex1",
                     "make_opposites_ex2",
                     "make_simdata_simspat_ex3",
                     "make_diffdata_simspat_ex4",
                     "make_simdata_diffspat_ex5",
                     "make_diffdata_diffspat_ex6",
                     "make_randoms_ex7",
                     "make_somedata_diffloc_diffsize_ex8",
                     "make_somesim_spatunique_ex9",
                     "make_somesim_spatsim_ex10",
                     "make_somesim_somespat_ex11",
                     "make_manydiff_somesim_ex12")
  
  spatial <- get(pattern_names[pattern_id])()
  
  if (!is.data.frame(clusterdata)) {  # randomly generates cluster data in the case that it is not passed in to clusterdata
    if (pattern_id == 1) {
      clusterno <- nlevels(as.factor(spatial[[1]]$data))
    }
    else {
      clusterno <- nlevels(as.factor(spatial[[1]]$data)) + 1
    }
    
    clusters <- generate_cluster_types(sc_data, clusterno, metadata_id, clusterdata)  
    
    clusterdata = data.frame(matrix(NA, nrow = length(clusters), ncol = length(clusters[[1]])))
    
    for (i in (1:clusterno)) {
      clusterdata[i, ] <- clusters[[i]]
    }
  }
  
  seurat_1 <- create_synth_seurat(spatial[[1]], sc_data, metadata_id, clusterdata[c(1,2,4,5,6),])
  
  seurat_2 <- create_synth_seurat(spatial[[2]], sc_data, metadata_id, clusterdata[c(1,3,4,5,6),])
  
  out <- c(seurat_1, seurat_2)
  
  return(out)
}
