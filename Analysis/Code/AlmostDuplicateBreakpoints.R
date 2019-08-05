# almost_duplicate_breakpoints accepts two vectors and a distance threshold
# it then returns indices for which the values in both vectors are closer 
# to each other than the distance threshold allows.
almost_duplicate_breakpoints = function(v1, v2, dist){
  if (length(v1) != length(v2)) stop("v1 and v2 must be of the same length")
  if (dist<0) stop("dist must be a positive number")
  dups = c()
  for (i in 1:(length(v1)-1)){
    simv1 = abs(v1[(i+1):length(v1)] - v1[i]) < dist
    simv2 = abs(v2[(i+1):length(v2)] - v2[i]) < dist
    dups = c(dups, i+which(simv1 & simv2))
  }
  unique(dups)
}
v1 = c(10, 11, 23, 25, 30, 40, 23, 24, 10, 0, 9, 124)
v2 = c(110, 111, 123, 125, 130, 140, 123, 124, 130, 0, 111, 24)
dist = 3
almost_duplicate_breakpoints(v1, v2, dist)

v1 = c()
v2 = c()
dist = 10
almost_duplicate_breakpoints(v1, v2, dist) == c()

v1 = c(NA, 3, 5, 10)
v2 = c(NA, 2, 4, 10)
dist = 10
almost_duplicate_breakpoints(v1, v2, dist) == c()


remove_all_nonunique_elements = function(v){
  v[!(duplicated(v, fromLast=FALSE) | duplicated(v, fromLast=TRUE))]
}

all_almost_nonunique_elements = function(v1, dist){
  if (any(is.na(v1))) stop("input must not contain NaNs")
  if (is.na(dist)) stop("distance must not be NaN")
  if (dist<0) stop("dist must be a positive number")
  dups = list()
  for (i in 1:(length(v1)-1)){
    if (is.na(v1[i])) next
    simv1 = abs(v1[(i+1):length(v1)] - v1[i]) <= dist
    if (any(simv1)) dups[[i]] = c(i, i+which(simv1))
  }
  dups = do.call(c, dups)
  unique(dups[!is.na(dups)])
}
# v = c(110, 111, 123, 125, 130, 140, 123, 124, 130, 0, 111, 24)
# v[-c(all_almost_nonunique_elements(v,5))]
# v = c(110, 111, 123, 125, 130, 140, 123, 124, 130, 0, 111, 24)
# v[-c(all_almost_nonunique_elements(v,0))]
# v = c(10, NA, 20)
# v[-c(all_almost_nonunique_elements(v,0))]

# defines which chromosomes are taken into account, may be used to exclude e.g. chrM
isdefchrom = function(v){
  v == "chr1" | v == "chr2" | v == "chr3" | v == "chr4" | v == "chr5" |  v == "chr6" |  v == "chr7" | v == "chr8" | v == "chr9" | v == "chr10" | v == "chr11" | v == "chr12" | v == "chr13" | v == "chr14" | v == "chr15" | v == "chr16" | v == "chr17" | v == "chr18" | v == "chr19" | v == "chr20" | v == "chr21" | v == "chr22" | v == "chrX" |v == "chrY" 
}
isdefchrom2 = function(v){
  v %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
}


