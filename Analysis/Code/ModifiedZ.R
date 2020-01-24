modifiedz = function(v){
  if (mad(v) == 0){
    return(0)
  } else {
    return(median(v) / mad(v))
  }
}