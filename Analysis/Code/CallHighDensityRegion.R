call_high_density_region = function (v, kernellength, resol, threshold=0){
  if (length(v) == 0) return(data.frame())
  v = v[!is.na(v)]
  a = c()
  for (i in 1:length(v)){ # for each data coordinate...
    a = c(a, seq(v[i]-(kernellength-1), v[i]+(kernellength-1), resol)) #... enumerate all coordinates (in a low-resolution coordinate system) that are within a distance of windowlength of that point
  }
  a = floor(a/resol)*resol
  df = data.frame(a = unique(a), count = 0)
  for (i in 1:length(a)){
    df[df$a == a[i], "count"] = df[df$a == a[i], "count"] + 1
  }
  df = df[df$count >= threshold,]
  if (nrow(df)==0) return(data.frame())
  df$windowlength = kernellength
  df$resol = resol
  df$countperwindowlength = df$count / df$windowlength
  df = df[order(df$a),]
  df$IntervalID = as.factor(cumsum(c(1,diff(df$a)>resol)))
  
  df = df %>% group_by(IntervalID) %>% summarise(Start = min(as.numeric(a)), End = max(as.numeric(a)))  # change windowlength, ist missverständlich
  
  #df$StartContribs = df$End - (kernellength - 1)
  #df$EndContribs = df$Start + (kernellength - 1)
  
  df$StartContribs = df$Start - (kernellength - 1)
  df$EndContribs = df$End + (kernellength - 1)

  df$n = NA
  df$FirstElement=NA
  df$LastElement=NA
  
  for (i in 1:nrow(df)){
    
    #df[i, "StartContribs"] = min(df[i, "StartContribs"], df[i, "Start"])
    #df[i, "EndContribs"] = max(df[i, "EndContribs"], df[i, "End"])
    
    df[i,"FirstElement"] = min(v[v>=as.numeric(df[i, "StartContribs"]) & v<=as.numeric(df[i, "EndContribs"])])
    df[i,"LastElement"] = max(v[v>=as.numeric(df[i, "StartContribs"]) & v<=as.numeric(df[i, "EndContribs"])])
    df[i,"n"] = length(v[v>=as.numeric(df[i, "StartContribs"]) & v<=as.numeric(df[i, "EndContribs"])])
  }
  
  df$ElementsPerMb = df$n / ((df$LastElement-df$FirstElement)/1000000)
  
  return(df)
}

#call_high_density_region(c(rep(10,20), 40, 24, 40, 25, 24, 28), 1, 1) 
# call_high_density_region(c(rep(10,20), 40, 24, 40, 25, 24, 28), 5, 1, 1)
# call_high_density_region(c(24, 24, 25, 28), 5, 1, threshold = 1)
# call_high_density_region(c(rep(10,20), 40, 24, 40, 25, 24, 28), 5, 1, threshold = 2)
# call_high_density_region(c(6, 10), 5, 1, threshold=1) # there are many intervals of length 5 next to each other that contain at least one data point: [2,6], [3,7], [4,8], ..., [10,14] ---> merge all of them to [2,14]
# call_high_density_region(c(6, 11), 5, 1, threshold=1) # there are many intervals of length 5 next to each other that contain at least one data point: [2,6], [3,7], [4,8], ..., [11,15] ---> merge all of them to [2,15]
# call_high_density_region(c(6, 15), 5, 1, threshold=1) # there are many intervals of length 5 next to each other that contain at least one data point: [2,6], [3,7], [4,8], ..., [11,15] ---> merge all of them to [2,15]
# call_high_density_region(c(6, 16), 5, 1, threshold=1) 
# call_high_density_region(c(6, 10), 5, 1, threshold=2) # there is only one interval of length 5 that contains at least two data points: [6,10] 
# call_high_density_region(c(6, 11), 5, 1, threshold=2) # there is no interval of length 5 that contains at least two data points: [6,10] 

call_high_density_region(c(2,8,9,10,14,16,19),3,1,3)
