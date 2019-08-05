hasOverlap = function(AStart, AEnd, BStart, BEnd){
  if (min(AEnd - AStart)<0 | min(BEnd-BStart)<0) stop("Start must be smaller than End")
  !(
    ((AStart < BStart) & (AEnd < BStart)) |
      ((BStart < AStart) & (BEnd < AStart))
  )
}
# hasOverlap(10,20,5,8)
# hasOverlap(10,20,21,21)
# hasOverlap(10,20,8,11)
# hasOverlap(10,20,12,15)
# hasOverlap(10,20,18,20)
# hasOverlap(10,20,5,35)
# hasOverlap(5,35,10,20)

hasOverlap_withChr = function(AChr, AStart, AEnd, BChr, BStart, BEnd){
  if (min(AEnd - AStart)<0 | min(BEnd-BStart)<0) stop("Start must be smaller than End")
  !(
    ((AStart < BStart) & (AEnd < BStart)) |
      ((BStart < AStart) & (BEnd < AStart))
  ) &
    (AChr == BChr)
}
#hasOverlap_withChr("chr1", 5,35, "chr1", 10,20)
#hasOverlap_withChr("chr1", 5,35, "chr2", 10,20)
#hasOverlap_withChr("chr1", 5,35, "chr1", 40,50)

hasOverlap_withSampleChr = function(ASample, AChr, AStart, AEnd, BSample, BChr, BStart, BEnd){
  if (min(AEnd - AStart)<0 | min(BEnd-BStart)<0) stop("Start must be smaller than End")
  (!(
    ((AStart < BStart) & (AEnd < BStart)) |
      ((BStart < AStart) & (BEnd < AStart))
  )) &
    (AChr == BChr) &
    (ASample == BSample)
}

lengthOverlap = function(AStart, AEnd, BStart, BEnd){
  if (min(AEnd - AStart)<0 | min(BEnd-BStart)<0) stop("Start must be smaller than End")
  return((pmin(AEnd, BEnd) - pmax(AStart, BStart) + 1)*hasOverlap(AStart, AEnd, BStart, BEnd))
}
# lengthOverlap(10,20,5,8)
# lengthOverlap(10,20,21,21)
# lengthOverlap(10,20,8,11)
# lengthOverlap(10,20,12,15)
# lengthOverlap(10,20,18,20)
# lengthOverlap(10,20,5,35)
# lengthOverlap(5,35,10,20)
# lengthOverlap(5,35,10,8)

FirstIncludesSecond = function(AStart, AEnd, BStart, BEnd){
  if (min(AEnd - AStart)<0 | min(BEnd-BStart)<0) stop("Start must be smaller than End")
  (AStart <= BStart) & (BEnd <= AEnd)
}

FirstIncludesSecond_withChr = function(AChr, AStart, AEnd, BChr, BStart, BEnd){
  if (min(AEnd - AStart)<0 | min(BEnd-BStart)<0) stop("Start must be smaller than End")
  (AChr == BChr) & (AStart <= BStart) & (BEnd <= AEnd)
}
