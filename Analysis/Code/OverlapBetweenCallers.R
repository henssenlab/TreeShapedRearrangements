source("~/Desktop/PalmTrees/Analysis/Code/hasOverlap.R")

load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")

palmtrees$Length = palmtrees$LastElement - palmtrees$FirstElement

samples = unique(palmtrees$Sample)
svcallers = unique(palmtrees$SVCaller)
print(svcallers)
results = data.frame()
for (sample_i in 1:length(samples)){
  print(samples[sample_i])
  for (svcaller_i1 in 1:length(svcallers)){
    for (svcaller_i2 in 1:length(svcallers)){
      
      thispts_caller1 = palmtrees %>% filter(Sample == samples[sample_i], SVCaller == svcallers[svcaller_i1])
      thispts_caller2 = palmtrees %>% filter(Sample == samples[sample_i], SVCaller == svcallers[svcaller_i1])
      
      if (nrow(thispts_caller1) == 0 | nrow(thispts_caller2) == 0){
        if (nrow(results) == 0){
          results = data.frame(
            Sample = samples[sample_i],
            SV1 = svcallers[svcaller_i1],
            SV2 = svcallers[svcaller_i2],
            OverlapPercentSV1 = 0
        }else{
          results = rbind(results,
                          data.frame(
                            Sample = samples[sample_i],
                            SV1 = svcallers[svcaller_i1],
                            SV2 = svcallers[svcaller_i2],
                            OverlapPercentSV1 = 0))
        }
        next
      }
      
      length_overlap = 0
      chrs = unique(thispts_caller1$Chr)
      
      length_overlap = 0 * 1:length(chrs)
      for (chr_i in 1:length(chrs)){
        length_overlap[chr_i] = 
          lengthOverlap(thispts_caller1 %>% filter(Chr == chrs[chr_i]) %>% .$FirstElement,
                        thispts_caller1 %>% filter(Chr == chrs[chr_i]) %>% .$LastElement,
                        thispts_caller2 %>% filter(Chr == chrs[chr_i]) %>% .$FirstElement,
                        thispts_caller2 %>% filter(Chr == chrs[chr_i]) %>% .$LastElement)
      }
      
      if (nrow(results) == 0){
        results = data.frame(
          Sample = samples[sample_i],
          SV1 = svcallers[svcaller_i1],
          SV2 = svcallers[svcaller_i2],
          OverlapPercentSV1 = sum(length_overlap) / sum(thispts_caller1$Length))
      }else{
        results = rbind(results,
                        data.frame(
                          Sample = samples[sample_i],
                          SV1 = svcallers[svcaller_i1],
                          SV2 = svcallers[svcaller_i2],
                          OverlapPercentSV1 = sum(length_overlap) / sum(thispts_caller1$Length)))
      }
      
    }
  }
}

