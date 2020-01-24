########## With this script we can get the intrasv whith bkps falling into circles. It makes a window of all the circles and count how many circles overlap with the bkp ##########

#python /gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/merge_circularized_regions_vs_intraSVs/script_merge_circularized_regions_with_results_merge_VCs_intraSV_2.0_ALL_CIRCLES_resultscirculome_bothbkpsinsamecircle.py /gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/merge_intraSV_diffvariantcallers/results_11.6.19/merge_intraSVs_diffvariantcallers_smufin_svaba_X_delly2_novobreak_patient_NBL34_noVCFdups_collapsedups_nofilters_infostrands.txt.FILTER2CALLERS.svtype NBL34

#Libraries
import sys
import subprocess
import re
import operator
import statistics
from itertools import chain
import numpy as np
import copy
import itertools

#Declare lists
results_circl_regions_FINAL = []
results_circl_VCs_FINAL = []
circular_regions_recurrent = []
intraSV_R = []
INTRASV_uniq = []
group_palmtree = []
FINAL_INTRASV_COLLAPSED = []
results_circl_VCs_FINAL_FILTERPASS = []
circular_regions_recurrent_FILTERPASS = []
INTRASV_uniq_FILTERPASS = []
group_palmtree_FILTERPASS = []
FINAL_INTRASV_COLLAPSED_FILTERPASS = []
list_circls_bkp1 = []
list_circls_bkp2 = []
list_eccDNA_bkp1 = []
list_eccDNA_bkp2 = []
list_ecDNA_bkp1 = []
list_ecDNA_bkp2 = []

circle_identity = 0
intrasv_type = 0
coeff_overlap = 0
coeff_overlap_eccDNA = 0
coeff_overlap_ecDNA = 0
coord_circl = 0

#Declare arguments
input_file_merge_VCs = sys.argv[1]
patient_id = sys.argv[2]

#File with all the amplified regions for all the patients
#input_file_circl_regions = '/gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/circulome_results_analysis/ecDNA_eccDNA_classification_circleseq/ecDNA_eccDNA_unique.bed'
input_file_circl_regions = '/gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/circulome_results_analysis/ecDNA_eccDNA_classification_WGS/ecDNA_eccDNA_WGS_inferredcircles.txt'

########### CREATE LIST OF INTRASVS OVERLAPING CIRCULAR REGIONS ###########

#Open and read the variant callers and amplified regions results files per patient
merge_VCs_file = open(input_file_merge_VCs,'r')
circl_regions_file = open(input_file_circl_regions,'r')

merge_VCs_results = re.split('\n',merge_VCs_file.read())
circl_regions_results = re.split('\n',circl_regions_file.read())

#List with circular regions for your patient
for res in range(0,len(circl_regions_results)-1):
 res_line = re.split('\t',circl_regions_results[res])
 if res_line[3] == patient_id: #Select by patient id
  results_circl_regions_FINAL.append(res_line)

#Make list with all the INTRASVs OVERLAPING the circular regions
#for res2 in range(0,len(results_circl_regions_FINAL)):
for res1 in range(0,len(merge_VCs_results)-1):
 res1_line = re.split('\t',merge_VCs_results[res1])
 for res2 in range(0,len(results_circl_regions_FINAL)):
  #Each bkp inside circles?
  if res1_line[1]==re.split('chr',results_circl_regions_FINAL[res2][0])[1] and int(results_circl_regions_FINAL[res2][1])<int(res1_line[2])<int(results_circl_regions_FINAL[res2][2]):
   if results_circl_regions_FINAL[res2][4]=='eccDNA':
    list_eccDNA_bkp1.append(results_circl_regions_FINAL[res2])
   elif results_circl_regions_FINAL[res2][4]=='ecDNA':
    list_ecDNA_bkp1.append(results_circl_regions_FINAL[res2])
   else:
    pass
   list_circls_bkp1.append(results_circl_regions_FINAL[res2])   
  if res1_line[3]==re.split('chr',results_circl_regions_FINAL[res2][0])[1] and int(results_circl_regions_FINAL[res2][1])<int(res1_line[4])<int(results_circl_regions_FINAL[res2][2]):
   if results_circl_regions_FINAL[res2][4]=='eccDNA':
    list_eccDNA_bkp2.append(results_circl_regions_FINAL[res2])
   elif results_circl_regions_FINAL[res2][4]=='ecDNA':
    list_ecDNA_bkp2.append(results_circl_regions_FINAL[res2])
   else:
    pass
   list_circls_bkp2.append(results_circl_regions_FINAL[res2])
  else:
   pass
 if len(list_circls_bkp1)>0 and len(list_circls_bkp2)>0:
  intrasv_type = "circle-circle" #both bkps on circles
  list_circls_bkp1.sort()
  list_circls_bkp2.sort()
  for l in range(0,len(list_circls_bkp1)):
   for L in range(0,len(list_circls_bkp2)):
    if list_circls_bkp1[l]==list_circls_bkp2[L]:
     intrasv_type = "intracircular"
    else:
     pass
  #Coordinates of the circles window
  chrom_circl1 = max(map(lambda x: x[0], list_circls_bkp1))
  pos1_circl1 = min(map(lambda x: x[1], list_circls_bkp1))
  pos2_circl1 = max(map(lambda x: x[2], list_circls_bkp1))
  chrom_circ2 = max(map(lambda x: x[0], list_circls_bkp2))
  pos1_circl2 = min(map(lambda x: x[1], list_circls_bkp2))
  pos2_circl2 = max(map(lambda x: x[2], list_circls_bkp2))
  coord_circl = chrom_circl1 + ":" + pos1_circl1 + "-" + pos2_circl1 + ";" + chrom_circ2 + ":" + pos1_circl2 + "-" + pos2_circl2
  #Number of circles overlapping with the bkp
  coeff_overlap1 = len(list(list_circls_bkp1 for list_circls_bkp1,_ in itertools.groupby(list_circls_bkp1)))
  coeff_overlap2 = len(list(list_circls_bkp2 for list_circls_bkp2,_ in itertools.groupby(list_circls_bkp2)))
  coeff_overlap = str(coeff_overlap1) + "-" + str(coeff_overlap2)
  coeff_overlap_eccDNA1 = len(list(list_eccDNA_bkp1 for list_eccDNA_bkp1,_ in itertools.groupby(list_eccDNA_bkp1)))
  coeff_overlap_eccDNA2 = len(list(list_eccDNA_bkp2 for list_eccDNA_bkp2,_ in itertools.groupby(list_eccDNA_bkp2)))
  coeff_overlap_ecDNA1 = len(list(list_ecDNA_bkp1 for list_ecDNA_bkp1,_ in itertools.groupby(list_ecDNA_bkp1)))
  coeff_overlap_ecDNA2 = len(list(list_ecDNA_bkp2 for list_ecDNA_bkp2,_ in itertools.groupby(list_ecDNA_bkp2)))
  coeff_overlap_eccDNA = 'eccDNA:' + str(coeff_overlap_eccDNA1) + '-' + str(coeff_overlap_eccDNA2)
  coeff_overlap_ecDNA = 'ecDNA:' + str(coeff_overlap_ecDNA1) + '-' + str(coeff_overlap_ecDNA2)
 elif len(list_circls_bkp1)>0 and len(list_circls_bkp2)==0:
  intrasv_type = "circle-genome" #bkp1 in circles
  list_circls_bkp1.sort()
  #Coordinates of the circles window
  chrom_circl1 = max(map(lambda x: x[0], list_circls_bkp1))
  pos1_circl1 = min(map(lambda x: x[1], list_circls_bkp1))
  pos2_circl1 = max(map(lambda x: x[2], list_circls_bkp1))
  coord_circl = chrom_circl1 + ":" + pos1_circl1 + "-" + pos2_circl1
  #Number of circles overlapping with the bkp
  coeff_overlap1 = len(list(list_circls_bkp1 for list_circls_bkp1,_ in itertools.groupby(list_circls_bkp1)))
  coeff_overlap = str(coeff_overlap1) + "-0"
  coeff_overlap_eccDNA1 = len(list(list_eccDNA_bkp1 for list_eccDNA_bkp1,_ in itertools.groupby(list_eccDNA_bkp1)))
  coeff_overlap_ecDNA1 = len(list(list_ecDNA_bkp1 for list_ecDNA_bkp1,_ in itertools.groupby(list_ecDNA_bkp1)))
  coeff_overlap_eccDNA = 'eccDNA:' + str(coeff_overlap_eccDNA1) + '-0'
  coeff_overlap_ecDNA = 'ecDNA:' + str(coeff_overlap_ecDNA1) + '-0'
 elif len(list_circls_bkp1)==0 and len(list_circls_bkp2)>0:
  intrasv_type = "circle-genome" #bkp2 in circles
  list_circls_bkp2.sort()
  #Coordinates of the circles window
  chrom_circ2 = max(map(lambda x: x[0], list_circls_bkp2))
  pos1_circl2 = min(map(lambda x: x[1], list_circls_bkp2))
  pos2_circl2 = max(map(lambda x: x[2], list_circls_bkp2))
  coord_circl = chrom_circ2 + ":" + pos1_circl2 + "-" + pos2_circl2
  #Number of circles overlapping with the bkp
  coeff_overlap2 = len(list(list_circls_bkp2 for list_circls_bkp2,_ in itertools.groupby(list_circls_bkp2)))
  coeff_overlap = "0-" + str(coeff_overlap2)
  coeff_overlap_eccDNA2 = len(list(list_eccDNA_bkp2 for list_eccDNA_bkp2,_ in itertools.groupby(list_eccDNA_bkp2)))
  coeff_overlap_ecDNA2 = len(list(list_ecDNA_bkp2 for list_ecDNA_bkp2,_ in itertools.groupby(list_ecDNA_bkp2)))
  coeff_overlap_eccDNA = 'eccDNA:0-' + str(coeff_overlap_eccDNA2)
  coeff_overlap_ecDNA = 'ecDNA:0-' + str(coeff_overlap_ecDNA2)  
 else:
  intrasv_type = "genome-genome"
  coord_circl = "None"
  coeff_overlap = "None"
  coeff_overlap_eccDNA = "None"
  coeff_overlap_ecDNA = "None"
 res1_line[2] = int(res1_line[2])
 res1_line[4] = int(res1_line[4])
 res1_line.append(intrasv_type)
 res1_line.append(coord_circl)
 res1_line.append(coeff_overlap)
 res1_line.append(coeff_overlap_eccDNA)
 res1_line.append(coeff_overlap_ecDNA)
 results_circl_VCs_FINAL.append(res1_line)
 list_circls_bkp1 = []
 list_circls_bkp2 = []
 list_eccDNA_bkp1 = []
 list_eccDNA_bkp2 = []
 list_ecDNA_bkp1 = []
 list_ecDNA_bkp2 = []
 intrasv_type = 0
 coeff_overlap = 0
 coeff_overlap_eccDNA = 0
 coeff_overlap_ecDNA = 0 
 coord_circl = 0


#######################

#Open the file to write in
results_patient = open('/gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/merge_circularized_regions_vs_intraSVs/INTRASVS_in_ALL_CIRCLES_WGS_patient_' + patient_id + '_nofilters_infostrands_precisioninfo_FILTER2CALLERS.txt.classification_intratype_circletype','w')

#Print results
for r in range(0,len(results_circl_VCs_FINAL)):
 if results_circl_VCs_FINAL[r][7]!=0:
  results_patient.write('\t'.join(str(v) for v in results_circl_VCs_FINAL[r]) + '\n')

results_patient.close()

print(patient_id + ' completed')

