#python /gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/merge_circularized_or_amplified_regions_vs_translocations_VCs/script_merge_circularized_regions_with_results_merge_VCs_new_5.0_ALL_CIRCLES_circleseq_CLASSIFICATION_TRANS_TYPE_infostrands_infocircletype.py /gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/merge_translocations_diffvariantcallers/results_22.5.19_strandinfo_finalres/merge_translocations_diffvariantcallers_smufin_svaba_brass_delly2_novobreak_patient_NB2013_noVCFdups_collapsedups_PASS_BAS99_infostrands.txt.MAPQfiltering.mapqcutoff_20 NB2013

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
results_amp_regions_FINAL = []
results_amp_VCs_FINAL = []
amplified_regions_recurrent = []
transloc_R = []
TRANSLOCATION_uniq = []
group_palmtree = []
FINAL_TRANSLOCATION_COLLAPSED = []
results_amp_VCs_FINAL_FILTERPASS = []
amplified_regions_recurrent_FILTERPASS = []
TRANSLOCATION_uniq_FILTERPASS = []
group_palmtree_FILTERPASS = []
FINAL_TRANSLOCATION_COLLAPSED_FILTERPASS = []
list_amps_bkp1 = []
list_amps_bkp2 = []
list_eccDNA_bkp1 = []
list_eccDNA_bkp2 = []
list_ecDNA_bkp1 = []
list_ecDNA_bkp2 = []

trans_type = 0
coeff_overlap = 0
coeff_overlap_eccDNA = 0
coeff_overlap_ecDNA = 0
coord_amp = 0

#Declare arguments
input_file_merge_VCs = sys.argv[1]
patient_id = sys.argv[2]

#File with all the amplified regions for all the patients
input_file_amp_regions = '/gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/circulome_results_analysis/ecDNA_eccDNA_classification_circleseq/ecDNA_eccDNA_unique.bed'
#input_file_amp_regions = '/gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/circulome_results_analysis/ecDNA_eccDNA_classification_WGS/ecDNA_eccDNA_WGS_inferredcircles.txt'

########### CREATE LIST OF TRANSLOCATIONS OVERLAPING AMPLIFIED REGIONS ###########

#Open and read the variant callers and amplified regions results files per patient
merge_VCs_file = open(input_file_merge_VCs,'r')
amp_regions_file = open(input_file_amp_regions,'r')

merge_VCs_results = re.split('\n',merge_VCs_file.read())
amp_regions_results = re.split('\n',amp_regions_file.read())

#List with amplified regions for your patient
for res in range(0,len(amp_regions_results)-1):
 res_line = re.split('\t',amp_regions_results[res])
 if res_line[3] == patient_id: #Select by patient id
  results_amp_regions_FINAL.append(res_line)

#Make list with all the translocations OVERLAPING the amplified regions
#for res2 in range(0,len(results_amp_regions_FINAL)):
for res1 in range(0,len(merge_VCs_results)-1):
 res1_line = re.split('\t',merge_VCs_results[res1])
 for res2 in range(0,len(results_amp_regions_FINAL)):
  if res1_line[1]==re.split('chr',results_amp_regions_FINAL[res2][0])[1] and int(results_amp_regions_FINAL[res2][1])<int(res1_line[2])<int(results_amp_regions_FINAL[res2][2]):
   if results_amp_regions_FINAL[res2][4]=='eccDNA':
    list_eccDNA_bkp1.append(results_amp_regions_FINAL[res2])
   elif results_amp_regions_FINAL[res2][4]=='ecDNA':
    list_ecDNA_bkp1.append(results_amp_regions_FINAL[res2])
   else:
    pass
   list_amps_bkp1.append(results_amp_regions_FINAL[res2])
  if res1_line[3]==re.split('chr',results_amp_regions_FINAL[res2][0])[1] and int(results_amp_regions_FINAL[res2][1])<int(res1_line[4])<int(results_amp_regions_FINAL[res2][2]):
   if results_amp_regions_FINAL[res2][4]=='eccDNA':
    list_eccDNA_bkp2.append(results_amp_regions_FINAL[res2])
   elif results_amp_regions_FINAL[res2][4]=='ecDNA':
    list_ecDNA_bkp2.append(results_amp_regions_FINAL[res2])
   else:
    pass    
   list_amps_bkp2.append(results_amp_regions_FINAL[res2])
  else:
   pass
 if len(list_amps_bkp1)>0 and len(list_amps_bkp2)>0:
  trans_type = "circle-circle" #2bkps on circles
  list_amps_bkp1.sort()
  list_amps_bkp2.sort()
  #Coordinates of the circles window
  chrom_amp1 = max(map(lambda x: x[0], list_amps_bkp1))
  pos1_amp1 = min(map(lambda x: x[1], list_amps_bkp1))
  pos2_amp1 = max(map(lambda x: x[2], list_amps_bkp1))
  chrom_amp2 = max(map(lambda x: x[0], list_amps_bkp2))
  pos1_amp2 = min(map(lambda x: x[1], list_amps_bkp2))
  pos2_amp2 = max(map(lambda x: x[2], list_amps_bkp2))
  coord_amp = chrom_amp1 + ":" + pos1_amp1 + "-" + pos2_amp1 + ";" + chrom_amp2 + ":" + pos1_amp2 + "-" + pos2_amp2
  #Number of circles overlapping with the bkp
  coeff_overlap1 = len(list(list_amps_bkp1 for list_amps_bkp1,_ in itertools.groupby(list_amps_bkp1)))
  coeff_overlap2 = len(list(list_amps_bkp2 for list_amps_bkp2,_ in itertools.groupby(list_amps_bkp2)))
  coeff_overlap = str(coeff_overlap1) + "-" + str(coeff_overlap2)
  coeff_overlap_eccDNA1 = len(list(list_eccDNA_bkp1 for list_eccDNA_bkp1,_ in itertools.groupby(list_eccDNA_bkp1)))
  coeff_overlap_eccDNA2 = len(list(list_eccDNA_bkp2 for list_eccDNA_bkp2,_ in itertools.groupby(list_eccDNA_bkp2)))
  coeff_overlap_ecDNA1 = len(list(list_ecDNA_bkp1 for list_ecDNA_bkp1,_ in itertools.groupby(list_ecDNA_bkp1)))
  coeff_overlap_ecDNA2 = len(list(list_ecDNA_bkp2 for list_ecDNA_bkp2,_ in itertools.groupby(list_ecDNA_bkp2)))
  coeff_overlap_eccDNA = 'eccDNA:' + str(coeff_overlap_eccDNA1) + '-' + str(coeff_overlap_eccDNA2)
  coeff_overlap_ecDNA = 'ecDNA:' + str(coeff_overlap_ecDNA1) + '-' + str(coeff_overlap_ecDNA2)
 elif len(list_amps_bkp1)>0 and len(list_amps_bkp2)==0:
  trans_type = "circle-genome" #1bkp on circle, 1bkp on genome
  list_amps_bkp1.sort()
  #Coordinates of the circles window
  chrom_amp1 = max(map(lambda x: x[0], list_amps_bkp1))
  pos1_amp1 = min(map(lambda x: x[1], list_amps_bkp1))
  pos2_amp1 = max(map(lambda x: x[2], list_amps_bkp1))
  coord_amp = chrom_amp1 + ":" + pos1_amp1 + "-" + pos2_amp1
  #Number of circles overlapping with the bkp
  coeff_overlap1 = len(list(list_amps_bkp1 for list_amps_bkp1,_ in itertools.groupby(list_amps_bkp1)))
  coeff_overlap = str(coeff_overlap1) + "-0"
  coeff_overlap_eccDNA1 = len(list(list_eccDNA_bkp1 for list_eccDNA_bkp1,_ in itertools.groupby(list_eccDNA_bkp1)))
  coeff_overlap_ecDNA1 = len(list(list_ecDNA_bkp1 for list_ecDNA_bkp1,_ in itertools.groupby(list_ecDNA_bkp1)))
  coeff_overlap_eccDNA = 'eccDNA:' + str(coeff_overlap_eccDNA1) + '-0'
  coeff_overlap_ecDNA = 'ecDNA:' + str(coeff_overlap_ecDNA1) + '-0'
 elif len(list_amps_bkp1)==0 and len(list_amps_bkp2)>0:
  trans_type = "circle-genome" #1bkp on circle, 1bkp on genome
  list_amps_bkp2.sort()
  #Coordinates of the circles window
  chrom_amp2 = max(map(lambda x: x[0], list_amps_bkp2))
  pos1_amp2 = min(map(lambda x: x[1], list_amps_bkp2))
  pos2_amp2 = max(map(lambda x: x[2], list_amps_bkp2))
  coord_amp = chrom_amp2 + ":" + pos1_amp2 + "-" + pos2_amp2
  #Number of circles overlapping with the bkp
  coeff_overlap2 = len(list(list_amps_bkp2 for list_amps_bkp2,_ in itertools.groupby(list_amps_bkp2)))
  coeff_overlap = str(coeff_overlap2) + "-0"
  coeff_overlap_eccDNA2 = len(list(list_eccDNA_bkp2 for list_eccDNA_bkp2,_ in itertools.groupby(list_eccDNA_bkp2)))
  coeff_overlap_ecDNA2 = len(list(list_ecDNA_bkp2 for list_ecDNA_bkp2,_ in itertools.groupby(list_ecDNA_bkp2)))
  coeff_overlap_eccDNA = 'eccDNA:' + str(coeff_overlap_eccDNA2) + '-0'
  coeff_overlap_ecDNA = 'ecDNA:' + str(coeff_overlap_ecDNA2) + '-0'
  idx = [0,3,4,1,2,5,6]
  res1_line = [res1_line[i] for i in idx]
 else:
  trans_type = "genome-genome"
  coord_amp = "None"
  coeff_overlap = "None"
  coeff_overlap_eccDNA = "None"
  coeff_overlap_ecDNA = "None"
 res1_line[2] = int(res1_line[2])
 res1_line[4] = int(res1_line[4])
 res1_line.append(trans_type)
 res1_line.append(coord_amp)
 res1_line.append(coeff_overlap)
 res1_line.append(coeff_overlap_eccDNA)
 res1_line.append(coeff_overlap_ecDNA)
 results_amp_VCs_FINAL.append(res1_line)
 list_amps_bkp1 = []
 list_amps_bkp2 = []
 list_eccDNA_bkp1 = []
 list_eccDNA_bkp2 = []
 list_ecDNA_bkp1 = []
 list_ecDNA_bkp2 = []
 trans_type = 0
 coeff_overlap = 0
 coeff_overlap_eccDNA = 0
 coeff_overlap_ecDNA = 0
 coord_amp = 0


########### CREATE TWO GROUPS OF RESULTS - ALL TRANSLOCATIONS AND FILTERED TRANSLOCATIONS (exclude LowMAPQ_filter ones) ###########

for trans in range(0,len(results_amp_VCs_FINAL)):
 if results_amp_VCs_FINAL[trans][6]!='LowMAPQ_filter':
  results_amp_VCs_FINAL_FILTERPASS.append(results_amp_VCs_FINAL[trans])


#######################

#Open the file to write in
results_patient = open('/gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/merge_circularized_or_amplified_regions_vs_translocations_VCs/TRANSLOCATIONS_in_ALL_CIRCLES_circleseq_circletype_infostrands_patient_' + patient_id + '_filter_PASS.txt','w')

results_patient_FILTERPASS = open('/gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/merge_circularized_or_amplified_regions_vs_translocations_VCs/TRANSLOCATIONS_in_ALL_CIRCLES_circleseq_circletype_infostrands_patient_' + patient_id + '_filter_PASS+MAPQ60.txt','w')

#Print results
for r in range(0,len(results_amp_VCs_FINAL)):
 if results_amp_VCs_FINAL[r][7]!=0:
  results_patient.write('\t'.join(str(v) for v in results_amp_VCs_FINAL[r]) + '\n')

for r in range(0,len(results_amp_VCs_FINAL_FILTERPASS)):
 if results_amp_VCs_FINAL_FILTERPASS[r][7]!=0:
  results_patient_FILTERPASS.write('\t'.join(str(v) for v in results_amp_VCs_FINAL_FILTERPASS[r]) + '\n')

results_patient.close()
results_patient_FILTERPASS.close()

print(patient_id + ' completed')

