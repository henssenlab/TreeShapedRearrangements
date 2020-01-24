#python /gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/merge_translocations_diffvariantcallers/script_merge_translocations_smufin_svaba_brass_delly2_novobreak_3.0.py /gpfs/projects/bsc05/elias/Anton_NEUROBLASTOMA/results_smufin_per_patient/results_smufin_large_intraRECcorrected_filterFINAL/somatic_large_SVs.sample_NB2009_intraRECcorrected_filterFINAL.txt /gpfs/projects/bsc05/montse/collaborations/neuroblastoma_pipeline/results_berlin_cohort1/NB2009/svaba/svaba-sv.PASS.vcf.gz /gpfs/projects/bsc05/montse/collaborations/neuroblastoma_pipeline/results_berlin_cohort1/NB2009/brass/NB2009-T1-DNA1-WGS1_vs_NB2009-N1-DNA1-WGS1.annot.vcf.gz /gpfs/projects/bsc05/montse/collaborations/neuroblastoma_pipeline/results_berlin_cohort1/NB2009/delly2/delly2-svs.pre.vcf /gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/novobreak_results/BerlinCohort_novobreak_complete_results.txt NBL2009

#Libraries
import sys
import subprocess
import re
import operator
import statistics
from itertools import chain
import gzip

#Declare lists
results_trans_smufin = []
results_trans_svaba = []
results_trans_brass = []
results_trans_delly2 = []
results_trans_novobreak = []
results_trans_TOTAL = []
results_TRANS_FINAL = []
transloc_R = []
TRANSLOCATION_uniq = []
TRANSLOCATION_same = []
rec_block = []
rec_block_final = []
VCs_used = []
results_TRANS_FINAL_2 = []
results_TRANS_FINAL_3 = []

#Declare arguments
input_file_smufin = sys.argv[1]
input_file_svaba = sys.argv[2]
input_file_brass = sys.argv[3]
input_file_delly2 = sys.argv[4]
input_file_novobreak = sys.argv[5]
patient_id = sys.argv[6]



########## RETRIEVE TRANSLOCATIONS FROM ALL THE CALLERS ##########

#Open and read the different variant callers results files

#Smufin
if input_file_smufin != ".":
 smufin_file = open(input_file_smufin,'r')
 smufin_results = re.split('\n',smufin_file.read())
 for res in range(1,len(smufin_results)-1):
  res_line = re.split('\t',smufin_results[res])
  #transform smufin info into strand info
  aux1_strand = re.split(',',res_line[6])[0]
  aux2_strand = re.split(',',res_line[6])[1]
  if len(re.findall(r'\+',aux1_strand))>0 and len(re.findall(r'\+',aux2_strand))>0:
      strands = '++'
  elif len(re.findall(r'-',aux1_strand))>0 and len(re.findall(r'-',aux2_strand))>0:
      strands = '--'
  elif len(re.findall(r'\+',aux1_strand))>0 and len(re.findall(r'-',aux2_strand))>0:
      strands = '+-'
  elif len(re.findall(r'-',aux1_strand))>0 and len(re.findall(r'\+',aux2_strand))>0:
      strands = '-+'
  else:
      pass
  #Precision SV - smufin allways precise because it uses split-reads
  precision = "PRECISE"
  #select translocations by different chromosome
  if res_line[2] != res_line[4]:
   res_line_trans = list(res_line[2:6])
   res_line_trans.append('SMUFIN=PASS,' + res_line[0]) #info caller
   res_line_trans.append(strands) #info strands
   res_line_trans.append(precision)
   res_line_trans.insert(0, patient_id) #info patient
   results_trans_smufin.append(res_line_trans)
else:
 pass

#Svaba - using the PASS file, all SVs are PASS (we don't have to re-filter)
if input_file_svaba != ".":
 with gzip.open(input_file_svaba, 'rb') as f:
  svaba_file = f.read().decode('utf-8')
  svaba_results = re.split('\n',svaba_file)
  for res1 in range(0,len(svaba_results)-1):
   if len(re.findall(r'#',svaba_results[res1]))>0:
    pass
   else:
    res1_line = re.split('\t',svaba_results[res1])
    #pass vcf to smufin style
    res1_line_chrom2 = re.split('\[|\]',res1_line[4])[1].split(':')[0]
    res1_line_pos2 = re.split('\[|\]',res1_line[4])[1].split(':')[1]
    #get strand info
    if len(re.findall(r'\[',res1_line[4]))==2 and len(re.split('\[',res1_line[4])[0])==0:
     strands = '-+'
    elif len(re.findall(r'\[',res1_line[4]))==2 and len(re.split('\[',res1_line[4])[0])>0:
     strands = '++'
    elif len(re.findall(r'\]',res1_line[4]))==2 and len(re.split('\]',res1_line[4])[0])==0:
     strands = '--'
    elif len(re.findall(r'\]',res1_line[4]))==2 and len(re.split('\]',res1_line[4])[0])>0:
     strands = '+-'
    else:
     pass
    #Precision SV - Svaba indicates if IMPRECISE
    if len(re.findall(r'IMPRECISE',res1_line[7]))>0:
     precision = "IMPRECISE"
    else:
     precision = "PRECISE" 
    #select translocations by different chromosome
    if res1_line[0] != res1_line_chrom2:
     res1_line_trans = list(res1_line[0:2])
     res1_line_trans.append(res1_line_chrom2)
     res1_line_trans.append(res1_line_pos2)
     res1_line_trans.append('SVABA=' + res1_line[6] + ',' + re.split(':',res1_line[2])[0])
     res1_line_trans.append(strands)
     res1_line_trans.append(precision)
     res1_line_trans.insert(0, patient_id)
     results_trans_svaba.append(res1_line_trans)
else:
 pass

#Brass
if input_file_brass != ".":
 with gzip.open(input_file_brass, 'rb') as f1:
  brass_file = f1.read().decode('utf-8')
  brass_results = re.split('\n',brass_file)
  for res2 in range(0,len(brass_results)-1):
   if len(re.findall(r'#',brass_results[res2]))>0:
    pass
   else:
    res2_line = re.split('\t',brass_results[res2])
    #pass vcf to smufin style
    res2_line_chrom2 = re.split('\[|\]',res2_line[4])[1].split(':')[0]
    res2_line_pos2 = re.split('\[|\]',res2_line[4])[1].split(':')[1]
    #get strand info
    if len(re.findall(r'\[',res2_line[4]))==2 and len(re.split('\[',res2_line[4])[0])==0:
     strands = '-+'
    elif len(re.findall(r'\[',res2_line[4]))==2 and len(re.split('\[',res2_line[4])[0])>0:
     strands = '++'
    elif len(re.findall(r'\]',res2_line[4]))==2 and len(re.split('\]',res2_line[4])[0])==0:
     strands = '--'
    elif len(re.findall(r'\]',res2_line[4]))==2 and len(re.split('\]',res2_line[4])[0])>0:
     strands = '+-'
    else:
     pass
    #Precision SV - Svaba indicates if IMPRECISE
    if len(re.findall(r'IMPRECISE',res2_line[7]))>0:
     precision = "IMPRECISE" 
    else:
     precision = "PRECISE"
    info2 = re.split(';',res2_line[7])
    for m2 in range(0,len(info2)):
     if len(re.findall("^BAS=",info2[m2]))>0:
      bas = info2[m2]
      bas_score = int(re.split('=',bas)[1])
      res2_line[7] = bas_score
     else:
      res2_line[7] = "NO_BAS"
      #pass
    #select translocations and filter by bas score>98
     if res2_line[0] != res2_line_chrom2 and isinstance(res2_line[7],int) and res2_line[7]>98: #filter by bas score 99
     #if res2_line[0] != res2_line_chrom2: #line to not filter by BAS score
      res2_line_trans = list(res2_line[0:2])
      res2_line_trans.append(res2_line_chrom2)
      res2_line_trans.append(res2_line_pos2)
      res2_line_trans.append('BRASS=PASS,' + re.split('_',res2_line[2])[0])
      res2_line_trans.append(strands)
      res2_line_trans.append(precision)
      res2_line_trans.insert(0, patient_id)
      results_trans_brass.append(res2_line_trans)
else:
 pass

#Delly2
if input_file_delly2 != ".":
 delly2_file = open(input_file_delly2,'r')
 delly2_results = re.split('\n',delly2_file.read())
 for res3 in range(0,len(delly2_results)-1):
  if len(re.findall(r'#',delly2_results[res3]))>0:
   pass
  else:
   res3_line = re.split('\t',delly2_results[res3])
   #filter for translocations
   if len(re.findall('BND',res3_line[2]))>0:
    #pass vcf to smufin style
    res3_line_chrom2 = re.split('\[|\]',res3_line[4])[1].split(':')[0]
    res3_line_pos2 = re.split('\[|\]',res3_line[4])[1].split(':')[1]
    #get strand info
    if len(re.findall(r'\[',res3_line[4]))==2 and len(re.split('\[',res3_line[4])[0])==0:
     strands = '-+'
    elif len(re.findall(r'\[',res3_line[4]))==2 and len(re.split('\[',res3_line[4])[0])>0:
     strands = '++'
    elif len(re.findall(r'\]',res3_line[4]))==2 and len(re.split('\]',res3_line[4])[0])==0:
     strands = '--'
    elif len(re.findall(r'\]',res3_line[4]))==2 and len(re.split('\]',res3_line[4])[0])>0:
     strands = '+-'
    else:
     pass
    if len(re.findall(r'IMPRECISE',res3_line[7]))>0:
     precision = "IMPRECISE"
    else:
     precision = "PRECISE"
    #select translocations
    if res3_line[0] != res3_line_chrom2 and res3_line[6]=="PASS": #filter for PASS
     res3_line_trans = list(res3_line[0:2])
     res3_line_trans.append(res3_line_chrom2)
     res3_line_trans.append(res3_line_pos2)
     res3_line_trans.append('DELLY2=' + res3_line[6] + ',' + res3_line[2])
     res3_line_trans.append(strands)
     res3_line_trans.append(precision)
     res3_line_trans.insert(0, patient_id)
     results_trans_delly2.append(res3_line_trans)
else:
 pass

#Novobreak - using the only PASS file (don't need to re-filter)
if input_file_novobreak != ".":
 if len(re.findall('NBL', patient_id))==0 and len(re.findall('N', patient_id))>0:
  patient_id_conv = 'C'+re.split('N',patient_id)[1]
 else:
  patient_id_conv = patient_id
 novobreak_file = open(input_file_novobreak,'r')
 novobreak_results = re.split('\n',novobreak_file.read())
 for res4 in range(0,len(novobreak_results)-1):
  #filter for patient (all patients in the same file)
  if len(re.findall(patient_id_conv, novobreak_results[res4]))>0:
   res4_line = re.split('\t',novobreak_results[res4])
   #select translocations
   if res4_line[6] == '<TRA>':
    if res4_line[1] != res4_line[3]:
     res4_line_trans = list(res4_line[1:6])
     res4_line_trans[0] = re.split('chr',res4_line_trans[0])[1]
     res4_line_trans[2] = re.split('chr',res4_line_trans[2])[1]
     res4_line_trans[4] = re.split(';',res4_line_trans[4])[0] + ',' + re.split(';',res4_line_trans[4])[1] #novobreak don't give us the strands 
     res4_line_trans.append('..')
     res4_line_trans.append('..')
     res4_line_trans.insert(0, patient_id)
     results_trans_novobreak.append(res4_line_trans)
else:
 pass

#Remove VCF duplicates from lists (function from https://stackoverflow.com/questions/34630772/removing-duplicate-elements-by-their-attributes-in-python)
#Duplicates are selected by the 5th element, the info which includes the SV id
def deduplicate(items):
 seen = set()
 for item in items:
  if not item[5] in seen:
   seen.add(item[5])
   yield item

#Join all lists and Remove the duplicates using deduplicate function, except smufin and novobreak - we don't have SV ids provided by those callers
if len(results_trans_smufin)>0:
 results_trans_TOTAL.append(results_trans_smufin)

if len(results_trans_svaba)>0:
 results_trans_TOTAL.append(list(deduplicate(results_trans_svaba)))

if len(results_trans_brass)>0:
 results_trans_TOTAL.append(list(deduplicate(results_trans_brass)))

if len(results_trans_delly2)>0:
 results_trans_TOTAL.append(list(deduplicate(results_trans_delly2)))

if len(results_trans_novobreak)>0:
 results_trans_TOTAL.append(results_trans_novobreak)

#Unnest the lists
results_TRANS_FINAL = list(chain.from_iterable(results_trans_TOTAL))
results_TRANS_FINAL = sorted(results_TRANS_FINAL, key=operator.itemgetter(1,3,2,4),  reverse=True)



########## COLLAPSE DUPLICATED TRANSLOCATIONS FROM DIFFERENT/SAME CALLERS ##########

#Change chrX and chrY for 23 and 24 to order them as integers
for T in range(0,len(results_TRANS_FINAL)-1):
 if results_TRANS_FINAL[T][1]=='X' and results_TRANS_FINAL[T][3]=='Y':
  results_TRANS_FINAL[T][1] = '23'
  results_TRANS_FINAL[T][3] = '24'
 elif results_TRANS_FINAL[T][1]=='Y' and results_TRANS_FINAL[T][3]=='X':
  results_TRANS_FINAL[T][1] = '24'
  results_TRANS_FINAL[T][3] = '23'
 elif results_TRANS_FINAL[T][1]=='X':
  results_TRANS_FINAL[T][1] = '23'
 elif results_TRANS_FINAL[T][1]=='Y':
  results_TRANS_FINAL[T][1] = '24'
 elif results_TRANS_FINAL[T][3]=='X':
  results_TRANS_FINAL[T][3] = '23'
 elif results_TRANS_FINAL[T][3]=='Y':
  results_TRANS_FINAL[T][3] = '24'
 else:
  pass

#Order the translocations this way: first lowest chromosome and position ex. 2 67587578 5 13231331 (original: 5 13231331 2 67587578)
for R in range(0,len(results_TRANS_FINAL)):
 #I remove the chromosomes that are not 1:22,X,Y
 if len(re.findall(r'hs37d5',results_TRANS_FINAL[R][1]))>0 or len(re.findall(r'_random',results_TRANS_FINAL[R][1]))>0 or len(re.findall(r'Un_gl',results_TRANS_FINAL[R][1]))>0 or results_TRANS_FINAL[R][1] == 'M' or len(re.findall(r'hs37d5',results_TRANS_FINAL[R][3]))>0 or len(re.findall(r'_random',results_TRANS_FINAL[R][3]))>0 or len(re.findall(r'Un_gl',results_TRANS_FINAL[R][3]))>0 or results_TRANS_FINAL[R][3] == 'M' or len(re.findall(r'GL',results_TRANS_FINAL[R][1]))>0 or len(re.findall(r'GL',results_TRANS_FINAL[R][3]))>0 or len(re.findall(r'MT',results_TRANS_FINAL[R][1]))>0 or len(re.findall(r'MT',results_TRANS_FINAL[R][3]))>0:
  pass
 else:
  results_TRANS_FINAL_2.append(results_TRANS_FINAL[R])

#results_TRANS_FINAL without rare chromosomes
results_TRANS_FINAL = []
results_TRANS_FINAL = results_TRANS_FINAL_2

for R in range(0,len(results_TRANS_FINAL)-1):
 #Change the order of chr-pos
 if int(results_TRANS_FINAL[R][1])>int(results_TRANS_FINAL[R][3]):
  new_chr_1 = results_TRANS_FINAL[R][3]
  new_pos_1 = results_TRANS_FINAL[R][4]
  new_chr_2 = results_TRANS_FINAL[R][1]
  new_pos_2 = results_TRANS_FINAL[R][2]
  if results_TRANS_FINAL[R][6]=='++':
   new_strand = '--'
  elif results_TRANS_FINAL[R][6]=='--':
   new_strand = '++'
  else:
   new_strand = results_TRANS_FINAL[R][6][::-1]
  results_TRANS_FINAL[R][1] = new_chr_1
  results_TRANS_FINAL[R][2] = new_pos_1
  results_TRANS_FINAL[R][3] = new_chr_2
  results_TRANS_FINAL[R][4] = new_pos_2
  results_TRANS_FINAL[R][6] = new_strand
 else:
  pass

#Re-sort list
results_TRANS_FINAL = sorted(results_TRANS_FINAL, key=operator.itemgetter(1,3,2,4),  reverse=False)

#insert strand info on info field
for s in range(0,len(results_TRANS_FINAL)):
  collapse_info_strand = results_TRANS_FINAL[s][5] + ',' + results_TRANS_FINAL[s][6] + ',' + results_TRANS_FINAL[s][7]
  results_TRANS_FINAL[s][5] = collapse_info_strand
  results_TRANS_FINAL_3.append(results_TRANS_FINAL[s])

#results_TRANS_FINAL without rare chromosomes and with strands in info
results_TRANS_FINAL = []
results_TRANS_FINAL = results_TRANS_FINAL_3

#Collapse translocations that have same chr and positions +/-500bp (read length) - both positions
for b in range(0,len(results_TRANS_FINAL)):
 sv = b
 num = sv + 1
 #Save duplicated translocations in TRANSLOCATION_same
 if num < len(results_TRANS_FINAL) and sv < len(results_TRANS_FINAL):
  if results_TRANS_FINAL[sv][3] == results_TRANS_FINAL[num][3] and results_TRANS_FINAL[sv][1] == results_TRANS_FINAL[num][1] and abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))<500 and abs(int(results_TRANS_FINAL[sv][4])-int(results_TRANS_FINAL[num][4]))<500:
   if num < len(results_TRANS_FINAL):
    transloc_R = list(results_TRANS_FINAL[sv])
    TRANSLOCATION_same.append(transloc_R)
   if num == len(results_TRANS_FINAL):
     break
  #Corrects redundancy and the loss of results when we pass from a group of duplicates to a unique
  elif (abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))>500 or abs(int(results_TRANS_FINAL[sv][4])-int(results_TRANS_FINAL[num][4]))>500) and len(transloc_R)>0:
   transloc_R = list(results_TRANS_FINAL[sv])
   TRANSLOCATION_same.append(transloc_R)
   transloc_R = []
  else:
   transloc_U = list(results_TRANS_FINAL[sv])
   TRANSLOCATION_uniq.append(transloc_U) #Save unique translocations in TRANSLOCATION_uniq
  #I do this step to prevent losing the last line on the file
 elif num == len(results_TRANS_FINAL) and len(transloc_R)>0:
  transloc_R = list(results_TRANS_FINAL[len(results_TRANS_FINAL)-1])
  TRANSLOCATION_same.append(transloc_R)
 elif num == len(results_TRANS_FINAL) and len(transloc_R)==0:
  transloc_R = list(results_TRANS_FINAL[len(results_TRANS_FINAL)-1])
  TRANSLOCATION_uniq.append(transloc_R)
 else:
  break

#List rec_block with the coordinates for the duplicate block change (each block contains the duplicates of a specific translocation)
for r in range(0,len(TRANSLOCATION_same)-1):
 if r < len(TRANSLOCATION_same) and r+1 < len(TRANSLOCATION_same):
  if (abs(int(TRANSLOCATION_same[r][2])-int(TRANSLOCATION_same[r+1][2]))>500 or abs(int(TRANSLOCATION_same[r][4])-int(TRANSLOCATION_same[r+1][4]))>500) or (TRANSLOCATION_same[r][1]!=TRANSLOCATION_same[r+1][1] or TRANSLOCATION_same[r][3]!=TRANSLOCATION_same[r+1][3]):
   end_rec_block = r
   rec_block.append(end_rec_block)
  else:
   pass
 else:
  break

rec_block.append(len(TRANSLOCATION_same)-1)
rec_block.insert(0,0)

#Create rec_block_final is the list with all the coord for the blocks
for a in range(0,len(rec_block)-1):
 rec_block_final.append(a)

for rec in range(1,len(rec_block)):
 if rec<len(rec_block):
  rec_block_final[rec-1] = [int(rec_block[rec-1])+1,int(rec_block[rec])]
 else:
  break

rec_block_final[0][0] = 0

list_recurrents_blocks = rec_block_final

#Make final list with the duplicated translocations separated by blocks
for b in range(0,len(rec_block_final)):
 list_recurrents_blocks[b] = TRANSLOCATION_same[list_recurrents_blocks[b][0]:list_recurrents_blocks[b][1]+1]

#Collapse the duplicated TRANSLOCATIONS keeping the coordinates from the first one in the block and appending all info from the different callers
if len(list_recurrents_blocks[0])>0:
 for c in range(0,len(list_recurrents_blocks)):
  for i in range(0,len(list_recurrents_blocks[c])):
   if len(re.findall(r',PRECISE',list_recurrents_blocks[c][i][5]))>0:
    collapse_rec_block = list_recurrents_blocks[c][i]
   elif len(re.findall(r'SMUFIN',list_recurrents_blocks[c][i][5]))>0:
    collapse_rec_block = list_recurrents_blocks[c][i]
   else:
    collapse_rec_block = list_recurrents_blocks[c][0]
   VCs_used.append(list_recurrents_blocks[c][i][5])
  VCs_string = ';'.join(set(VCs_used))
  collapse_rec_block[5] = VCs_string
  TRANSLOCATION_uniq.append(collapse_rec_block) #Save the unique TRANS (from the duplicated blocks) in TRANSLOCATION_uniq - ALL TRANS ARE NOW UNIQUE
  VCs_used = []
else:
 pass


########## SAVE RESULTS IN FILE ##########

#Print callers names for each patient in the file name
if input_file_smufin != ".":
 f_name_smufin = "smufin"
else:
 f_name_smufin = "X"

if input_file_svaba != ".":
 f_name_svaba = "svaba"
else:
 f_name_svaba = "X"

if input_file_brass != ".":
 f_name_brass = "brass"
else:
 f_name_brass = "X"

if input_file_delly2 != ".":
 f_name_delly2 = "delly2"
else:
 f_name_delly2 = "X"

if input_file_novobreak != ".":
 f_name_novobreak = "novobreak"
else:
 f_name_novobreak = "X"

#Open the file to write in
results_patient = open('/gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/merge_translocations_diffvariantcallers/merge_translocations_diffvariantcallers_' + f_name_smufin + '_' + f_name_svaba + '_' + f_name_brass + '_' + f_name_delly2 + '_' + f_name_novobreak + '_' + 'patient_' + patient_id + '_noVCFdups_collapsedups_PASS_BAS99_infostrands_precisioninfo.txt','w')

#Re-sort list
TRANSLOCATION_uniq = sorted(TRANSLOCATION_uniq, key=operator.itemgetter(1,3,2,4),  reverse=False)

#Change chr 23 and 24 for X and Y again
for T in range(0,len(TRANSLOCATION_uniq)):
 if TRANSLOCATION_uniq[T][1]=='23' and TRANSLOCATION_uniq[T][3]=='24':
  TRANSLOCATION_uniq[T][1] = 'X'
  TRANSLOCATION_uniq[T][3] = 'Y'
 elif TRANSLOCATION_uniq[T][1]=='24' and results_TRANS_FINAL[T][3]=='23':
  TRANSLOCATION_uniq[T][1] = 'Y'
  TRANSLOCATION_uniq[T][3] = 'X'
 elif TRANSLOCATION_uniq[T][1]=='23':
  TRANSLOCATION_uniq[T][1] = 'X'
 elif TRANSLOCATION_uniq[T][1]=='24':
  TRANSLOCATION_uniq[T][1] = 'Y'
 elif TRANSLOCATION_uniq[T][3]=='23':
  TRANSLOCATION_uniq[T][3] = 'X'
 elif TRANSLOCATION_uniq[T][3]=='24':
  TRANSLOCATION_uniq[T][3] = 'Y'
 else:
  pass

#Print results
for r in range(0,len(TRANSLOCATION_uniq)):
 results_patient.write('\t'.join(str(v) for v in TRANSLOCATION_uniq[r][0:6]) + '\n')

results_patient.close()

print(patient_id + ' completed')

