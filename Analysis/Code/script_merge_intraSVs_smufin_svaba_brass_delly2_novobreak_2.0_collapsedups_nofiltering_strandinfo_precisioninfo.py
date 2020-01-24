#python /gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/merge_intraSV_diffvariantcallers/script_intras.py /gpfs/projects/bsc05/elias/Anton_NEUROBLASTOMA/results_smufin_per_patient/results_smufin_large_intraRECcorrected_filterFINAL/somatic_large_SVs.sample_NB2009_intraRECcorrected_filterFINAL.txt /gpfs/projects/bsc05/montse/collaborations/neuroblastoma_pipeline/results_berlin_cohort1/NB2009/svaba/svaba-sv.PASS.vcf.gz /gpfs/projects/bsc05/montse/collaborations/neuroblastoma_pipeline/results_berlin_cohort1/NB2009/brass/NB2009-T1-DNA1-WGS1_vs_NB2009-N1-DNA1-WGS1.annot.vcf.gz /gpfs/projects/bsc05/montse/collaborations/neuroblastoma_pipeline/results_berlin_cohort1/NB2009/delly2/delly2-svs.pre.vcf /gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/novobreak_results/MERGE_BerlinCohort_PeiferWGS_NovobreakResults_filter_hg19_goodformat.txt NBL2009

#Libraries
import sys
import subprocess
import re
import operator
import statistics
from itertools import chain
import gzip

#Declare lists
results_intra_smufin = []
results_intra_svaba = []
results_intra_brass = []
results_intra_delly2 = []
results_intra_novobreak = []
results_intra_TOTAL = []
results_INTRA_FINAL = []
intra_R = []
INTRASV_uniq = []
INTRASV_same = []
rec_block = []
rec_block_final = []
VCs_used = []
results_INTRA_FINAL_2 = []
results_INTRA_FINAL_3 = []

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
  #select intrachromosomal rearrangements by same chromosome
  if res_line[2] == res_line[4]:
   sv_type = "unknown"
   length = abs(int(res_line[3])-int(res_line[5]))
   res_line_intra = list(res_line[2:6])
   res_line_intra.append('SMUFIN=PASS,' + res_line[0]) #info caller
   res_line_intra.append(strands) #info strands
   res_line_intra.append(precision)
   if strands == '++':
    sv_type = "deletion"
   res_line_intra.append(sv_type)
   res_line_intra.append(length)
   res_line_intra.insert(0, patient_id) #info patient
   results_intra_smufin.append(res_line_intra)
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
    #select intrachromosomal rearrangements by same chromosome
    if res1_line[0] == res1_line_chrom2:
     sv_type = "unknown"
     length = abs(int(res1_line[1])-int(res1_line_pos2))
     res1_line_intra = list(res1_line[0:2])
     res1_line_intra.append(res1_line_chrom2)
     res1_line_intra.append(res1_line_pos2)
     res1_line_intra.append('SVABA=' + res1_line[6] + ',' + re.split(':',res1_line[2])[0])
     res1_line_intra.append(strands)
     res1_line_intra.append(precision)
     if strands == '++':
      sv_type = "deletion"
     res1_line_intra.append(sv_type)
     res1_line_intra.append(length)
     res1_line_intra.insert(0, patient_id)
     results_intra_svaba.append(res1_line_intra)
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
    for m2 in range(0,len(info2)): #retrieve BAS score
     if len(re.findall("^BAS=",info2[m2]))>0:
      bas = info2[m2]
      bas_score = int(re.split('=',bas)[1])
      res2_line[7] = bas_score
     else:
      res2_line[7] = "NO_BAS"
      #pass
     for m3 in range(0,len(info2)): #retrieve SVCLASS
      if len(re.findall("^SVCLASS=",info2[m3]))>0:
       sv_class = info2[m3]
       sv_type = re.split('=',sv_class)[1]
       res2_line[8] = sv_type
      else:
       res2_line[8] = "unknown"
    #select intrachromosomal rearrangements by same chromosome and filter by bas score>98
    #if res2_line[0] == res2_line_chrom2 and isinstance(res2_line[7],int) and res2_line[7]>98: #filter by bas score 99
    if res2_line[0] == res2_line_chrom2: #line to not filter by BAS score
     length = abs(int(res2_line[1])-int(res2_line_pos2))
     res2_line_intra = list(res2_line[0:2])
     res2_line_intra.append(res2_line_chrom2)
     res2_line_intra.append(res2_line_pos2)
     res2_line_intra.append('BRASS=BAS:'+ res2_line[7] + ',' + re.split('_',res2_line[2])[0])
     res2_line_intra.append(strands)
     res2_line_intra.append(precision)
     res2_line_intra.append(res2_line[8])
     res2_line_intra.append(length)
     res2_line_intra.insert(0, patient_id)
     results_intra_brass.append(res2_line_intra)
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
   #pass vcf to smufin style
   res3_line_chrom2 = re.split('CHR2=',res3_line[7])[1].split(';')[0]
   res3_line_pos2 = re.split('END=',res3_line[7])[1].split(';')[0]
   if len(re.findall(r'IMPRECISE',res3_line[7]))>0:
    precision = "IMPRECISE"
   else:
    precision = "PRECISE"
   #select intrachromosomal rearrangements by same chromosome
   #if res3_line[0] == res3_line_chrom2 and res3_line[6]=="PASS": #filter for PASS
   if res3_line[0] == res3_line_chrom2: #don't filter for PASS
    sv_type = "unknown"
    length = abs(int(res3_line[1])-int(res3_line_pos2))
    res3_line_intra = list(res3_line[0:2])
    res3_line_intra.append(res3_line_chrom2)
    res3_line_intra.append(res3_line_pos2)
    res3_line_intra.append('DELLY2=' + res3_line[6] + ',' + res3_line[2])
    res3_line_intra.append('..')
    res3_line_intra.append(precision)
    if len(re.findall(r'DEL', res3_line[2]))>0:
     sv_type = "deletion"
    elif len(re.findall(r'DUP', res3_line[2]))>0:
     sv_type = "tandem-duplication"
    elif len(re.findall(r'INV', res3_line[2]))>0:
     sv_type = "inversion"
    elif len(re.findall(r'INS', res3_line[2]))>0:
     sv_type = "insertion"
    res3_line_intra.append(sv_type)
    res3_line_intra.append(length)
    res3_line_intra.insert(0, patient_id)
    results_intra_delly2.append(res3_line_intra)
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
   #select intrachromosomal rearrangements by same chromosome
   if res4_line[1] == res4_line[3]:
    length = abs(int(res4_line[2])-int(res4_line[4]))
    res4_line_intra = list(res4_line[1:6])
    res4_line_intra[0] = re.split('chr',res4_line_intra[0])[1]
    res4_line_intra[2] = re.split('chr',res4_line_intra[2])[1]
    res4_line_intra[4] = re.split(';',res4_line_intra[4])[0] + ',' + re.split(';',res4_line_intra[4])[1] #novobreak don't give us the strands
    res4_line_intra.append('..')
    res4_line_intra.append('..')
    if res4_line[6] == "<DEL>":
     sv_type = "deletion"
    elif res4_line[6] == "<DUP>":
     sv_type = "tandem-duplication"
    elif res4_line[6] == "<INV>":
     sv_type = "inversion"
    res4_line_intra.append(sv_type)
    res4_line_intra.append(length)
    res4_line_intra.insert(0, patient_id)
    results_intra_novobreak.append(res4_line_intra)
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
if len(results_intra_smufin)>0:
 results_intra_TOTAL.append(results_intra_smufin)

if len(results_intra_svaba)>0:
 results_intra_TOTAL.append(list(deduplicate(results_intra_svaba)))

if len(results_intra_brass)>0:
 results_intra_TOTAL.append(list(deduplicate(results_intra_brass)))

if len(results_intra_delly2)>0:
 results_intra_TOTAL.append(list(deduplicate(results_intra_delly2)))

if len(results_intra_novobreak)>0:
 results_intra_TOTAL.append(results_intra_novobreak)

#Unnest the lists
results_INTRA_FINAL = list(chain.from_iterable(results_intra_TOTAL))
results_INTRA_FINAL = sorted(results_INTRA_FINAL, key=operator.itemgetter(1,2,4),  reverse=True)

########## COLLAPSE DUPLICATED INTRASVS FROM DIFFERENT/SAME CALLERS ##########

#Change chrX and chrY for 23 and 24 to order them as integers
for T in range(0,len(results_INTRA_FINAL)-1):
 if results_INTRA_FINAL[T][1]=='X' and results_INTRA_FINAL[T][3]=='Y':
  results_INTRA_FINAL[T][1] = '23'
  results_INTRA_FINAL[T][3] = '24'
 elif results_INTRA_FINAL[T][1]=='Y' and results_INTRA_FINAL[T][3]=='X':
  results_INTRA_FINAL[T][1] = '24'
  results_INTRA_FINAL[T][3] = '23'
 elif results_INTRA_FINAL[T][1]=='X':
  results_INTRA_FINAL[T][1] = '23'
 elif results_INTRA_FINAL[T][1]=='Y':
  results_INTRA_FINAL[T][1] = '24'
 elif results_INTRA_FINAL[T][3]=='X':
  results_INTRA_FINAL[T][3] = '23'
 elif results_INTRA_FINAL[T][3]=='Y':
  results_INTRA_FINAL[T][3] = '24'
 else:
  pass

#Order the intraSVs this way: first lowest position ex. 2 67587578 2 13231331 (original: 2 13231331 2 67587578)
for R in range(0,len(results_INTRA_FINAL)):
 #I remove the chromosomes that are not 1:22,X,Y and remove results with length<100bp
 if len(re.findall(r'hs37d5',results_INTRA_FINAL[R][1]))>0 or len(re.findall(r'_random',results_INTRA_FINAL[R][1]))>0 or len(re.findall(r'Un_gl',results_INTRA_FINAL[R][1]))>0 or results_INTRA_FINAL[R][1] == 'M' or len(re.findall(r'hs37d5',results_INTRA_FINAL[R][3]))>0 or len(re.findall(r'_random',results_INTRA_FINAL[R][3]))>0 or len(re.findall(r'Un_gl',results_INTRA_FINAL[R][3]))>0 or results_INTRA_FINAL[R][3] == 'M' or len(re.findall(r'GL',results_INTRA_FINAL[R][1]))>0 or len(re.findall(r'GL',results_INTRA_FINAL[R][3]))>0 or len(re.findall(r'MT',results_INTRA_FINAL[R][1]))>0 or len(re.findall(r'MT',results_INTRA_FINAL[R][3]))>0 or int(results_INTRA_FINAL[R][9])<50:
  pass
 else:
  results_INTRA_FINAL_2.append(results_INTRA_FINAL[R])

#results_INTRA_FINAL without rare chromosomes
results_INTRA_FINAL = []
results_INTRA_FINAL = results_INTRA_FINAL_2


for R in range(0,len(results_INTRA_FINAL)-1):
 #Change the order of chr-pos
 if int(results_INTRA_FINAL[R][2])>int(results_INTRA_FINAL[R][4]):
  new_pos_1 = results_INTRA_FINAL[R][4]
  new_pos_2 = results_INTRA_FINAL[R][2]
  if results_INTRA_FINAL[R][6]=='++':
   new_strand = '--'
  elif results_INTRA_FINAL[R][6]=='--':
   new_strand = '++'
  else:
   new_strand = results_INTRA_FINAL[R][6][::-1]
  results_INTRA_FINAL[R][2] = new_pos_1
  results_INTRA_FINAL[R][4] = new_pos_2
  results_INTRA_FINAL[R][6] = new_strand
 else:
  pass

#Re-sort list
results_INTRA_FINAL = sorted(results_INTRA_FINAL, key=operator.itemgetter(1,2,4),  reverse=False)

#insert strand info and svtype on info field
for s in range(0,len(results_INTRA_FINAL)):
  collapse_info_strand = results_INTRA_FINAL[s][5] + ',' + results_INTRA_FINAL[s][6] + ',' + results_INTRA_FINAL[s][7] + ',' + results_INTRA_FINAL[s][8]
  results_INTRA_FINAL[s][5] = collapse_info_strand
  results_INTRA_FINAL_3.append(results_INTRA_FINAL[s][0:6])

#results_INTRA_FINAL without rare chromosomes and with strands in info
results_INTRA_FINAL = []
results_INTRA_FINAL = results_INTRA_FINAL_3

#Collapse intraSVs that have same chr and positions +/-300bp (read length) - both positions
for b in range(0,len(results_INTRA_FINAL)):
 sv = b
 num = sv + 1
 #Save duplicated intraSVs in INTRASV_same
 if num < len(results_INTRA_FINAL) and sv < len(results_INTRA_FINAL):
  if results_INTRA_FINAL[sv][3] == results_INTRA_FINAL[num][3] and results_INTRA_FINAL[sv][1] == results_INTRA_FINAL[num][1] and abs(int(results_INTRA_FINAL[sv][2])-int(results_INTRA_FINAL[num][2]))<=300 and abs(int(results_INTRA_FINAL[sv][4])-int(results_INTRA_FINAL[num][4]))<=300:
   if num < len(results_INTRA_FINAL):
    intra_R = list(results_INTRA_FINAL[sv])
    INTRASV_same.append(intra_R)
   if num == len(results_INTRA_FINAL):
     break
  #Corrects redundancy and the loss of results when we pass from a group of duplicates to a unique
  elif (abs(int(results_INTRA_FINAL[sv][2])-int(results_INTRA_FINAL[num][2]))>300 or abs(int(results_INTRA_FINAL[sv][4])-int(results_INTRA_FINAL[num][4]))>300) and len(intra_R)>0:
   intra_R = list(results_INTRA_FINAL[sv])
   INTRASV_same.append(intra_R)
   intra_R = []
  else:
   transloc_U = list(results_INTRA_FINAL[sv])
   INTRASV_uniq.append(transloc_U) #Save unique intraSV in INTRASV_uniq
  #I do this step to prevent losing the last line on the file
 elif num == len(results_INTRA_FINAL) and len(intra_R)>0:
  intra_R = list(results_INTRA_FINAL[len(results_INTRA_FINAL)-1])
  INTRASV_same.append(intra_R)
 elif num == len(results_INTRA_FINAL) and len(intra_R)==0:
  intra_R = list(results_INTRA_FINAL[len(results_INTRA_FINAL)-1])
  INTRASV_uniq.append(intra_R)
 else:
  break

#List rec_block with the coordinates for the duplicate block change (each block contains the duplicates of a specific intraSV)
for r in range(0,len(INTRASV_same)-1):
 if r < len(INTRASV_same) and r+1 < len(INTRASV_same):
  if (abs(int(INTRASV_same[r][2])-int(INTRASV_same[r+1][2]))>300 or abs(int(INTRASV_same[r][4])-int(INTRASV_same[r+1][4]))>300):
   end_rec_block = r
   rec_block.append(end_rec_block)
  else:
   pass
 else:
  break

rec_block.append(len(INTRASV_same)-1)
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

#Make final list with the duplicated intraSVs separated by blocks
for b in range(0,len(rec_block_final)):
 list_recurrents_blocks[b] = INTRASV_same[list_recurrents_blocks[b][0]:list_recurrents_blocks[b][1]+1]

#Collapse the duplicated INTRASVS keeping the coordinates from the first one in the block and appending all info from the different callers
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
  #VCs_string = ';'.join(set(VCs_used))
  VCs_string = ';'.join(VCs_used) #do not collapse NOVOBREAK info
  collapse_rec_block[5] = VCs_string
  INTRASV_uniq.append(collapse_rec_block) #Save the unique INTRASVS (from the duplicated blocks) in INTRASV_uniq - ALL INTRASVS ARE NOW UNIQUE
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
results_patient = open('/gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/merge_intraSV_diffvariantcallers/merge_intraSVs_diffvariantcallers_' + f_name_smufin + '_' + f_name_svaba + '_' + f_name_brass + '_' + f_name_delly2 + '_' + f_name_novobreak + '_' + 'patient_' + patient_id + '_noVCFdups_collapsedups_nofilters_infostrands_precisioninfo.txt','w')

#Re-sort list
INTRASV_uniq = sorted(INTRASV_uniq, key=operator.itemgetter(1,3,2,4),  reverse=False)

#Change chr 23 and 24 for X and Y again
for T in range(0,len(INTRASV_uniq)):
 if INTRASV_uniq[T][1]=='23' and INTRASV_uniq[T][3]=='24':
  INTRASV_uniq[T][1] = 'X'
  INTRASV_uniq[T][3] = 'Y'
 elif INTRASV_uniq[T][1]=='24' and results_INTRA_FINAL[T][3]=='23':
  INTRASV_uniq[T][1] = 'Y'
  INTRASV_uniq[T][3] = 'X'
 elif INTRASV_uniq[T][1]=='23':
  INTRASV_uniq[T][1] = 'X'
 elif INTRASV_uniq[T][1]=='24':
  INTRASV_uniq[T][1] = 'Y'
 elif INTRASV_uniq[T][3]=='23':
  INTRASV_uniq[T][3] = 'X'
 elif INTRASV_uniq[T][3]=='24':
  INTRASV_uniq[T][3] = 'Y'
 else:
  pass

#Print results
for r in range(0,len(INTRASV_uniq)):
 results_patient.write('\t'.join(str(v) for v in INTRASV_uniq[r][0:6]) + '\n')

results_patient.close()

print(patient_id + ' completed')

