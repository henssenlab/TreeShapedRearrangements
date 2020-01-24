#ex. python /gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/merge_translocations_diffvariantcallers/script_filter_mapqBAM.py /gpfs/scratch/bsc05/bsc05096/Anton_NEUROBLASTOMA/merge_translocations_diffvariantcallers/merge_translocations_diffvariantcallers_smufin_svaba_brass_delly2_novobreak_patient_NB2013_noVCFdups_collapsedups_PASS.txt /gpfs/scratch/bsc05/bsc05638/executions_pipeline_PERIS/data_donors/NB2013/alignment/sanger-mem/NB2013-T1-DNA1-WGS1/NB2013-T1-DNA1-WGS1.bam 20 NB2013

#Libraries
import sys
import subprocess
import re
import operator
import statistics
from itertools import chain
import gzip

#Declare lists
FINAL_vcfs_results = []
final_mapq_list = []

#Declare arguments
input_file_vcfs = sys.argv[1]
input_bam = sys.argv[2]
mapq_cutoff = sys.argv[3]
patient_id = sys.argv[4]

########## 

#Open VCFs collapsed results file

file_vcfs = open(input_file_vcfs,'r')
vcfs_results = re.split('\n',file_vcfs.read())
for res in range(0,len(vcfs_results)-1):
 res_line = re.split('\t',vcfs_results[res])
 chr1_bkp = res_line[1]
 pos1_bkp = int(res_line[2])
 chr2_bkp = res_line[3]
 pos2_bkp = int(res_line[4])
 #Get read counts that translocate (go to other chromosomes) within +/-1000bp around the vcf breakpoints in BAM
 command1 = "samtools view -q " + str(mapq_cutoff) + " " + str(input_bam) + " " + str(chr1_bkp) + ":" + str(pos1_bkp-1000) + "-" + str(pos1_bkp+1000) + " | cut -f1-11 | grep -vP "+"\"\\"+"t="+"\\"+"t\" | grep -P "+"\"\\"+"t"+chr2_bkp+"\\"+"t"+str(pos2_bkp)[0:len(str(pos2_bkp))-4]+"...."+"\\"+"t\" | wc -l" #first bkp
 command1_1 = "samtools view " + str(input_bam) + " " + str(chr1_bkp) + ":" + str(pos1_bkp-1000) + "-" + str(pos1_bkp+1000) + " | grep -P "+"\"\\"+"t="+"\\"+"t\" | cut -f17-30" #first bkp SA+XA: multiple mapping + same chrom =
 command1_2 = "samtools view " + str(input_bam) + " " + str(chr1_bkp) + ":" + str(pos1_bkp-1000) + "-" + str(pos1_bkp+1000) + " | grep -vP "+"\"\\"+"t="+"\\"+"t\" | cut -f17-30" #first bkp SA+XA: multiple mapping + diff chrom NO =
 command2 = "samtools view -q " + str(mapq_cutoff) + " " + str(input_bam) + " " + str(chr2_bkp) + ":" + str(pos2_bkp-1000) + "-" + str(pos2_bkp+1000) + " | cut -f1-11 | grep -vP "+"\"\\"+"t="+"\\"+"t\" | grep -P "+"\"\\"+"t"+chr1_bkp+"\\"+"t"+str(pos1_bkp)[0:len(str(pos1_bkp))-4]+"...."+"\\"+"t\" | wc -l" #second bkp
 command2_1 = "samtools view " + str(input_bam) + " " + str(chr2_bkp) + ":" + str(pos2_bkp-1000) + "-" + str(pos2_bkp+1000) + " | grep -P "+"\"\\"+"t="+"\\"+"t\" | cut -f17-30" #second bkp SA+XA: multiple mapping + same chrom =
 command2_2 = "samtools view " + str(input_bam) + " " + str(chr2_bkp) + ":" + str(pos2_bkp-1000) + "-" + str(pos2_bkp+1000) + " | grep -vP "+"\"\\"+"t="+"\\"+"t\" | cut -f17-30" #second bkp SA+XA: multiple mapping + diff chrom NO =
 get_readcount_process1 = subprocess.Popen(command1, stdout=subprocess.PIPE, shell=True)
 get_readcount_process1_1 = subprocess.Popen(command1_1, stdout=subprocess.PIPE, shell=True)
 get_readcount_process1_2 = subprocess.Popen(command1_2, stdout=subprocess.PIPE, shell=True)
 get_readcount_process2 = subprocess.Popen(command2, stdout=subprocess.PIPE, shell=True)
 get_readcount_process2_1 = subprocess.Popen(command2_1, stdout=subprocess.PIPE, shell=True)
 get_readcount_process2_2 = subprocess.Popen(command2_2, stdout=subprocess.PIPE, shell=True)
 out1, err1 = get_readcount_process1.communicate() #out is the output from samtools process - command1
 out1_1, err1_1 = get_readcount_process1_1.communicate() #out is the output from samtools process - command1_1
 out1_2, err1_2 = get_readcount_process1_2.communicate() #out is the output from samtools process - command1_2
 #Do this to get multiple mapping info v1
 out1_1_result = re.split('\n|\t',out1_1.decode('utf-8'))
 out1_2_result = re.split('\n|\t',out1_2.decode('utf-8'))
 for o in range(0,len(out1_1_result)-1):
  if len(re.findall(r'SA:Z',out1_1_result[o]))>0:
   sep_results = re.split('SA:Z:|;',out1_1_result[o])
   list_sep_results = list(filter(None,sep_results))
   for l in list_sep_results:
    pattern = re.compile(chr2_bkp+","+str(pos2_bkp)[0:len(str(pos2_bkp))-4]+"....,")
    if pattern.match(str(l)):
     mapq = re.split(',',l)[4]
     if int(mapq)>=int(mapq_cutoff):
      final_mapq_list.append(mapq)
 out_1_1_count = len(final_mapq_list)
 final_mapq_list = []
 for o in range(0,len(out1_2_result)-1):
  if len(re.findall(r'SA:Z',out1_2_result[o]))>0:
   sep_results = re.split('SA:Z:|;',out1_2_result[o])
   list_sep_results = list(filter(None,sep_results))
   for l in list_sep_results:
    pattern = re.compile(chr2_bkp+","+str(pos2_bkp)[0:len(str(pos2_bkp))-4]+"....,")
    if pattern.match(str(l)):
     mapq = re.split(',',l)[4]
     if int(mapq)>=int(mapq_cutoff):
      final_mapq_list.append(mapq)
 out_1_2_count = len(final_mapq_list)
 final_mapq_list = []
 #Do this to get multiple mapping info v2
 out2, err2 = get_readcount_process2.communicate() #out is the output from samtools process - command2
 out2_1, err2_1 = get_readcount_process2_1.communicate() #out is the output from samtools process - command2_1
 out2_2, err2_2 = get_readcount_process2_2.communicate() #out is the output from samtools process - command2_2
 out_2_1_result = re.split('\n|\t',out2_1.decode('utf-8'))
 out_2_2_result = re.split('\n|\t',out2_2.decode('utf-8'))
 for o in range(0,len(out_2_1_result)-1):
  if len(re.findall(r'SA:Z',out_2_1_result[o]))>0:
   sep_results = re.split('SA:Z:|;',out_2_1_result[o])
   list_sep_results = list(filter(None,sep_results))
   for l in list_sep_results:
    pattern = re.compile(chr2_bkp+","+str(pos2_bkp)[0:len(str(pos2_bkp))-4]+"....,")
    if pattern.match(str(l)):
     mapq = re.split(',',l)[4]
     if int(mapq)>=int(mapq_cutoff):
      final_mapq_list.append(mapq)
 out_2_1_count = len(final_mapq_list)
 final_mapq_list = []
 for o in range(0,len(out_2_2_result)-1):
  if len(re.findall(r'SA:Z',out_2_2_result[o]))>0:
   sep_results = re.split('SA:Z:|;',out_2_2_result[o])
   list_sep_results = list(filter(None,sep_results))
   for l in list_sep_results:
    pattern = re.compile(chr2_bkp+","+str(pos2_bkp)[0:len(str(pos2_bkp))-4]+"....,")
    if pattern.match(str(l)):
     mapq = re.split(',',l)[4]
     if int(mapq)>=int(mapq_cutoff):
      final_mapq_list.append(mapq)
 out_2_2_count = len(final_mapq_list)
 final_mapq_list = []
 read_counts = [int(out1)+int(out_1_1_count)+int(out_1_2_count),int(out2)+int(out_2_1_count)+int(out_2_2_count)]
 #Append flag for read counts and if the translocation PASS or not (LowMAPQ_filter) our MAPQ filtering
 if read_counts[0]>=12 and read_counts[1]>=12:
  res_line.append("MAPQ_filter_PASS=" + str(mapq_cutoff) + ",2bkps," + str(read_counts[0]+read_counts[1]))
 elif (read_counts[0]>=12 and read_counts[1]<12) or (read_counts[0]<12 and read_counts[1]>=12):
  res_line.append("MAPQ_filter_PASS=" + str(mapq_cutoff) + ",1bkp," + str(read_counts[0]+read_counts[1]))
 else:
  res_line.append("LowMAPQ_filter")
 FINAL_vcfs_results.append(res_line)
  

#Open the file to write in
results_patient = open(input_file_vcfs + '.MAPQfiltering' + '.mapqcutoff_' + mapq_cutoff,'w')

#Print results
for r in range(0,len(FINAL_vcfs_results)):
 results_patient.write('\t'.join(str(v) for v in FINAL_vcfs_results[r]) + '\n')

results_patient.close()

print(patient_id + ' completed')


