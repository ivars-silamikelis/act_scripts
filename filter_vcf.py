import vcf
import numpy as np
import pandas as pd
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
input_name="act.raw.ug.7.vcf" #"test.vcf" #"act.raw.ug.5.vcf",
ad_threshold = 0.6
dp_threshold = 10
vcf_reader= vcf.Reader(open(input_name,"r"))
for sampl in vcf_reader:
	break
vcf_reader= vcf.Reader(open(input_name,"r"))

sequences=dict()
seqs_string=dict()
dps = dict()
dps_stats=dict()
for sample in sampl:
	sam=re.sub(".variant\d*","",sample.sample)
	sequences[sam]=[]
        dps[sam]=[]
        dps_stats[sam]=dict()
is_ok=1
j=1
# getRefIdx :: [Vcf_record] -> Int
def write_vcf(name1,name2,good_vcf):
    writer=vcf.Writer(open(name1,"w"), vcf.Reader(open(name2,"r")))
    for i in range(0, len(good_vcf)):
        writer.write_record(good_vcf[i])
    writer.close()
good_vcfs = []
def get_ref_idx (vcf_record):
	for i in zip (range(0,len (vcf_record.samples)), vcf_record.samples):
		if i[1].sample == "TCDC-AB0715":
			return i[0]

# check_for_artifacts :: [Vcf_record] -> Int -> Bool
def check_for_artifacts(vcf_record, idx):
	if vcf_record.samples[idx]["GT"]==None or vcf_record.samples[idx]["GT"]!='0':
		return True
	else:
		return False



#def check_constraints(vcf_record, dp, ad, gt_idx):
#	if gtype_idx > 0:
		
counter=1		
ref_idx = get_ref_idx(sampl)
for rec in vcf_reader:
	for sample_rec in rec.samples:
            if sample_rec["DP"]==None:
                break
            else:
                sam=re.sub(".variant\d*","",sample_rec.sample)
                dps[sam].append(sample_rec["DP"])
	print counter
	counter+=1

for cov_sample in dps.keys():
        dps_stats[cov_sample]["median"]=np.median([dps[cov_sample]])
        dps_stats[cov_sample]["std"]=np.std([dps[cov_sample]])


vcf_reader= vcf.Reader(open(input_name,"r"))

i=1
for rec in vcf_reader:
	#print i
	i+=1
#        print i
	#nonem artifaktus
	if rec.samples[ref_idx]["GT"]==None or rec.samples[ref_idx]["GT"]!='0':
		continue
	if (check_for_artifacts(rec, ref_idx)):
		continue
	is_ok=1
	for sample_rec in rec.samples:
		if sample_rec["GT"]==None:
			is_ok=0
			break
		gtype_idx=int(sample_rec["GT"])
                #print sample_rec["AD"][gtype_idx],sample_rec["AD"][0]
                sam=re.sub(".variant\d*","",sample_rec.sample)
                #dp filter or fraction filter # or std filter
		if sample_rec["DP"]<dp_threshold or sample_rec["AD"][gtype_idx]/float(np.sum(sample_rec["AD"]))<ad_threshold: # or ((np.abs(dps_stats[sam]["median"] - sample_rec["DP"])/dps_stats[sam]["std"] > 2)):
                        #print ((np.abs(dps_stats[sam]["median"] - sample_rec["DP"])/dps_stats[sam]["std"]))
                        is_ok=0
			break
			
	if is_ok==0:
		continue
	if is_ok==1:
		good_vcfs.append(rec)
	for sample_rec in rec.samples:
		sam=re.sub(".variant\d*","",sample_rec.sample)
		if sample_rec["GT"]=="0" and is_ok==1:
			sequences[sam].append(rec.REF)
		else:
			if is_ok==1:
				gt=int(sample_rec["GT"])-1
				sequences[sam].append(str(rec.ALT[gt]))
for sample in sequences.keys():
	seqs_string[sample]="".join(sequences[sample])
records=[]
write_vcf("output_all_f"+str(ad_threshold)+"_c"+str(dp_threshold)+".7.vcf", input_name, good_vcfs)
#for sample in seqs_string:
#	records.append(SeqRecord(Seq(seqs_string[sample]), id=sample, description=""))
#SeqIO.write(records, "output_all_f"+str(ad_threshold)+"_c"+str(dp_threshold)+".7.fasta", "fasta")
