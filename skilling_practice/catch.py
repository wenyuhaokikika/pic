import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

step=0#每次增加的步长
resultfile='result.fa'#结果文件路径
gfffile='ASM14920v2-_outfmt6-perl.gff'#blast后perl的转化gff文件
fnafile='GCF_000149205.2_ASM14920v2_genomic.fna'#物种的基因组

def catch(step,resultfile,gfffile,fnafile):
	data=pd.read_csv(gfffile,sep="\t",skiprows=1,header=None)
	data=data.loc[data[2]=='gene',:]
	data[8]=data[8]+";chrom="+data[0]+";start="+data[3].astype(str)+';end='+data[4].astype(str)+';strand='+data[6]
	with open(resultfile,'a+') as sf:
		for seq_record in SeqIO.parse(fnafile, "fasta"):
			da=data.loc[data[0]==seq_record.id,:]
			da[3],da[4]=da[3]-step,da[4]+step
			for start,end,des in zip(da[3],da[4],data[8]):
				seq=seq_record.seq[start:end+1]
				re=SeqRecord(Seq(str(seq)),id=des,description='')
				SeqIO.write(re, sf, "fasta")
if __name__=='__main__':
	if(os.path.exists(resultfile)):
		os.remove(resultfile)
	catch(step,resultfile,gfffile,fnafile)