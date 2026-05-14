import argparse
parser=argparse.ArgumentParser( description='Program calculates bacterial chromosome coordinates from given .gff file')
parser.add_argument("gff",help="Name of input .gff file downloaded from GenBank ")
parser.add_argument("-o",help="output file")
parser.parse_args()
args = parser.parse_args()
gff=args.gff
if args.o != None:
	bcc=args.o
else:
	bcc="bcc.bcc"
def sgrp(gl,ori,dnaAsign,lct):
# returns gene coordinates
# gl chromosome length
# ori centre of ori
# +/- of dnaA gene from .gff file
# columns 3,4 and 6 from the record in the .gff file
	strand={True:"D",False:"G"}
	g=sum([int(x) for x in lct[:2]])/2
	if g>ori:
		gori=g-ori
	else:
		gori=g-ori+gl
	if abs(gori)>gl/2:
		gori-=gl/2
	if abs(gori)>gl:
		gori+=gl
	dnaAstrand=dnaAsign=="+"
	gstrand=lct[2]=="+"
	if dnaAstrand:
		if ori>g:
			gori=gl-ori+g
		else:
			gori=g-ori
	else:
		gstrand=not gstrand
		if ori>g:
			gori=ori-g
		else:
			gori=ori+gl-g
	if gori>gl/2:
		gori=gori-gl
		gstrand= not gstrand
	pom=gori/gl
	val=round(gori/gl*200,3)
	return(str(val)+strand[gstrand])
def orientation(gff):
#  function returns dictionary {"error":error,"ori":ori,"gl":gl,"dnaAsign":dnaAsign,"genes":[(gene,start,end,strand,locus_tag)]} if error is "--" ori is centre of ori, gl is length of the chromosome and dnaAsign is strand of the dnaA in the sequence, genes are gff coordinates of genes from the chromosome 
	error="--"
	r=open(gff,"r")
	genes=[]
	region="none"
	gl=0
	chromosome="*"
	ori=0
	dnaAsign="--"
# reading gff coordinates of CDSs from the chromosome
	for line in r:
		if line[0]!="#":
			cols=[x.strip().strip("\"") for x in line.split("\t")]
			cigar={y[0]:y[1] for y in [x.split("=") for x in cols[8].split(";")]}
			gn="--"
			if "gene" in cigar:
				gn=cigar["gene"]
			if cols[2]=="region":
				if "genome" in cigar:
					region=cigar["genome"]
					if cigar["genome"]=="genomic":
						region="chromosome"
				elif "chromosome" in cigar:
					if cigar["chromosome"]=="1":
						region="chromosome"
					else:
						region="other"
				else:
					region="none"
				if region=="chromosome":
					if chromosome=="*" and gl==0:
						gl=int(cols[4])-int(cols[3])
						chromosome=cols[0]
					else:
						error="more than one chromosome"
			if cols[2]=="gene" and region=="chromosome":
				genes.append((gn,int(cols[3]),int(cols[4]),cols[6],cigar["locus_tag"]))
# checking strand of genes surrouding dnaA
	dnaA=[genes.index(x) for x in genes if x[0][:4]=="dnaA" ]
	if error!="--":
		pass
	elif len(genes)==0:
		error="no region assigned as chromosome"
	elif len(dnaA)==0:
		error="dnaA not found in the chromosome "
	elif len(dnaA)>1:
		error="more than one dnaA gene in the chromosome"
	else: 
		dnaAsign=genes[dnaA[0]][3]
		rng=6
		#for i in range(1,dnaA[0]-rng):
		if genes[dnaA[0]][3]=="+":
		    before=[genes[dnaA[0]-i] for i in range(rng,0,-1) if dnaA[0]-i>=0]+[genes[-i] for i in range(rng,0,-1) if dnaA[0]-i<0]   
		    after=[genes[dnaA[0]+i] for i in range(0,rng) if dnaA[0]+i<len(genes)]+[genes[i] for i in range(0,rng) if dnaA[0]+i>len(genes)]  
		else:
		    after=[genes[dnaA[0]-i] for i in range(rng,-1,-1) if dnaA[0]-i>=0]+[genes[-i] for i in range(rng,0,-1) if dnaA[0]-i<0]   
		    before=[genes[dnaA[0]+i] for i in range(1,rng) if dnaA[0]+i<len(genes)]+[genes[i] for i in range(1,rng) if dnaA[0]+i>len(genes)]  
		before_leading=len([x for x in before if x[3]==genes[dnaA[0]][3]])
		after_leading=len([x for x in after if x[3]==genes[dnaA[0]][3]])
		if before_leading<=rng/2 and after_leading>=rng/2:
#			print("***************************")	
			if genes[dnaA[0]][3]=="+":
				begin=before[-1][2]
				end=after[0][1]    
				if begin<end:
					ori=begin+int((end-begin)/2)-1
				else:
					ori=begin+int((gl-begin+end)/2)-1
				if ori>gl:
					ori-=gl
			else:
				begin=before[0][1]
				end=after[-1][2]
				if begin>end:
					ori=begin+int((end-begin)/2)-1
				else:
					ori=begin+int((gl-begin+end)/2)-1
				if ori>gl:
					ori-=gl
		else:
			error="strands don't change near dnaA"
	rt={"error":error}
	rt["dnaAsign"]=dnaAsign
	rt["gl"]=gl
	rt["ori"]=ori
	rt["genes"]=[x for x in genes]
#	todo={x:[y for y in genes if y[0].split("_")[0]==x] for x in genes}
#	for td in todo:
#		rt[td]=[sgrp(gl,ori,dnaAsign,y[1:4]) for y in todo[td]] 		
	return rt
ori=orientation(gff)
if ori["error"]!="--":
	exit("BCC not computed: "+ori["error"])
w=open(bcc,"w")
for gene in ori["genes"]:
	w.write("\t".join([gene[4],str(sgrp(ori["gl"],ori["ori"],ori["dnaAsign"],gene[1:4]))])+"\n")
print("Bacterial chromosome coordinates were written to the file "+bcc)
