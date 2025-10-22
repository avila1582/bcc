import argparse
parser=argparse.ArgumentParser( description='Program calculates bacterial chromosome coordinates from given .gff file')
parser.add_argument("gff",help="Name of input .gff file downloaded from GenBank or created by prokka")
parser.add_argument("ori",help="ori as set by orifinder (https://tubic.org/Ori-Finder2022/) e.g \"4,647,931 ... 204\"")
parser.add_argument("dnaA",help="locus tag of dnaA gene in the .gff file")
parser.add_argument("-o",help="output file")
parser.parse_args()
args = parser.parse_args()
gff=args.gff
ori=args.ori
dnaA=args.dnaA
if args.o != None:
	bcc=args.o
else:
	bcc="bcc.bcc"
def sgrp(gl,ori,dnaAsign,gi):
# returns gene coordinates
# gl genome length
# ori centre of ori
# +/- of dnaA gene from .gff file
# gi columns 3,4,5  and 6 from .gff file
	strand={True:"D",False:"G"}
	g=sum([int(x) for x in gi[:2]])/2
	if g>ori:
		gori=g-ori
	else:
		gori=g-ori+gl
	if abs(gori)>gl/2:
		gori-=gl/2
	if abs(gori)>gl:
		gori+=gl
	dnaAstrand=dnaAsign=="+"
	gstrand=gi[3]=="+"
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
	val=round(gori/gl*200,2)
	return(str(val)+strand[gstrand])
# Reading regions and chromosome finding
r=open(gff,"r")
regions={}
for line in r:
	if line.find("region") >-1 and line[0]=="#":
		cols=[x.strip().strip("\"") for x in line.split(" ")]
		regions[cols[1]]=int(cols[3])
chrom=[k for k,v in sorted(regions.items(),key=lambda x:x[1],reverse=True)][0]
gl=regions[chrom]
# parsing .gff file	
r=open(gff,"r")
sites={}
for line in r:	
	if line[0]==">":
		break
	if line[0]!="#":
		cols=[x.strip().strip("\"") for x in line.split("\t")]
		cigar=cols[8].split(";")
		vals=[x.split("=")[1] for x in cigar if x.find("=")!=-1]
		keys=[x.split("=")[0] for x in cigar if x.find("=")!=-1]			
		if cols[0]==chrom and (cols[2]=="gene" or cols[2]=="pseudogene" ):
			sites[vals[keys.index("locus_tag")]]=cols[3:7]
# Setting centre of ori
try:
	oril=[int(x.replace(",","")) for x in ori.split(" ... ")]
except:
	exit("ori should be in this format: \"4,647,931 ... 204\"")
if len(oril)!=2 :
	exit("ori should be in this format: \"4,647,931 ... 204\"")
if oril[1]<oril[0]:
	ori=oril[0]+int(((gl-oril[0])+oril[1])/2)
	if ori>=gl:
		ori-=gl
else:
	ori= int(sum(oril)/2)	
#Setting gene coordinates			
coor={}
if dnaA not in sites :
	exit("Given locus tag is not included in the .gff file")
dnaAsign=sites[dnaA][3]
for site in sites:
	coor[site]=sgrp(gl,ori,dnaAsign,sites[site])
dst={k:float(coor[k][:-1]) for k in coor}
dnaApt=[k for k, v in sorted({k:dst[k] for k in dst if dst[k]>0}.items(),key=lambda x:x[1])][0]
if dnaApt!=dnaA:
	exit("Data are not suitable for setting gene coordinates.")
w=open(bcc,"w")
w.write("# Bacterial chromosome coordinates of "+chrom+"\n")
for site in sites:
	w.write("\t".join([site,coor[site]])+"\n")

