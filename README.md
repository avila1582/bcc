# bcc
Program calculates bacterial chromosome coordinates from given .gff file

usage: bcc.py [-h] [-o O] gff ori dnaA

positional arguments:

  gff         Name of input .gff file downloaded from GenBank or created by prokka
              
  ori         ori as set by orifinder (https://tubic.org/Ori-Finder2022/) e.g "4,647,931 ... 204"
              
  dnaA        locus tag of dnaA gene in the .gff file

optional arguments:

  -h, --help  show this help message and exit
  -o O        output file


Example:

bcc.py ATBA_ACICU.gff "6,768 ... 7,437" ACICU_RS00025 -o ATBA_ACICU.bcc
