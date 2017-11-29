#!/usr/bin/env python
docstring='''
generate_bnumber_mapping.py UP000000625.txt.gz > UP000000625.tsv
    Read UniProt annotation file UP000000625.txt.gz, generate mapping 
    among UniProt accession, Gene name, b number, and sequence
'''

import sys,os
import subprocess
import re

if len(sys.argv)<=1:
    sys.stderr.write(docstring)
    exit()

cmd="%scat %s|grep -P '(^AC)|(^GN)|(^//)|(^SQ)|(^\s+)'"%(
    'z'*sys.argv[1].endswith(".gz"),sys.argv[1])
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
stdout,stderr=p.communicate()

accession_pat=re.compile(
    "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})")
name_pat=re.compile("Name=(\S+)")
synonyms_pat=re.compile("Synonyms=([\s\S]+?);")
syn_name_pat=re.compile("(\S+)")
locus_pat=re.compile("OrderedLocusNames=([\s\S]+?);")
orf_pat=re.compile("ORFNames=([\s\S]+?);")
bnumber_pat=re.compile("^b\d+$")

txt=''
for block in stdout.split("//\n"):
    lines=block.strip().splitlines();
    if len(lines)<=1:
        continue
    accession_list=accession_pat.findall(lines[0])
    GN=''.join([l[len("GN   "):] for l in lines[1:] if l.startswith("GN   ")])
    name_list=[n.rstrip(';') for n in name_pat.findall(GN)]
    for m in synonyms_pat.findall(GN):
        for name in m.split(', '):
            name_list+=list(set([n for n in name.split() if \
                not n.startswith('{') and not n in name_list]))
    bnumber_list=[]
    for m in locus_pat.findall(GN):
        for name in m.split(", "):
            bnumber_list+=[b for b in name.split() if bnumber_pat.match(b)]
    for m in orf_pat.findall(GN):
        for name in m.split(", "):
            for sub_name in name.split():
                bnumber_list+=[b for b in sub_name.split('/'
                    ) if bnumber_pat.match(b)]

    sequence=''.join(block.split("SQ   ")[1].splitlines()[1:]).replace(' ','')

    txt+='\t'.join([','.join(accession_list),','.join(name_list),
                    ','.join(bnumber_list),sequence])+'\n'

sys.stdout.write(txt)
