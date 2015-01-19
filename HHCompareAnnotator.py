'''
Created on 27 Nov 2014

@author:     Ranko Gacesa

@copyright:  2014 King's College London. All rights reserved.

@license:    Attribution-NonCommercial-ShareAlike 4.0 International (http://creativecommons.org/licenses/by-nc-sa/4.0/)

@contact:    ranko.gacesa@kcl.ac.uk
@deffield    updated: Updated

- basically: grabs sequences from inputSeqs
and BLAST annotates them against inputDB
and then adds these annotations to results file
'''
import re
import os
from Bio.Blast import NCBIXML
from RBioTools import loadAllFastas
from RBioTools import hashFastas
import csv

class BRec():
	def __init__(self):
		self.query = ''
		self.target = ''
		self.eV = -1
		self.a_query = ''
		self.a_match = ''
		self.a_subject = ''

def parseBLAST(inFile,parseEVmax):	
	bhits = {}
	#print " --> PARSING BLAST INPUT <--"
	with open(inFile) as resultHandle:       
		blastRecords = NCBIXML.parse(resultHandle)                                    
		for blastRecord in blastRecords:    
			cnt = 0
			for alignment in blastRecord.alignments:
				for hsp in alignment.hsps:
					if hsp.expect <= parseEVmax:
						cnt += 1
						brek = BRec()
						
						brek.query = blastRecord.query
						#brek.query_from = blastRecord.query_start	
						#brek.query_to = blastRecord.query_end
						brek.target = alignment.title					
						brek.eV = hsp.expect
						brek.a_query = hsp.query
						brek.a_match = hsp.match
						brek.a_target = hsp.sbjct															
						try: #blastRecord.query in bhits.keys():
							bhits[brek.query].append(brek)
						except: 
							bhits[brek.query] = [brek]
	return bhits

BLAST = '~/Development/tools/blast+/blastp'

inputRes = 'AD_supplement_5_clustering_e-10_hid.txt'
inputDB = '/home/ranko/Development/DBs/toxprot/toxprot.fa'
inputF =  '/home/ranko/Development/DBs/toxprot/toxprot.fa'
outputek = './__unex_annot_out_e-10.txt'
outputekTree = './__unex_annot_out_tree_e-10.txt'
outputekTreeH = './__unex_nnot_out_tree_ho_e-10.txt'
inputUniProtIDs = True

orfans = []
groups = []  
allhdrs = []
annotated_all = {}	  # {hdr -> annot}
annotated_groups = [] # [[group11[hdr,annot]],...]
annotated_orfans = {}

allText = ''

with open(inputRes) as iF:
	for l in iF: 
		allText = allText+l;
		#print ' --> ',l.strip()
		if l[0] == '>':
			orfans.append(l.strip())
			allhdrs.append(l.strip())
		if l[0] == '(':
			m = re.findall ('\((.+?)\:',l) 
			for m2 in re.findall('\,(.+?)\:',l):
				m.append(m2)	
			gug = set()		
			for i in m:
				j = i.replace('(','').strip()
				if not j[0] == '>': 
					j = '>'+j
					gug.add(j)
					allhdrs.append(j)
			groups.append(gug)

hesh = hashFastas(loadAllFastas(inputF,uniProtIDs=True),fixSpaces=True)
#hesh = hashFastas(inputF)
allTextFull = str(allText)
allTextHrds = str(allText)

#for k in hesh.keys(): 
#	print k,hesh[k]
gc = 0
c = 0
sc = 0
print ' --> ANNOTATING ORFANS:'
# do orfans
for h in orfans:	
	if inputUniProtIDs: 
		h = h[:h[4:].find('|')+5]
	sc+=1
	c+=1
	fO = './tmp/seq_'+str(sc).zfill(4)+'.fa'
	with open(fO,'w') as oF:
		oF.write(h+'\n'+hesh[h])
	#print '   -> blasting'
	os.system(BLAST+' -query '+fO+ ' -db '+inputDB+' -num_threads 4 '+'-outfmt 5 '+' -evalue 1.0e-10 '+' -max_target_seqs 5 -out '+fO.replace('.fa','.xml'))
	#print '      -> done'
	#print '   -> reading '	
	hits = parseBLAST(fO.replace('.fa','.xml'),1.0e-10)
	annots = []
	ano = []
	anoful = []
	for q in hits.keys():
		for h in hits[q]:
			if not h.target in ano:
				anoful.append((h.target,h.eV))
	anoful = sorted(anoful,key=lambda x: x[1])
	print c,q,' -> ',anoful
	if not q in annotated_orfans: 
		annotated_orfans[q] = anoful
		annotated_all[q] = anoful
		print '    -> replacing',q,'with',q+' -> '+anoful[0][0]
		allTextFull = allTextFull.replace(q,q+' -> '+anoful[0][0])
		allTextHrds = allTextHrds.replace(q,q+' -> '+anoful[0][0][:anoful[0][0].find(' ')])		

with open(outputek,'w') as oF:
	vrit = csv.writer(oF,delimiter=',',quotechar='"')
	vrit.writerow(['--------> ANNOTATIONS FOR UNGROUPED SEQUENCES <---------- '])
	c = 0
	for q in annotated_orfans.keys():
		c+=1
		vril = []
		vril.append(c)
		vril.append(q)
		for a in annotated_orfans[q]:
			vril.append(a[0])
		vrit.writerow(vril)
		
# do groups
with open(outputek,'a') as oF:
	vrit = csv.writer(oF,delimiter=',',quotechar='"')
	vrit.writerow(['--------> ANNOTATIONS FOR GROUPS <---------- '])
	print ' annotating groups!'
	gc = 0
	for grp in groups:
		c = 0
		gc+=1
		for h in grp:
			if inputUniProtIDs: 
				h = h[:h[4:].find('|')+5]			
			sc+=1
			c+=1 
			fO = './tmp/seq_'+str(sc).zfill(4)+'.fa'
			with open(fO,'w') as oF:
				oF.write(h+'\n'+hesh[h])
			#print '   -> blasting'
			os.system(BLAST+' -query '+fO+ ' -db '+inputDB+' -num_threads 4 '+'-outfmt 5 '+' -evalue 1.0e-10 '+' -max_target_seqs 5 -out '+fO.replace('.fa','.xml'))
			#print '      -> done'
			#print '   -> reading '	
			hits = parseBLAST(fO.replace('.fa','.xml'),1.0e-10)
			annots = []
			ano = []
			anoful = []
			for q in hits.keys():
				for h in hits[q]:
					if not h.target in ano:
						anoful.append((h.target,h.eV))
			anoful = sorted(anoful,key=lambda x: x[1])
			print gc,c,q,' -> ',anoful
			annotated_all[q] = anoful
			vril = []
			vril.append(gc)
			vril.append(c)
			vril.append(q)
			for a in anoful: 
				vril.append(a[0])
			print '    -> replacing',q,'with',q+' -> '+anoful[0][0]
			allTextFull = allTextFull.replace(q,q+' -> '+anoful[0][0])
			allTextHrds = allTextHrds.replace(q,q+' -> '+anoful[0][0][:anoful[0][0].find(' ')])						
			vrit.writerow(vril)

with open(outputekTree,'w') as oF:
	oF.write(allTextFull)
with open(outputekTreeH,'w') as oF:
	oF.write(allTextHrds)