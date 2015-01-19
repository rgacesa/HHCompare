#!/usr/bin/python2.7
# encoding: utf-8
'''
HHCompare -- code for HMM based paralogues analysis

orthoHH is a pipeline for HMM-HMM comparison based hierarchial clustering 
and analysis of potential paralogues in sequence set.
It is based on premise that paralogues sequences have very high similarity
(below certain e-value cutoff for HMM vs HMM alignment). 

code generates HMM models from input (set of sequences), compares them
all vs all, groups similar ones and then repeats; simplified workflow is as following:

1) generate HMMs from sequences or groups
2) do HMM-HMM comparison between HMMs generated under 1)
3) evalute pairs for similarity, merge pair(s) with similarity below CUTOFF in group(s)
4) repeat 1) if any merge was done, finish run if no merge was done
5) parse results and generate hierarchial trees of groups 

Notes: 
- input: multiple protein sequences in FASTA format; not tested for nucleic acid sequences
- output: set of trees respresenting potential paralogue clusters and unclustered sequences
- note that code DOES NOT provide graphical representation of results, it generates newick outputs 
for detected groups of clustered sequences; various tools can convert it into tree-like graphics
(example: http://etetoolkit.org/treeview/)

Dependencies: 
- RBioTools.py : various bioinformatics functions; should be in same folder as code
- HH-suite2.0, installed locally; make sure paths are set property when running
    download at ftp://toolkit.genzentrum.lmu.de/pub/HH-suite/, follow appropriate instructions for install
- ClustalW2, installed locally;make sure paths are set property when running
    - download at ftp://ftp.ebi.ac.uk/pub/software/clustalw2/
- python2.7

More notes: 
- last update: 19/01/2015
- if python is not in /usr/bin/python2.7, code will not self-execute; run it as <python> <codename> or edit 1st line
- in order to avoid entering paths each time code is run, change lines 329-331 default= to appropriate default
- tested under Ubuntu 12.04.5 LTS
- to annotate results, use HHCompareAnnotator.py (note: script is in early prototype version, and might need tweaking)
- run ./HHCompare -h for help

Example run: 
- ./HHCompare.py -I ./test.fa -O __test_out --HHMAKE <path> --HHALIGN <path> --CLUSTALW2 <path> 


@author:     Ranko Gacesa

@copyright:  2014 King's College London. All rights reserved.

@license:    Attribution-NonCommercial-ShareAlike 4.0 International (http://creativecommons.org/licenses/by-nc-sa/4.0/)
	Distributed on an "AS IS" basis without warranties or conditions of any kind, either express or implied.	

@contact:    ranko.gacesa@kcl.ac.uk
@deffield    updated: Updated
'''

import sys
import os
import csv
from math import log10

# parallel processor for shell
import sh

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

from RBioTools import loadAllFastas
from RBioTools import hashFastas
from RBioTools import getFaHeader
from RBioTools import getFaSeq
from RBioTools import saveFastas

import copy

__all__ = []
__version__ = 0.2
__date__ = '2014-11-18'
__updated__ = '2015-01-19'

'''-------  timeLineToNewick ---------------------------------------- 
Parser for output
-> simply grab output file
and return which fastas are paired (best pairs only)
pairing is defined as having match < cutoff, in case
of multiple matches < cutoff, take best one
'''

'''  
-> basically do following: 
-> 1) check last line, find all groups with at least 2 members
   -> we have that many trees in our set => trees
-> 1b) check last line, find all groups with only 1 member
   -> there are "orphan trees", basically genes without orthologues
      -> record them, but don't bother afterwards
      
-> 2) iterate over trees (groups detected in 1) )
   -> go start -> end, starting with point where we have at least
   2 members in the group: this is (A,B); for each change, add another
   member in format (<previous>,X)
   
   returns: [[orfan_names],[trees]]
            ex: [['A','B','H'],['(C,D)','((G,I),J)']  
---------------------------------------------------------------------  '''
def timeLineToNewick(tl,mfasta,hdrtype,hdrFix=0):		
	inOrfID = set() # ids in orfans
	inTreeID = set() # ids in trees
	
	orfN = [] # orfan names, when headers are extracted from mfasta
	treeGID = [] # all tree group coordinates (group #, timeline row]
	''' 
	first find orfans and proper trees
	'''
	tln = len(tl)	
	for tlline in reversed(tl):
		gc = -1
		tln-=1		
		#print tln,tlline		
		for g in tlline:
			gc+=1
			# orfan
			if len(g) == 4:
				#try:
					if not g[3] in inOrfID and not g[3] in inTreeID:
						inOrfID.add(g[3])
						if hdrFix==1:
							fH = getFaHeader(mfasta[g[3]]).replace(' ','_')
							orfN.append(fH)						
						else:
							orfN.append(getFaHeader(mfasta[g[3]]))
				#except:
				#	print "ERRROER", g
				#	exit()				
			# tree
			if len(g) > 4:
				for i in g[3:]:
					inTreeID.add(i)
				treeGID.append([gc,tln])
		#treeGID.a(treeGIDl)
	#print 'all orfans: ', inOrfID
	#print 'in trees :', inTreeID
	#print 'total: ',len(inOrfID)+len(inTreeID)
	#print 'groups with potential trees:',treeGID 
	#exit()
	''' 
	now pick up all actual trees 
	'''
	trec = []
	for g in treeGID:
		pt = tl[g[1]][g[0]]
		trec.append(pt)
	'''
	trec = []
	for g in range(0,len(tl)):
		trecl = []
#		print 'tree ',g
		for ti in range(1,len(tl)):
			if not tl[ti][g] == tl[ti-1][g] and not tl[ti][g] == [False]:
				trecl.append(tl[ti][g])  
		if not trecl == []:
			trec.append(trecl)	
	tn = 0
	'''	
	''' trec: holds all "potential" trees
	-> now just remove trees entirely contained within other trees!
	'''
	trec_fin = []
	#print ' -- all potential trees -- '
	#for i in trec: 
	#	print i
	#print ' ------- resolving ------ '
	# remove duplicates
	trec_hesh = set()
	for n in range(0,len(trec)):
		trec_hesh.add(tuple(trec[n][3:]))
	#for t in trec_hesh:
	#	print t
	trec_nd = []
	for t in trec_hesh: 
		for n in range(0,len(trec)):
			if tuple(trec[n][3:]) == t: 
				trec_nd.append(trec[n])
				break
	#print ' --- not sup --- '
	#print 'd',len(trec),'nd',len(trec_nd),'hesh',len(trec_hesh)
	#for i in trec_nd: 
	#	print i
		
	# remove trees contained in other trees
	trec = copy.deepcopy(trec_nd)
	for n1 in range(0,len(trec)):
		doAdd = True
		for n2 in range(0,len(trec)):
			if not n1 == n2:				
				if set(tuple(trec[n1][3:])) < set(tuple(trec[n2][3:])):
					doAdd = False
					
					#print n1,trec[n1][3:],' < ',n2,trec[n2][3:]
					break
		if doAdd: 
			trec_fin.append(trec[n1]) 								
	finalT = []
	#print ' --- resolved trees --- '
	#for t in trec_fin: 
	#	print t
	#print ' --- final check: ---- '
	tit = 0
	for i in trec_fin: 
		tit+=len(i)-3
	#print '  -> total in orfans: ',len(orfN),'in trees',tit,'total',tit+len(orfN),' should be ',len(mfasta)
	#print '  -> total DIs in orfans: ',len(inOrfID),'in trees',len(inTreeID),'total',len(inOrfID)+len(inTreeID),' should be ',len(mfasta)
	#print ' --- looking for overlaps --- '
	#for n1 in range(0,len(trec_fin)):
	#	for n2 in range(0,len(trec_fin)):
	#		if not n1 == n2:
	#			for i1 in trec_fin[n1][3:]:
	#				if i1 in trec_fin[n2][3:]:
	#					print 'overlap: ',trec_fin[n1],trec_fin[n2]
	#for n1 in range(0,len(trec_fin)):
	#	for i1 in trec_fin[n1][3:]:
	#		if i1 in inOrfID: 
	#			print 'overlap',trec_fin[n1],inOrfID
	
	# replace ids with actual names
	for tc in trec_fin:		
		#d print tc		
		ntw = tc[2]
		ids = tc[3:]
		for i in ids:
			A = [('(',','),(',',')'),(',',':'),('):',','),('(',':'),(':',',')]
			for (a,b) in A: 
				#d print 'checking',a,'->',b
				if not ntw.find(a+str(i)+b) == -1:
					chks = ntw.find(a+str(i)+b)+1
					chke = chks+len(str(i))
					if hdrtype == 2:
						if hdrFix==1:
							replacement =  getFaHeader(mfasta[i])[1:].replace('(','_').replace(')','_').replace(':','_').replace(' ','_')
						else:												
							replacement = getFaHeader(mfasta[i])[1:].replace('(','_').replace(')','_').replace(':','_')						
					elif hdrtype == 1:
						replacement = getFaHeader(mfasta[i])[1:getFaHeader(mfasta[i]).find(' ')].replace('(','_').replace(')','_').replace(':','_')
					ntw = ntw[:chks]+replacement+ntw[chke:]
		ntw = ntw+';'
		finalT.append(ntw)
		#print 'tree',tn,'->',trecl	
	return [orfN,finalT]
''' --------------- OTHER SMALL FUNCTIONS ------------ '''
def decodePair(pair):
#	('_out_iter1_seq_0008.hhm', '_out_iter1_seq_0003.hhm', 4.3e-130)
	p0 = pair[0]
	p0 = int(p0[p0.find('_seq_')+5:-4])
	p1 = pair[1]
	p1 = int(p1[p1.find('_seq_')+5:-4])
	return [p0,p1,pair[2]]

''' --------------------------------------------
HH comparison parser
------------------------------------------------ '''
def parseHHcomparison(inFile,cutoff=1.0e-10,bestNr=1000):
	cutoff = float(cutoff)
	cutoff_safe = cutoff
	above_cutoff_safe = []     # (Q,T,eV)
	checked = set()
	getRes = False
	#ret = []	
	with open(inFile) as iF:
		c = 0
		for l in iF:
			l = l.strip()
			c+=1
			if c % 2 == 1 and ' VS ' in l:
#				print c,l
				inQ = l.split(' VS ')[0].strip()
				inT = l.split(' VS ')[1].strip()
				if not inQ == inT:
					getRes = True
			elif c % 2 == 0 and getRes:
				s = sorted([inQ,inT])[0] + sorted([inQ,inT])[1]
				if not s in checked:
					checked.add(s)
					pV = l[l.find('P-value = ')+10:]					
					pV = float(pV)
					#d print pV,cutoff_safe  
					cn = 0
					if pV <= cutoff_safe:
						#d print pV,'<',cutoff_safe
						for i in above_cutoff_safe:
							if i[0] == inQ:
								cn+=1													
						if cn < bestNr:
							above_cutoff_safe.append((inQ,inT,pV))						
					else:
						pass 				
				getRes = False

	#print ' cutoff ',cutoff_safe
	return above_cutoff_safe

''' --------------------------------------------
# HH vs HH comparison function via HHAlign, 
implemented using Thread Pool (from multiprocessing)
---------------------------------------------- '''
from multiprocessing import Pool
def doHHComparison(hhFile1,hhFile2,outFile,HHALIGN,tmpFolder,resFile):
	hhAlign = sh.Command(HHALIGN) 
	grep = sh.Command("grep")
	#print a,'vs',b,';',
	af = tmpFolder+'/'+hhFile1
	bf = tmpFolder+'/'+hhFile2		
	#print '      ->', hhFile1,'VS',hhFile2								
	job = ([grep(hhAlign("-i",af,"-v","2","-t",bf,_bg=False),"P-value",_bg=False),hhFile1,hhFile2])
	job[0].wait()
	return [str(job[0]),str(job[1]),str(job[2])]
#		str(job[1]+' VS '+job[2]),str(job[0])]
''' ----------------------------------------
# MAIN
---------------------------------------------- '''
def main(argv=None):
	'''Command line options.'''

	if argv is None:
		argv = sys.argv
	else:
		sys.argv.extend(argv)
		
	program_version = "v%s" % __version__
	program_build_date = str(__updated__)
	program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
	program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
	program_license = '''%s

	Created by Ranko Gacesa on %s.
	Copyright 2014 King's College London. All rights reserved.

 	Licence: Attribution-NonCommercial-ShareAlike 4.0 International 
 	(http://creativecommons.org/licenses/by-nc-sa/4.0/)
	Distributed on an "AS IS" basis without warranties 
	or conditions of any kind, either express or implied.
	see README orHHCompare.py for details	

USAGE
''' % (program_shortdesc, str(__date__))

	# Setup argument parser
	parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
	parser.add_argument('-V', '--version', action='version', version=program_version_message)
	parser.add_argument('-I','--input', required=True, help='input file (set of seqs in fasta file)')
	parser.add_argument('-O','--output', default='_out', help='output file [def: out]') # name is _out_iter_x_seq_y.fa 
	parser.add_argument('--cutoff',default=1.0e-10,help='HMM-HMM comparison e-value [def: 1.0e-10]')
	parser.add_argument('--HHMAKE', default='/home/ranko/Development/tools/HHsuite/bin/hhmake', help='full path to hhmake of HHsuite2.0!')
	parser.add_argument('--HHALIGN', default='/home/ranko/Development/tools/HHsuite/bin/hhalign', help='full path to hhalign of HHsuite2.0!')
	parser.add_argument('--CLUSTALW2', default='/home/ranko/Development/tools/clustal/clustalw2', help='full path to clustalw2!')
	parser.add_argument('--tmpFolder', default='./hh_tmp',help='temporary files folder; def: ./hh_tmp')
	parser.add_argument('--nrThreads', default=20,help='max, number of threads to open; recommended = 2 x CPUs; def: 20')
	parser.add_argument('--timeLineOut', default='__TIMELINE.txt',help='output file for timeline; def: __TIMELINE.txt')
	parser.add_argument('--resultsOut', default='__results.txt',help='output file for results; def: __results.txt; note: placed in current folder, not in tmp')
	parser.add_argument('--headersOut', default=2,type=int,help='how to label trees: 1: ID only, 2: whole header; def: 2')
	parser.add_argument('--doDistances', default=1,type=int,help='if 1, writes distances to trees; distance = 1/(1/log(eV/100))*10; def = 1')
	parser.add_argument('--doEValues', default=1,type=int,help='if 1, writes e-values to trees; def = 1')
	parser.add_argument('--hdrFix',default=1,type=int,help='if 1, replaces spaces in header with underscores for Newick compatibility; def = 1')
	parser.add_argument('--maxMerges',default=1,type=int,help='if  > 1, uses fast mode, merging more then one pair per run; reduces # of iterations, but might reduce accuracy a bit; def = 1')
	# Process arguments
	args = parser.parse_args()
	# cleanup tmp
	os.system('rm -r -f ./'+args.tmpFolder)		
	os.system('mkdir ./'+args.tmpFolder)
	# load input		
	''' generate master fasta list '''
	mFastas = loadAllFastas(args.input)
	#d for f in mFastas:
	#d 		print getFaHeader(f),getFaSeq(f)
	''' generate timeline 
	basically timeline tells us which fastas are grouped in which way
	during which iteration
	format: [[T/F,eV,NEW,n1,n2,...,nx],[],[],...,[]] where each sublist is list of fasta numeric IDs followed by true
	or false tag which describes if this group is done with or not; if False, we reached root; 
	after T/F flag is evalue of current pair addition (if any); -1 if none
	NEW = newick file of current shape of given tree - used to mix, splice and merge trees
	'''
	timeline = []
	timeline_line = []
	for f in range(0,len(mFastas)): 
		timeline_line.append([True,-1,str(f),f])
	timeline.append(timeline_line)
		
	#d print timeline
		
	''' start iterations '''
	doGo = True	
	iteration = 0
	while doGo:
		doGo = False
		''' 1): dump to files 
		- go through timeline of appropriate iteration, 
		dump to file every fasta with assigned indexes from groups
		that have True on [0]
		'''
		gc = -1
		print ' --------------------------'		
		print '    ---> ITERATION ',iteration,'<---'
		print ' --------------------------'
		print '    --> ',iteration,': saving FASTAS '
		saved = 0		
		for g in timeline[iteration]:			
			print g,
			gc+=1
			if g[0]: 				
				fatos = [] # fastas to save list				
				for f in g[3:]:
					fatos.append(mFastas[f])
				if len(g) > 4:
					filename = args.tmpFolder+'/'+'_i_'+str(iteration)+'_g_'+str(gc)+'.fa'				
				else:
					filename = args.tmpFolder+'/'+'_i_'+str(iteration)+'_g_'+str(gc)+'.mfa'		
					#print fatos
				saveFastas(fatos, filename, 60)
				saved+=len(fatos)		
		print ''
		print '      --> Done! Saved total of ',saved,'fastas'
		
		''' 1b) generate Multiple Alignments (if needed)
		- check all outputs, generate MAs if needed
		'''
		print '    --> 1): making multiple alignments as needed '
		gc = -1
		for g in timeline[iteration]:
			gc+=1
			if len(g) > 4 and g[0] == True:
				ifa = args.tmpFolder+'/'+'_i_'+str(iteration)+'_g_'+str(gc)+'.fa'
				ofa = args.tmpFolder+'/'+'_i_'+str(iteration)+'_g_'+str(gc)+'.mfa'	
				print '      -> MA on group',gc
				os.system(args.CLUSTALW2+' -INFILE='+ifa+' -OUTFILE='+ofa +' -OUTPUT=FASTA'+' >> clustal_log')
		
		''' 2) generate HMMs 
		- just run HMM generator on all fastas of appropriate iteration
		(_iter_X_) in file name
		'''
		print '    --> 2): generating HMMs '
		v = 0 # debug: just verbosity stuff for HMM generator
		hhMake = sh.Command(args.HHMAKE)
		
		runningjobs = []		
		for fn in os.listdir(args.tmpFolder):			
			#os.system(args.HHMAKE+' -v '+str(v)+' -i '+ffa+' -M 50 -id 95 -o '+fhhm)							
			if '_i_'+str(iteration)+'_' in fn and '.mfa' in fn:								
				#print fn,
				ffa = args.tmpFolder+'/'+fn
				fhhm = ffa.replace('.mfa','.hhm')				
				runningjobs.append(hhMake("-v",str(v),"-i",ffa,"-M",50,"-id",95,"-o",fhhm,_bg=True))								
				if len(runningjobs) >= int(args.nrThreads):
					for job in runningjobs:
						job.wait()
						runningjobs = []
						
		for job in runningjobs:
			for job in runningjobs:
				job.wait()
				runningjobs = []				
		#print ''
		print '      --> Done! '		
		# verbosity of hhmake and hhalign
		v = '0'				
		''' 3) HMM - HMM compare 
		- run all vs all (without repeats) HMM-HMM comparison
		'''
		print '    --> 3): running HMM-HMM comparison <parallel>'		
		listA = set()
		for fn in os.listdir(args.tmpFolder):			
			if '_i_'+str(iteration)+'_' in fn and '.hhm' in fn:
				listA.add(fn)										
		listB = copy.copy(listA)
		outnameres = args.tmpFolder+'/'+'_i_'+str(iteration)+'_hhcomparison.txt'
		
		pool = Pool(processes=int(args.nrThreads))			  # start 4 worker processes						
		with open(outnameres,'w') as oFT:
			pass
		res = []
		for a in listA:
			#print '.',
			listB.remove(a)
			for b in listB:				
				#print a,'vs',b,';',
				af = args.tmpFolder+'/'+a
				bf = args.tmpFolder+'/'+b				
				#print '      ->', a,'VS',b
		#		doHHComparison(a, b, outnameres, args.HHALIGN, args.tmpFolder,outnameres)								
				res.append(pool.apply_async(doHHComparison,[a, b, outnameres, args.HHALIGN, args.tmpFolder,outnameres]))
		for r in res: 
			job = r.get(timeout=600)
#			print '    -> JOB: '+job[1]+' VS '+job[2]+' DONE!'
			print '.',
			with open(outnameres,'a') as oFT:		
				oFT.write(job[1]+' VS '+job[2]+'\n')
				oFT.write(str(job[0]))			
			
		print ''																		
		print '      --> Done! '	
		''' 4a) parse results '''
		#print args.cutoff
		print '    '
		spairs = []
		if os.path.isfile(outnameres):
			spairs = sorted(parseHHcomparison(outnameres, args.cutoff, 1000),key=lambda x: x[2])
		#for p in spairs: 
		#	print p
		#exit()
		''' 5) merge best pair 
		 -> if any pairs found
		 -> if not, we are done with the run!
		'''
		maxMerges = args.maxMerges
		merged = []
		nrMerges = 0
		doMerges = True
		tlchecker = copy.deepcopy(timeline[iteration])
		tlnew = copy.deepcopy(timeline[iteration])
		iteration+=1
		while nrMerges < maxMerges and doMerges:
			doneMerge = False
			if len(spairs) == 0: 
				doMerges = False
				break
			if len(spairs) > 0:
				# copy timeline into new line								
				doGo = True
				pc = -1			
				for p in spairs:
					doMrg = True
					#print p
					pc+=1
					# best pair - merge
					#if not doneMerge:
					if not doneMerge:	
						grp1 = int(p[0].replace('.hhm','').split('_g_')[1])
						grp2 = int(p[1].replace('.hhm','').split('_g_')[1])
						#print merged,grp1,grp2
						#print grp1,'check',tlchecker[grp1][3:]
						for g in tlchecker[grp1][3:]:
							#print 'checking',g,'vs',merged
							if g in merged:
								#print g, 'skip'
								doMrg = False
						#print grp2,'check',tlchecker[grp2][3:]
						for g in tlchecker[grp2][3:]:
							#print 'checking',g,'vs',merged							
							if g in merged:
								#print g, 'skip'
								doMrg = False
						#print p,doMrg		
						if doMrg:
							# not in merged already												
							doneMerge = True	
							for g in tlchecker[grp1][3:]:
								merged.append(g)
							for g in tlchecker[grp2][3:]:
								merged.append(g)
																		
							if grp2 < grp1:
								t = grp2
								grp2 = grp1
								grp1 = t
							
							print '      -> merging',p							
							#print '      -> appending',grp2,'to',grp1
							tlnew[grp1][0] = True
					
							ev = p[2]
							if p[2] <= 1.0e-300: 
								ev = 1.0e-300
							distance = str(round(1/log10(1/(ev/100))*10,4))
					
							# add one - make first pair(s)
							if len(tlnew[grp2]) == 4 and len(tlnew[grp1]) == 4:
								tlnew[grp1].append(tlnew[grp2][3])
								if args.doDistances == 0:
									tlnew[grp1][2] = '('+str(p[2])+','+str(tlnew[grp2][2])+')'
								elif args.doDistances == 1:
									tlnew[grp1][2] = '('+str(tlnew[grp1][2])+':'+distance+','+str(tlnew[grp2][2])+':'+distance+')'
							
								if args.doEValues == 1:
									tlnew[grp1][2] = tlnew[grp1][2]+''+str(p[2])												
					# add to already estabolished group
					
							elif len(tlnew[grp2]) == 4:
								tlnew[grp1].append(tlnew[grp2][3])
								if args.doDistances == 0:
									tlnew[grp1][2] = '('+str(p[2])+','+str(tlnew[grp2][2])+')'
								elif args.doDistances == 1:
									tlnew[grp1][2] = '('+str(tlnew[grp1][2])+':'+distance+','+str(tlnew[grp2][2])+':'+distance+')'
								#							tlnew[grp1][2] = '('+str(tlnew[grp1][2])+ ','+str(tlnew[grp2][2])+':'+distance+')'
								if args.doEValues == 1:
									tlnew[grp1][2] = tlnew[grp1][2]+''+str(p[2])																												
					# merge groups
							else:
								for i in tlnew[grp2][3:]:
									tlnew[grp1].append(i)
								if args.doDistances == 0:							
									tlnew[grp1][2] = '('+str(tlnew[grp1][2])+','+str(tlnew[grp2][2])+')'
								elif args.doDistances == 1:
									tlnew[grp1][2] = '('+str(tlnew[grp1][2])+':'+distance+','+str(tlnew[grp2][2])+':'+distance+')'							
								if args.doEValues == 1:
									tlnew[grp1][2] = tlnew[grp1][2]+''+str(p[2])								
							tlnew[grp1][1] = spairs[pc][2]
							tlnew[grp2] = [False]																								
				nrMerges+=1							
			# end of merging
		#print 'debug: tlnew',tlnew			
		for g in range(0,len(tlnew)):
			tlnew[g][0] = False		
		#print 'debug: tlnew',tlnew						
		for p in spairs:
			gr1 = int(p[0].replace('.hhm','').split('_g_')[1])
			gr2 = int(p[1].replace('.hhm','').split('_g_')[1])
			if len(tlnew[gr1]) > 1:
				tlnew[gr1][0] = True				
			if len(tlnew[gr2]) > 1:
				tlnew[gr2][0] = True
		
		# ------- temp solution --------		
		#for gc in range(0,len(tlnew)): 
		#	if len(tlnew[gc]) >= 4: 
		#		tlnew[gc][0] = True
		# ------------------------------
												
		print 'debug: tlnew',tlnew
		nt = 0
		for g in tlnew:
			if len(g) > 3: 
				nt+=len(g)-3
		print ' total seqs in timeline: ',nt 
			
		#for p in spairs:
		#	print p
		
		timeline.append(tlnew)				
#		print 'debug: printing entire timeline --'
#		ct = 0
#		for t in timeline:
#			ct+=1
#			print ct,t
				# rest of them - make True		
	# final cleanup
	
	print '     --> DUMPING TIMELINE '
	ct = 0
	with open('./'+args.timeLineOut,'w') as oF:
		writ = csv.writer(oF,delimiter=',',quotechar='"')
		for t in timeline:
			ct+=1
			#print ct,t
			writ.writerow(t)
			
	print ' --> CONSTUCTING TREES: '
	with open(args.resultsOut,'w') as oF:	
		ret = timeLineToNewick(timeline, mFastas,args.headersOut,hdrFix=args.hdrFix)
		print '      --> DONE! '
		print ' ---------------------------------------------- '
		print '   -------------- RESULTS OUTPUT ------------ '
		print ' ---------------------------------------------- '
		print '   --> SEQS WITHOUT PARALOGUES: '
		oF.write('   --> SEQS WITHOUT PARALOGUES: '+'\n')
		for n in ret[0]:
			print n
			if args.headersOut==1:
				oF.write(n[:n.find(' ')]+'\n')
			else:
				oF.write(n+'\n')
		print '   --> TREES: '
		oF.write('\n   --> TREES'+'\n')
		cn = 0
		for n in ret[1]:
			cn+=1			
			print cn,n
			print '---------------------'
			oF.write(n+'\n')
	print ' ------------- ANALYSIS DONE ! ---------------'
	print '  -> timeline saved to '+args.timeLineOut
	print '  -> trees and orfans saved to '+args.resultsOut
# -----------------------------------
# ------ DEFINE WHAT TO CALL --------
# -----------------------------------
if __name__ == "__main__":
	main()