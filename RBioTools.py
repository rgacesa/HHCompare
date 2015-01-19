''' NEW: Jul 30 2014 '''

''' grab GO annotation (function) '''
'''
Created on 31 Oct 2013

@author: ranko
'''

import tempfile
import pycurl
import os
import csv

def get_filename_parts_from_url(url):
	fullname = url.split('/')[-1].split('#')[0].split('?')[0]
	t = list(os.path.splitext(fullname))
	if t[1]:
		t[1] = t[1][1:]
	return t

def retrieveUrl(url, filename=None):
	if not filename:
		garbage, suffix = get_filename_parts_from_url(url)
		f = tempfile.NamedTemporaryFile(suffix = '.' + suffix, delete=False)
		filename = f.name
	else:
		f = open(filename, 'wb')
	c = pycurl.Curl()
	c.setopt(pycurl.URL, str(url))
	c.setopt(pycurl.WRITEFUNCTION, f.write)
	try:
		c.perform()
	except:
		filename = None
	finally:
		c.close()
		f.close()
	return filename


def getGoFunction(protID):
	funcs = []
	tmpFile = retrieveUrl("http://www.ebi.ac.uk/QuickGO/GAnnotation?protein="+protID+"&format=tsv")
	with open(tmpFile,'r') as iF:
		csvreader = csv.reader(iF,delimiter="\t",quotechar='"')
		for row in csvreader:
			#print row[11]
			if 'function' in row[11].lower():			
				funcs.append(row[7])
	return funcs

'''
Created on Jul 30, 2013

@author: ranko
'''

''' grabs FASTA from GI and blastDB '''

'''

various misc functions
'''

import commands
import string
import os
import random

def idGenerator(size=8, chars=string.ascii_uppercase + string.digits):	
	return ''.join(random.choice(chars) for _ in range(size))

'''
returns lists of fastas (one fasta = one string)
'''
def loadAllFastas(iF,fixSpaces=False,uniProtIDs=False):		
	fastas = []
	isFirstFasta = True
	currentFa = ''
	with open (iF,'r') as ifile:
		for a in ifile:
			if fixSpaces:
				a = a.replace(' ','_')
			if a[0] == '>':
				if uniProtIDs: 
					a = a[:a[4:].find('|')+5]+'\n'
				#do not write first fasta
				if isFirstFasta:
					isFirstFasta = False
				#write every other (but last) fasta
				else:
					fastas.append(currentFa.strip())
					currentFa = ''
			currentFa += a
	# write last fasta
		fastas.append(currentFa.strip())
	return fastas 

def hashFastas(fastas,fixSpaces=False):
	hFastas = {}
	for f in fastas:
		if fixSpaces:
			hFastas[getFaHeader(f).replace(' ','_').strip()] = getFaSeq(f).replace('\n','').replace('\r','').replace('\c','').replace(' ','_').strip()
		else:	
			hFastas[getFaHeader(f).strip()] = getFaSeq(f).replace('\n','').replace('\r','').replace('\c','').strip()
	return hFastas

def getFaHeader(fa):
	return fa[0:fa.find('\n')].replace('\n','').replace('\r','').replace('\c','').strip()

def getFaSeq(fa):
	return fa[fa.find('\n')+1:].replace('\n','').replace('\r','').replace('\c','').strip()

def splitToChunks(fa,step):
	return [fa[i:i+step] for i in range(0, len(fa), step)]

def saveFasta(fa, filename,maxlength=100000):
	with open(filename,'w') as oF:
		oF.write(getFaHeader(fa)+'\n')
		sp = splitToChunks(getFaSeq(fa), maxlength)
		for l in sp:
			oF.write(l+'\n')

def saveFastas(falist, filename,maxlength=100000):
	with open(filename,'w') as oF:
		for fa in falist:
			oF.write(getFaHeader(fa)+'\n')
			sp = splitToChunks(getFaSeq(fa), maxlength)
			for l in sp:
				oF.write(l+'\n')

def fastaToPir(fastaseq):
	hdr = getFaHeader(fastaseq)
	seq = getFaSeq(fastaseq)
	hdrT = hdr.strip().replace('>','')	
	hdrPir = '>P1;'+hdrT+'\n'
	hdrPir2 = 'sequence:'+hdrT+':::::::0.00:0.00\n'
	seqPir = seq
	qLength = len(seq)		
	if seq.strip()[-1] == '*':
		pass
	else:
		seqPir += '*'
	pir = hdrPir+hdrPir2+seqPir
	return pir

def argToList(argI):
    argI = argI.strip()
    ret = []
    if not (argI[0] == '[' and argI[-1] == ']'):
        #print 'warning: CL arg',argI,'not list!'
        pass
    elif argI == '[]':
        pass
    else: 
        for i in argI[1:-1].split(','):
            ret.append(i.strip())
    return ret

def argToListI(argI):
    argI = argI.strip()
    ret = []
    if not (argI[0] == '[' and argI[-1] == ']'):
        #print 'warning: CL arg',argI,'not list!'
        pass
    elif argI == '[]':
        pass
    else: 
        for i in argI[1:-1].split(','):
            ret.append(int(i.strip()))
    return ret

def determineSearchFileType(inFile,specific=True):
	outType = 'unknown'
	#print 'detecting file type'
	with open(inFile,'r') as inF:
		lc = 0
		isHHblits = False
		isProtF = True
		isDNAF = True
		for l in inF:
			#print len(l)
			if len(l) > 1:
				lc += 1
				l = string.lower(l).strip()
				# different types of BLAST
				# P or DELTA (notation is same)
				if l[0] == '>' and lc == 1:
					outType = 'fasta'
				if not l[0] == '>' and outType == 'fasta':
					#check for protein fasta				
					for a in l.lower():
						if a not in 'acgturykmswbdhvn*-':
							isDNAF = False
						if a not in 'abcdefghiklmnpqrstvwyxzju':							
							isProtF = False
						if not isDNAF and not isProtF:
							print 'not DNA and not PROTEIN:'
							print 'seq: ',l.lower()
							print 'problem: symbol:',a
							exit(1)
				if '<blastoutput_program>blastp' in l:
					if specific:
						outType = 'blastp'
					else: 
						outType = 'blast'
				elif '<blastoutput_program>psiblast' in l:
					if specific:
						outType = 'psiblast'
					else: 
						outType = 'psiblast'
			# blast N
			# note: mega blast and disc mega blast are the same as blastn
				elif '<blastoutput_program>blastn' in l:
					if specific:
						outType = 'blastn'
					else: 
						outType = 'blast'
			# blast X
				elif '<blastoutput_program>blastx' in l:
					if specific:
						outType = 'blastx'
					else: 
						outType = 'blast'
			# T blast N
				elif '<blastoutput_program>tblastn' in l:
					if specific:
						outType = 'tblastn'
					else: 
						outType = 'blast'                             
			# T blast X
				elif '<blastoutput_program>tblastx' in l:
					if specific:
						outType = 'tblastx'
					else:
						outType = 'blast'                      
			# HMMER
				elif '# jackhmmer' in l:
					if specific:
						outType = 'jackhmmer'
					else: 
						outType = 'hmmer'
				elif '# hmmsearch' in l or '# phmmer' in l:
					if specific:
						outType = 'hmmer'
					else: 
						outType = 'hmmer'
			#HHblits check    
				if lc == 1 and 'query' in l:
					isHHblits = True                                           
				elif lc == 2 and 'match_columns' in l and isHHblits:
					pass
				elif lc == 3 and 'no_of_seqs' in l and isHHblits:
					pass                                     
				elif lc == 4 and 'neff' in l and isHHblits:
					pass                                           
				else:
					isHHblits = False                                           
				if lc >=4 and not outType == 'fasta':
					if isHHblits:
						outType = 'hhblits'
					return outType
	if outType == 'fasta':
		if isProtF and isDNAF:
			print 'WARNING: cannot determine fasta type for ',inFile,'assuming protein!'
			outType = 'fastap'
		elif isProtF:
			outType = 'fastap'
		elif isDNAF:
			outType = 'fastan'
		else:
			print 'WARNING: fasta seems flawed for ',inFile
			print 'BREAKING!'
			print 
			exit (0)
	#print 'success'
	return outType           

def extractFirstGI(data):
    giStart = data.find('gi|')+3
    data = data[giStart:]
    giEnd = data.find('|')
    gi = data[0:giEnd]
    return gi

''' helper function for extracting
possible reference ID 
-> st = input header
-> idname = type of id ('gi' or 'ref' or whatever)
-> allIDs: if True, searches for multiple IDs
'''
def extractId(st,idname,allIDs=True):
    chunk = 'N/D'
    if allIDs:
        maxNR = 1000
    else: 
        maxNR = 1
    fndID = 0
    ret = []
    st = st.lower()
    sI = st.find(idname+'|')
    while fndID < maxNR and not sI == -1:
        sIC = st[sI+len(idname)+1:]
        sE = sIC.find('|')
        if sE == -1:
            chunk = sIC
            fndID +=1
        else:
            chunk = sIC[0:sE]
            fndID +=1
        ret.append(chunk)
        st = sIC[sE+1:]
        sI = st.find(idname+'|')
    return ret

''' parses header, grabs out all
reference IDs among following: 
 -> gi, ref, sp, tr, pdb, gb, emb 
 --> returns dictionary of format: 
    {'gi' : 'a', 'ref': ['b','c'],'pdb': [] ...} 
 --> if trim=True: remove .xxx extension from IDs
 (helps resolve some not found references) 
'''
def parseHeader(hdr,trim=False):
    idDict = {}    
    hdr = string.lower(hdr)
    possibleIDs = ('gi','ref','refseq','sp','swissprot','tr','trembl','emb','embl','pdb','gb','uniprot','uni')
    for i in possibleIDs:
        eID = extractId(hdr,i)      
        if trim:                        
            for iOneID in range(0,len(eID)):
                if not eID[iOneID].find('.') == -1:
                    eID[iOneID] = eID[iOneID][0:eID[iOneID].find('.')]
        if not len(eID) == 0:
            idDict[i] = eID
            
    # SPECIAL RULES FOR PDB PARSING (PDB SUCKS DONKEY ASS...)
    if '|pdbid' in hdr:
        pdbID = hdr[1:hdr.find('|pdbid')]
        if trim:
            if not pdbID.find('.') == -1:
                pdbID = pdbID[0:pdbID.find('.')]        
        idDict['pdb'] = [pdbID]  
    return idDict

def extractAllGIs(data):
    giList = []
    while data.find('gi|') > -1:
        giStart = data.find('gi|')+3
        data = data[giStart:]
        giEnd = data.find('|')        
        try:
            giList.append(int(str(data[0:giEnd])))
        except:
            pass                 
        #giList.append(int(str(data[0:giEnd])))
        data = data[giEnd:]
    return giList    

def extractBlastTitleAllTitles(data):    
    return data.split('>')

def extractBlastTitleGiTitleMap(data):
    giList = []
    titles = data.split('>')
    ret = {}
    while data.find('gi|') > -1:
        giStart = data.find('gi|')+3
        data = data[giStart:]
        giEnd = data.find('|')
        giList.append(int(str(data[0:giEnd])))
        data = data[giEnd:]
    
    for gi in giList:
        if gi not in ret.keys():
            for t in titles:
                t = str(t)
                if str(gi) in t:
                    ret[gi] = t
    return ret

def extractHMMERTitleGiTitleMap(data):
    giList = []
    titles = data.split('| ')
    ret = {}
    while data.find('gi|') > -1:
        giStart = data.find('gi|')+3
        data = data[giStart:]
        giEnd = data.find('|')
        giList.append(int(str(data[0:giEnd])))
        data = data[giEnd:]
    
    for gi in giList:
        if gi not in ret.keys():
            for t in titles:
                t = str(t)
                if str(gi) in t:
                    ret[gi] = t
    return ret

def fetchFastaFromBlastDB (id, blastDBCmd, dbPath):
	cmd = blastDBCmd+' -db '+dbPath+' -entry '+str(id)+' -target_only'
	#print cmd
	#print cmd
	co = commands.getstatusoutput(cmd)[1]	
	#print co
	if 'error' not in co.lower():
		return co
	else:
		print '     -> warning: DB ret [cmd:'+cmd+' {id:'+str(id)+'} got: '+co.lower()
		return 'error'

class WaterResults:
    def __init__(self):
        self.length = ''
        self.gaps = ''
        self.identity = ''
        self.score = ''
        self.identityV = 0.0
        self.similarity = ''
        self.similarityV = 0.0
        self.gapsV = 0.0
        self.longerL = 0
        self.shorterL = 0
    def toStr(self):
        return 'lt:'+self.length+'[L:'+str(self.longerL)+',S:'+str(self.shorterL)+'], gaps:'+self.gaps\
            +'['+str(self.gapsV)+']'+', identity:'+self.identity+'['+str(self.identityV)+']'\
            +', similarity:'+self.similarity+'['+str(self.similarityV)+']'+', score:'+str(self.score)            
    def getCoverage(self):
        return float(self.length)/float(self.shorterL)
    def getNonGapCoverage(self):
        return float(self.similarityV)
            
def doWaterman(waterPath, seqA, seqB, gapOpen='10.0', gapExtend='0.5', outfile='waterTMP.txt', delTmp=False):
    ret = WaterResults()
    L1 = 0
    L2 = 0    
    with open(seqA,'r') as seqAfile:
        lns = seqAfile.readlines()
        c = 0
        seq = ''
        for l in lns:
            c+=1                    
            if c > 1: 
                seq = seq + l.replace(' ','').strip()
        L1 = len(seq)
                
    with open(seqB,'r') as seqBfile:
        lns = seqBfile.readlines()
        c = 0
        seq = ''
        for l in lns:
            c+=1                    
            if c > 1: 
                seq = seq + l.replace(' ','').strip()                
        L2 = len(seq) 
        
    if L1 >= L2:
        ret.longerL = L1
        ret.shorterL = L2
    else:
        ret.longerL = L2
        ret.shorterL = L1
    cmd = waterPath +' -asequence ' +seqA + ' -bsequence '+seqB + ' -gapopen '+gapOpen+' -gapextend '+gapExtend+' -outfile '+outfile
    commands.getstatusoutput(cmd)
    with open(outfile) as toParse:
        rows = toParse.readlines()
        for r in rows:
            r = r.replace(' ','').strip()
            if '#Length:' in r:                
                ret.length = r.split(":")[1]
            elif '#Identity:' in r:
                ret.identity = r.split(":")[1]
                ret.identityV = float(ret.identity.split("(")[1].replace('%','').replace(')',''))/100.0                
            elif '#Similarity:' in r:
                ret.similarity = r.split(":")[1]
                ret.similarityV = float(ret.similarity.split("(")[1].replace('%','').replace(')',''))/100.0                
            elif '#Gaps:' in r:
                ret.gaps = r.split(":")[1]
                ret.gapsV = float(ret.gaps.split("(")[1].replace('%','').replace(')',''))/100.0
            elif '#Score:' in r:
                ret.score = r.split(":")[1]       
    if delTmp:
    	print 'deleting',outfile
    	os.remove(outfile)         
    return ret


class HmmResultRecord:
    def __init__(self, _gi = '', _number = 0, _eV = 0.0, _tax = '', _filename = '', _desc = '', _hmmType = ''):
        self.gi = _gi
        self.number = _number
        self.eVmin = float(_eV)
        self.eVmax = float(_eV)
        self.eVavg = float(_eV)
        self.eVtotal = float(_eV)
        self.eV = float(_eV)
        self.tax = _tax
        self.filename = _filename
        self.desc = _desc
        self.hmmType = _hmmType
        self.fasta = ''
        self.header = '' #fasta header
        self.seq = '' #fasta seq
        self.waterSim = 0.0  # waterman similarity
    def toStr(self): 
        ret = "gi: " + str(self.gi) + " nr: " + str(self.number) + " eV: " + str(self.eV) + " eVavg: "+str(self.eVavg)+" hmmType: "+str(self.hmmType)+ " tax: "+str(self.tax)+" file: "+str(self.filename)+" desc: "+str(self.desc)
        if self.fasta != '':
            ret = ret + "\n"+self.fasta
        return ret
    
    def toSeq(self):
        ret = []
        ret.append(self.gi)
        ret.append(self.number)
        ret.append(self.eVavg)
        ret.append(self.hmmType.strip())
        ret.append(self.tax.strip())
        ret.append(self.desc.strip())
        ret.append(self.fasta.strip())
        return ret       
    def toSeq2(self):
		ret = []
		ret.append(self.gi)
		ret.append(self.number)
		ret.append(self.eVavg)
		ret.append(self.header.strip())
		ret.append(self.seq.strip())
		return ret      	 
    def getCSVheader(self):
        ret = []
        ret.append('gi')
        ret.append('nr.hits')
        ret.append('eVavg')
        ret.append('hmmType')
        ret.append('taxID')
        ret.append('description')
        ret.append('fasta')
        return ret
       
    def getCSVheader2(self):
		ret = []
		ret.append('gi')
		ret.append('nr.hits')
		ret.append('eVavg')
		ret.append('header')
		ret.append('seq')
		return ret
	
    def loadFromCsv(self, csvEntry):
        self.gi = csvEntry[0]
        self.number = csvEntry[1]
        self.eVavg = csvEntry[2]
        self.header = csvEntry[3]
        self.seq = csvEntry[4]

    def loadFromCsv2(self, csvEntry):
		self.gi = csvEntry[0]
		self.number = csvEntry[1]
		self.eVavg = csvEntry[2]
		self.header = csvEntry[3]
		self.seq = csvEntry[4]

def getFastaLt(inFile):
	with open (inFile,'r') as faf:
		ln = 0
		seq = ''
		for l in faf:
			ln+=1
			if ln > 1:
				seq = seq+l.strip()			
		return len(seq)
			

def addEvToFaHdr(header, gi):
	gi = float(gi)
	fb = header.find('| ')+1
	if fb > -1:
		p1 = header[0:fb]
		p2 = '__eV{'
		p3 = str('%.1e' % gi)
		p4 = '}Ve__'
		p5 = header[fb:]
		hdrret = p1+p2+p3+p4+p5
		return hdrret
	else:
		return header

def getEvFromMFaHdr(header):
	eV = -1
	eVPS = header.find('__eV{')+5
	eVPE = header.find('}Ve__')
	eVPStr = header[eVPS:eVPE]
	#print '-H-'
	#print header
	#print '---'
	#print eVPStr
	try:
		f = float(eVPStr)
		return float(eVPStr)
	except: 
		return 10000.0

def mFAHdrToFaHdr(header):
	eVPS = header.find('__eV{')
	eVPE = header.find('}Ve__')+5
	eVPStr = header[0:eVPS]+header[eVPE:] 
	return eVPStr

def isMFAHdr(header):
	ret = False
	#print header.find('__eV{')
	#print header.find('}Ve__')
	if not header.find('__eV{') == -1 and not header.find('}Ve__') == -1:
		ret = True
	return ret

def getPDBID(title):
	at = title
	#print at
	atc = at[at.find('pdb|')+4:]
	#print '-> ',atc
	pdbID = atc[0:atc.find('|')]
	pdbChain = atc[atc.find('|')+1:]
	pdbChain = pdbChain[0:pdbChain.find(' ')]
	return (pdbID,pdbChain)


'''        
print 'fetching FASTA' 
out = fetchFastaFromBlastDB('12345','/home/ranko/bioinfo_tools/ncbi-blast-2.2.27+/bin/blastdbcmd','/home/ranko/bioinfo_tools/DBs/nr_db/nr')
print out
'''
'''
waterPath = 'water'
seqA = 'testA.fa'
seqB = 'testB.fa'
print 'doing WATER'
out = doWaterman(waterPath,seqA,seqB,delTmp=True)
print str(out.toStr())
print "shorter L = "+str(out.shorterL)," longer L = "+str(out.longerL)
out = doWaterman(waterPath,seqB,seqA)
print str(out.toStr())
print "shorter L = "+str(out.shorterL)," longer L = "+str(out.longerL)
print 'coverage',str(out.getCoverage()*100),'%'
print 'nongaps coverage',str(out.getNonGapCoverage()*100),'%'
for f in loadAllFastas('/home/ranko/Development/workspace/DHPipeline/results/POTENTIAL_HITS/testNrf2_allseqs_nr_archaea_merged_hhres_filtered_nd.hhr'):
	print f
	print 'H:',getFaHeader(f)
	print 'S:',getFaSeq(f)
iF = '/home/ranko/Development/workspace/DHPipeline/output/testNrf2_vs_allseqs_nr_archaea_merged-13-11-2013_16:18/POTENTIAL_HITS/testNrf2_vs_allseqs_nr_archaea_merged_p_hits_wf.fa'
print determineSearchFileType(iF)
BLAST_DB_CMD='/home/ranko/Development/tools/blast+/blastdbcmd'
NRDB='/home/ranko/Development/DBs/nr_fasta/nr.fa'
gi = 12837234
res_fasta = fetchFastaFromBlastDB(gi,BLAST_DB_CMD,NRDB)
res_header = getFaHeader(res_fasta)
res_seq = getFaSeq(res_fasta)
print res_fasta
print 'h:',res_header
print 's:',res_seq
mFH = addEvToFaHdr(res_header,1.0e-5) 		
print 'hm:',mFH
print getEvFromMFaHdr(mFH)
print mFAHdrToFaHdr(mFH)

exit(1)
tst = '/home/ranko/Development/workspace/DHPipeline/output/testNrf2_vs_allseqs_nr_archaea_merged-14-11-2013_12:02/POTENTIAL_HITS/testNrf2_allseqs_nr_archaea_merged_hhres_filtered_nd_wf.fa'
fas = loadAllFastas(tst)
for f in fas:
	h = getFaHeader(f)
	print h
	print isMFAHdr(h)
'''