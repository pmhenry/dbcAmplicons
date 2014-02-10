#py-abundance
#!/usr/bin/env python

# Copyright 2013, Institute for Bioninformatics and Evolutionary Studies
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys, os, traceback
import re
import glob
import time
from collections import OrderedDict
from dbcAmplicons import sampleTable
from collections import Counter

class fixrankLine:
    """
    Parse a single line from fixrank formated file, generated by dbcAmplicons classify
    HWI-M01380:50:000000000-A641U:1:2116:11332:23553|Foxtrot267:16S         Bacteria        domain  1.0     "Fusobacteria"  phylum  1.0     "Fusobacteria"  class   1.0     "Fusobacteriales"     order   1.0     "Fusobacteriaceae"      family  1.0     Clostridium XIX genus   0.59
    """
    def __init__(self,line, rank = 'genus', threshold = 0.5):
        """
        Initialize the fixrankLine object given a 'line'
        """
        self.goodRead = True
        self.sample = None
        self.primer = None
        self.taxon = OrderedDict()
        parse = line.split('\t')
        if len(parse) % 3 != 2:
            print("ERROR:[fixrankline] incorrect number of columns in parseing the line, %s" % line)
        name = parse[0].split('|')
        name = name[1].split(':')
        self.sample = name[0] 
        if (len(name) > 1):
            self.primer = name[1]
        levels = len(parse)/3
        self.call = ('Unknown|Root',1.0)
        for i in xrange(levels):
            #self.taxon[parse[i*3+3]] = (re.sub(r'["\']+', "",parse[i*3+2]),float(parse[i*3+4]))
            if float(parse[i*3+4]) >= threshold:
                self.call = (re.sub(r'["\']+', "",parse[i*3+2]) +'|'+ parse[i*3+3],float(parse[i*3+4]))
                if parse[i*3+3] == rank:
                    break
            else:
                break
    def getCall(self):
        """
        Retrieve the call
        """
        return self.call
    def isOk(self):
        """
        Return whether the read is 'good' True or 'bad' False 
        """
        return self.goodRead
    def getSampleID(self):
        """
        Return the reads sample ID
        """
        return self.sample
    def getProject(self):
        """
        Return the reads project ID
        """
        return self.project
    def assignRead(self, sTable):
        """
        Given a samplesTable object, assign a sample ID and project ID using the reads barcode and primer designation
        """
        self.project = sTable.getProjectID(self.sample,self.primer)
        self.sample = sTable.getSampleID(self.sample,self.primer)
        if self.project == None:
            self.goodRead = False
        return 0


class abundanceApp:
    """
    Generate an abundance table from a fixrank formated file
    Takes fixrank formated files from dbcAmplicons classify anda taxonomic rank to build table from, allowable values are (domain, phylum, class, order, family, genus, and species{if performed}
    and output an abundance and proportions table with taxon in rows and samples as columns.
    """ 
    def __init__(self):
        self.verbose=False
    def start(self, fixrank_file, samplesFile, output_prefix='table', rank = 'genus', threshold = 0.5, verbose=True, debug=False):
        """
            Start processing classification fixrank files
        """
        self.verbose = verbose
        evalSample = samplesFile != None
        try:
            lines = 0
            lasttime = time.time()
            self.ffixrank = []
            ## samples
            if evalSample:
                sTable = sampleTable(samplesFile)
                if verbose:
                    print "sample table length: %s, and %s projects." % (sTable.getSampleNumber(),len(sTable.getProjectList()))

            ## check input fixrank files
            for ffile in fixrank_file:
                self.ffixrank.extend(glob.glob(ffile))
                if len(self.ffixrank) == 0 or not all(os.path.exists(f) for f in self.ffixrank):
                    print ('ERROR:[abundance_app] fixrank file(s) not found')
                    raise

            abundanceTable = dict()
            bootscore = dict()
            sampleList = []
            sampleCounts = Counter()
            for ffile in self.ffixrank:
                with open(ffile, "rb") as infile:
                    for line in infile:
                        rank = fixrankLine(line.rstrip('\n'),rank,threshold)
                        if rank.getSampleID() not in sampleList:
                            sampleList.append(rank.getSampleID()) 
                        if evalSample:
                            rank.assignRead(sTable)
                        if rank.isOk():
                            sampleCounts[rank.getSampleID()] +=1
                            tax = rank.getCall()
                            if tax[0] in abundanceTable.keys():
                                abundanceTable[tax[0]][rank.getSampleID()] +=1
                                bootscore[tax[0]] += tax[1]
                            else:
                                abundanceTable[tax[0]] = Counter()
                                abundanceTable[tax[0]][rank.getSampleID()] +=1
                                bootscore[tax[0]] = tax[1]
                        lines += 1
            ## output file
            if evalSample:
                ab_name = output_prefix + '.abundance.txt'
                prop_name = output_prefix + '.proportions.txt'
                sampleList = sTable.getSampleList()
            else:
                ab_name = output_prefix + '.abundance.txt'
                prop_name = output_prefix + '.proportions.txt'
            try:
                abFile = open(ab_name, 'w')
                propFile = open(prop_name, 'w')
            except:
                print ("Can't open files (%s,%s) for writing" % (ab_name,prop_name))
            # write out header line
            sampleList = sorted(sampleList)
            txt = 'Taxon_Name\tLevel\tMeanBootstraptValue\t'+ '\t'.join(sampleList) + '\n'
            abFile.write(txt)
            propFile.write(txt)
            taxa_keys = sorted(abundanceTable.keys())
            for taxa in taxa_keys:
                txt1 = str(taxa.replace('|', '\t')) + '\t' + str(round(bootscore[taxa]/sum(abundanceTable[taxa].values()),3))
                txt2 = str(taxa.replace('|', '\t')) + '\t' + str(round(bootscore[taxa]/sum(abundanceTable[taxa].values()),3))
                for sample in sampleList:
                    txt1 = '\t'.join([txt1,str(abundanceTable[taxa][sample])])
                    txt2 = '\t'.join([txt1,str(round(abundanceTable[taxa][sample]/sampleCounts[sample],3))])
                abFile.write(txt1  + '\n')
                propFile.write(txt2 + '\n')
            cntFile = open(output_prefix + '.taxa_counts.txt', 'w')
            cntFile.write("Taxon_Name\tCount\n")
            for abt in abundanceTable:
                cntFile.write(str(abt) + '\t' + str(sum(abundanceTable[abt].values())) + '\n')
            if self.verbose:
                print "%s lines processed in %s minutes" % (lines,round((time.time()-lasttime)/(60),2))
            self.clean()
            return 0    
        except (KeyboardInterrupt, SystemExit):
            self.clean()
            print("%s unexpectedly terminated" % (__name__))
            return 1
        except:
            self.clean()
            print("A fatal error was encountered.")
            if debug:
                print "".join(traceback.format_exception(*sys.exc_info()))
            return 1

    def clean(self):
        if self.verbose:
            print("Cleaning up.")
        try:
            ## nothing to be done
            pass
        except:
            pass


