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
"""
sequenceReads.py stores and processes individual DNA sequence reads.
"""
try: 
    from dbcAmplicons import editdist
    editdist_loaded = True
except ImportError:
    print("Warning: editdist library not loaded, Insertion/Deletion detetion in barcodes and primers will not be performed")
    editdist_loaded = False 


#------------------- calculate distance between two barcode sequences ------------------------------
def barcodeDist(b_1, b_2):
    """
    gets edit distnace between two equal-length barcode strings
    """
    if editdist_loaded :
        return editdist.distance(b_1,b_2)
    elif len(b_1) == len(b_2) and len(b_1) > 0:
        return sum(map(lambda x: x[0] != x[1], zip(b_1, b_2) ))
    else:
        print "ERROR:[barcodeDist] lengths of barcodes and index read do not match!"
        print "Target", b_1
        print "Index read:", b_2
        raise Exception

#------------------- calculate distance between primer sequence and first part of read ------------------------------
def primerDist(primer,read, max_diff, end_match):
    """
    gets edit distance between primer and read sequence
    """
    if editdist_loaded:
        return editdist.bounded_distance(primer,read,max_diff, end_match)
    else:
        read = read[0:len(primer)]
        if end_match == 0:
            dist = sum(map(lambda x: x[0] != x[1], zip(read, primer) ))
        else:
            dist = sum(map(lambda x: x[0] != x[1], zip(read[:-end_match], primer[:-end_match]) ))
            for i in range(len(primer)-end_match,len(primer)):
                if read[i] != primer[i]:
                    dist += 100
                    return [dist,len(primer)]
        return [dist,len(primer)]


# ---------------- Class for 4 read sequence data with double barcode set ----------------
class FourSequenceReadSet:
    """ 
    Class to hold one Illumina four read (double barocodes) set, only one read names is used for all sequence reads
    Class processes a read by defining the barcode pair and primer pair [optional] and project [optional]. Finally
    class returns a paired read set for output
    """
    def __init__(self,name,read_1,qual_1,read_2,qual_2,bc_1,bc_2):
        """
        Initialize a FourSequenceReadSet with a name, two read sequences and cooresponding quality sequence and 
        two barcode sequences. A read is initially defined as 'not' a good read and requires processing before 
        being labeled as a good read.
        """
        self.barcode = [None,0,0] ## when filled, a vector of length 3 [barcodePairID, editdist1, editdist2]
        self.primer = [None,None,0,0,None,0,0] ## when filled, a vector of length 7 [primerPairID,primerID1,editdist1,readendpos1, primerID2,editdist2,readendpos2]
        self.sample = None
        self.project = None
        self.name = name
        self.read_1 = read_1
        self.qual_1 = qual_1
        self.read_2 = read_2
        self.qual_2 = qual_2
        self.bc_1 = bc_1
        self.bc_2 = bc_2
        self.goodRead = False
    def assignBarcode(self, bcTable, max_diff):
        """
        Given a barcodeTable object and the maximum number of allowed difference (mismatch, insertion, deletions)
        assign a barcode pair ID from the reads barcodes.
        """
        ### Barcode One Matching ###
        bc1 = None
        bc1Mismatch = len(self.bc_1)
        for key in bcTable.getP7():
            bcdist = barcodeDist(key, self.bc_1)
            if bcdist < bc1Mismatch:
                bc1 = key
                bc1Mismatch = bcdist
        ### Barcode Two Matching ###
        bc2 = None
        bc2Mismatch = len(self.bc_2)
        for key in bcTable.getP5():
            bcdist = barcodeDist(key, self.bc_2)
            if bcdist <= bc2Mismatch:
                bc2 = key
                bc2Mismatch = bcdist
        ### Barcode Pair Matching ###
        self.barcode = [bcTable.getMatch(bc1,bc2),bc1Mismatch,bc2Mismatch]
        if self.barcode[0] != None and bc1Mismatch <= max_diff and bc2Mismatch <= max_diff:
            self.goodRead = True
        return 0
    def assignPrimer(self, prTable, max_diff, endmatch):
        """
        Given a primerTable object, the maximum number of allowed difference (mismatch, insertion, deletions) and the
        required number of end match bases (final endmatch bases must match) assign a primer pair ID from the read
        sequences.
        """
        ### Barcode One Matching ###
        pr1 = None
        pr1Mismatch = max_diff+1
        pr1Position = 0
        for key in prTable.getP5sequences():
            prdist = primerDist(key, self.read_1,max_diff,endmatch)
            if prdist[0] < pr1Mismatch:
                pr1 = key
                pr1Mismatch = prdist[0]
                pr1Position = prdist[1]
        ### Barcode Two Matching ###
        pr2 = None
        pr2Mismatch = max_diff+1
        pr2Position = 0
        for key in prTable.getP7sequences():
            prdist = primerDist(key, self.read_2,max_diff,endmatch)
            if prdist[0] <= pr2Mismatch:
                pr2 = key
                pr2Mismatch = prdist[0]
                pr2Position = prdist[1]
        ### Barcode Pair Matching ###
        combined_pr = prTable.getMatch(pr1,pr2)
        self.goodRead = self.goodRead and combined_pr[0] != None and pr1Mismatch <= max_diff and pr2Mismatch <= max_diff
        self.primer = [combined_pr[0],combined_pr[1],pr2Mismatch,pr1Position,combined_pr[2],pr2Mismatch,pr2Position]
        return 0
    def assignRead(self, sTable):
        """
        Given a samplesTable object, assign a sample ID and project ID using the reads barcode and primer designation
        """
        self.project = sTable.getProjectID(self.barcode[0],self.primer[0])
        self.sample = sTable.getSampleID(self.barcode[0],self.primer[0])
        if self.project == None:
            self.goodRead = False
        return 0
    def getBarcode(self):
        """
        Return the reads barcode pair ID
        """
        return self.barcode[0]
    def getPrimer(self):
        """
        Return the reads primer pair ID
        """
        return self.primer[0]
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
    def getRead(self):
        """ 
        Create four line string ('\n' separator included) for the read pair, returning a length 2 vector (one for each read)
        """
        if self.primer[0] != None:
            read1_name = "%s 1:N:0:%s:%s %s|%s|%s|%s %s|%s|%s" % (self.name, self.barcode[0], self.primer[0], self.bc_1, self.barcode[1], self.bc_2 , self.barcode[2], self.primer[1], self.primer[2], self.primer[3])
            read2_name = "%s 2:N:0:%s:%s %s|%s|%s|%s %s|%s|%s" % (self.name, self.barcode[0], self.primer[0], self.bc_1, self.barcode[1], self.bc_2 , self.barcode[2], self.primer[4], self.primer[5], self.primer[6])
            r1 = '\n'.join([read1_name, self.read_1[self.primer[3]:],'+',self.qual_1[self.primer[3]:]])
            r2 = '\n'.join([read2_name, self.read_2[self.primer[6]:],'+',self.qual_2[self.primer[6]:]])
        else:
            read1_name = "%s 1:N:0:%s %s|%s|%s|%s" % (self.name, self.barcode[0], self.bc_1, self.barcode[1], self.bc_2 , self.barcode[2])
            read2_name = "%s 2:N:0:%s %s|%s|%s|%s" % (self.name, self.barcode[0], self.bc_1, self.barcode[1], self.bc_2 , self.barcode[2])
            r1 = '\n'.join([read1_name, self.read_1,'+',self.qual_1])
            r2 = '\n'.join([read2_name, self.read_2,'+',self.qual_2])
        return [r1,r2]

# ---------------- Class for 2 read sequence data processed with dbcAmplicons preprocess ----------------
class TwoSequenceReadSet:
    """ 
    Class to hold one Illumina two read set, assumed to have already been preprocessed with Barocode and Primers already
    identified. Class processes a read by defining sample and project ids. Finally class returns a paired read set for output.
    """
    def __init__(self,name_1,read_1,qual_1,name_2,read_2,qual_2):
        """
        Initialize a TwoSequenceReadSet with names, two read sequences and cooresponding quality sequence.
        Barcode and primer sequence is inferred by their placement in the read names.
        A read is initially defined as 'not' a good read and requires processing before being labeled as a good read.
        """
        try:
            split_name = name_1.split(" ")
            self.name = split_name[0]
            self.barcode = split_name[1].split(":")[3]
            self.barcode_string = split_name[2]
            self.sample = None
            if (len(split_name) == 4):
                self.primer_string1 = split_name[3]
                self.primer_string2 = name_2.split(" ")[3]
                self.primer = split_name[1].split(":")[4]
            else:
                self.primer_string1 = None
                self.primer_string2 = None
                self.primer = None
            self.read_1 = read_1
            self.qual_1 = qual_1
            self.read_2 = read_2
            self.qual_2 = qual_2
            self.goodRead = False
        except IndexError:
            print 'ERROR:[TwoSequenceReadSet] Read names are not formatted in the expected manner'
            raise
        except:
            print 'ERROR:[TwoSequenceReadSet] Unknown error occured initiating read'
            raise
    def assignRead(self, sTable):
        """
        Given a samplesTable object, assign a sample ID and project ID using the reads barcode and primer designation
        """
        self.project = sTable.getProjectID(self.barcode,self.primer)
        self.sample = sTable.getSampleID(self.barcode,self.primer)
        if self.project != None:
            self.goodRead = True
        return 0
    def getRead(self):
        """ 
        Create four line string ('\n' separator included) for the read pair, returning a length 2 vector (one for each read)
        """
        if self.primer != None:
            read1_name = "%s 1:N:0:%s:%s %s %s" % (self.name, self.barcode, self.primer, self.barcode_string, self.primer_string1)
            read2_name = "%s 2:N:0:%s:%s %s %s" % (self.name, self.barcode, self.primer, self.barcode_string, self.primer_string2)
        else:
            read1_name = "%s 1:N:0:%s %s" % (self.name, self.barcode, self.barcode_string)
            read2_name = "%s 2:N:0:%s %s" % (self.name, self.barcode, self.barcode_string)
        r1 = '\n'.join([read1_name, self.read_1,'+',self.qual_1])
        r2 = '\n'.join([read2_name, self.read_2,'+',self.qual_2])
        return [r1,r2]

# ---------------- Class for 2 read sequence data processed with dbcAmplicons preprocess ----------------
class OneSequenceReadSet:
    """ 
    Class to hold a one Illumina read set, assumes the paired reads produced by dbcAmplicons preprocess have been merged
    """
    def __init__(self,name_1,read_1,qual_1):
        """
        Initialize a OneSequenceReadSet with name, one read sequences and cooresponding quality sequence.
        A read is initially defined as 'not' a good read and requires processing before being labeled as a good read.
        """
        self.goodRead = False
        self.primer = None
        ## parse the read name for the sample and primer ids
        try:
            split_name = name_1.split(" ")
            self.name = split_name[0]
            self.barcode = split_name[1].split(":")[3]
            if (len(split_name) == 4):
                self.primer = split_name[1].split(":")[4]
        except IndexError:
            print 'ERROR:[OneSequenceReadSet] Read names are not formatted in the expected manner'
            raise
        except:
            print 'ERROR:[OneSequenceReadSet] Unknown error occured initiating read'
            raise            
        self.read_1 = read_1
        self.qual_1 = qual_1
    def getRead(self):
        """ 
        Create four line string ('\n' separator included) for the read pair, returning a length 1 vector (one read)
        """
        if self.primer != None:
            read1_name = "%s 1:N:0:%s:%s" % (self.name, self.barcode, self.primer)
        else:
            read1_name = "%s 1:N:0:%s" % (self.name, self.barcode)
        r1 = '\n'.join([read1_name, self.read_1,'+',self.qual_1])
        return [r1]


