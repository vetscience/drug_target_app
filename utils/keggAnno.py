#!/usr/bin/env python

import os, sys, optparse, errno
from Utils import KeggTree, Base

#################################################
def options():
    parser = optparse.OptionParser('usage: python %prog -i blastFile -d Kegg -p 1 -r 5 -g geneIdFile')
    parser.add_option('-i', '--bfile', dest='bfile', help='BLAST results file against KEGG or KOBAS database', metavar='BLAST', default='')
    parser.add_option('-d', '--dir', dest='dir', help='Directory for the results (default is the current directory)', metavar='DIR', default='.')
    parser.add_option('-p', '--prefix', dest='prefix', help='Prefix for the result-filess (default is "")', metavar='PREFIX', default='')
    parser.add_option('-e', '--evalue', dest='evalue', help='E-value acceptance limit for BLAST hits (default 1e-5)', metavar='EVALUE', default='1e-5')
    parser.add_option('-r', '--rank', dest='rank', help='Number of BLAST hits to consider i.e. best K-term hit selected (default 5). Switch off by 0', metavar='RANK', default='5')
    parser.add_option('-g', '--group', dest='grpFile', help='A subset of gene identifiers to process (defaults to gene ids found in the blast file)', metavar='GROUP', default='')
    options, args = parser.parse_args()
    if options.bfile == '':
        print '\nBlast file not given:'
        parser.print_help()
        #print options
        sys.exit(1)
    return options


#################################################
def createDir(mydir):
    '''Creates a directory for the assembly if one does not exist yet.
    '''
    try:
        os.makedirs(mydir)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise

#################################################
def printLevel(level, dCategories, handle, handleFull):
    '''
    '''
    for cat in sorted(dCategories.iterkeys()):
        myCatTrCnt, myCatKtCnt, mySubCats, mySubSubCats, mySubSubSubCats = set(), set(), {}, {}, {}
        for subCat in dCategories[cat].iterkeys():
            mySubCatTrCnt, mySubCatKtCnt = set(), set()
            #mySubCats[subCat] = {}
            mySubSubCats[subCat] = {}
            mySubSubSubCats[subCat] = {}
            for subSubCat in dCategories[cat][subCat].iterkeys():
                mySubSubCatTrCnt, mySubSubCatKtCnt = set(), set()
                #mySubSubCats[subCat][subSubCat] = {}
                mySubSubSubCats[subCat][subSubCat] = {}
                for subSubSubCat in dCategories[cat][subCat][subSubCat].iterkeys():
                    ktCnt = set([item[0] for item in dCategories[cat][subCat][subSubCat][subSubSubCat]])
                    trCnt = set([item[1] for item in dCategories[cat][subCat][subSubCat][subSubSubCat]])
                    mySubSubCatKtCnt |= ktCnt
                    mySubSubCatTrCnt |= trCnt
                    mySubCatKtCnt |= ktCnt
                    mySubCatTrCnt |= trCnt
                    myCatKtCnt |= ktCnt
                    myCatTrCnt |= trCnt
                    try:
                        mySubSubSubCats[subCat][subSubCat][subSubSubCat].append((len(ktCnt), len(trCnt), subSubSubCat, trCnt))
                    except KeyError:
                        mySubSubSubCats[subCat][subSubCat][subSubSubCat] = []
                        mySubSubSubCats[subCat][subSubCat][subSubSubCat].append((len(ktCnt), len(trCnt), subSubSubCat, trCnt))
                try:
                    mySubSubCats[subCat][subSubCat].append((len(mySubSubCatKtCnt), len(mySubSubCatTrCnt), subSubCat, mySubSubCatTrCnt))
                except KeyError:
                    mySubSubCats[subCat][subSubCat] = []
                    mySubSubCats[subCat][subSubCat].append((len(mySubSubCatKtCnt), len(mySubSubCatTrCnt), subSubCat, mySubSubCatTrCnt))
            try:
                mySubCats[subCat].append((len(mySubCatKtCnt), len(mySubCatTrCnt), subCat, mySubCatTrCnt))
            except KeyError:
                mySubCats[subCat] = []
                mySubCats[subCat].append((len(mySubCatKtCnt), len(mySubCatTrCnt), subCat, mySubCatTrCnt))
        handle.write("%s\t%d\t%d\n" %(cat, len(myCatKtCnt), len(myCatTrCnt)))
        handleFull.write("%s\t%d\t%d\t%s\n" %(cat, len(myCatKtCnt), len(myCatTrCnt), '\t'.join(myCatTrCnt)))
        if level > 1:
            subCats = sorted([mySubCats[key] for key in mySubCats.iterkeys()])[::-1]
            for subCat in subCats:
                if len(subCat) > 1:
                    print >> std.err, "Fatal error in printing categories (1). Exiting..."
                    sys.exit(0)
                subCat = subCat[0]
                handle.write("\t%s\t%d\t%d\n" %(subCat[2], subCat[0], subCat[1]))
                handleFull.write("\t%s\t%d\t%d\t%s\n" %(subCat[2], subCat[0], subCat[1], '\t'.join(subCat[-1])))
                if level > 2:
                    subSubCats = sorted([mySubSubCats[subCat[2]][key] for key in mySubSubCats[subCat[2]].iterkeys()])[::-1]
                    for subSubCat in subSubCats:
                        if len(subSubCat) > 1:
                            print >> std.err, "Fatal error in printing categories (2). Exiting..."
                            sys.exit(0)
                        subSubCat = subSubCat[0]
                        if subSubCat[2] != '-':
                            handle.write("\t\t%s\t%d\t%d\n" %(subSubCat[2], subSubCat[0], subSubCat[1]))
                            handleFull.write("\t\t%s\t%d\t%d\t%s\n" %(subSubCat[2], subSubCat[0], subSubCat[1], '\t'.join(subSubCat[-1])))
                            if level > 3:
                                subSubSubCats = sorted([mySubSubSubCats[subCat[2]][subSubCat[2]][key] for key in mySubSubSubCats[subCat[2]][subSubCat[2]].iterkeys()])[::-1]
                                for subSubSubCat in subSubSubCats:
                                    if len(subSubSubCat) > 1:
                                        print >> std.err, "Fatal error in printing categories (3). Exiting..."
                                        sys.exit(0)
                                    subSubSubCat = subSubSubCat[0]
                                    if subSubSubCat[2] != '-':
                                        handle.write("\t\t\t%s\t%d\t%d\n" %(subSubSubCat[2], subSubSubCat[0], subSubSubCat[1]))
                                        handleFull.write("\t\t\t%s\t%d\t%d\t%s\n" %(subSubSubCat[2], subSubSubCat[0], subSubSubCat[1], '\t'.join(subSubSubCat[-1])))


#################################################
def printBrTerms(keggBrDetailsFile, handlePrint, handleFullPrint):
    '''
    ['Tcan_14523', 'br:ko01000', 'K16722', 'spindle-defective protein 2', 'cel:CELE_F32H2.3', 'Enzymes', '', '1e-14']
    ['Tcan_14522', 'br:ko03036', 'K07375', 'tubulin beta', 'spo:SPBC26H8.07c', 'Chromosome', 'Chromosome and associated proteins|Eukaryotic Type|Centrosomal proteins|Microtubules|Other tubulins| TUBB; tubulin beta', '4e-17']
    ['Tcan_14522', 'br:ko04147', 'K07375', 'tubulin beta', 'spo:SPBC26H8.07c', 'Exosome', 'Exosome|Exosomal proteins|Exosomal proteins of colorectal cancer cells| TUBB; tubulin beta', '4e-17']
    ['Tcan_14522', 'br:ko04812', 'K07375', 'tubulin beta', 'spo:SPBC26H8.07c', 'Cytoskeleton proteins', 'Cytoskeleton|Eukaryotic cytoskeleton proteins|Microtubules|Tubulins|Tubulins| TUBB; tubulin beta', '4e-17']
    ['Tcan_14522', 'br:ko01000', 'K07375', 'tubulin beta', 'spo:SPBC26H8.07c', 'Enzymes', '', '4e-17']
    ['Tcan_14522', 'br:ko00001', 'K07375', 'tubulin beta', 'spo:SPBC26H8.07c', 'KEGG Orthology (KO)', 'KO|Human Diseases|Infectious diseases|05130 Pathogenic Escherichia coli infection [PATH:ko05130]| TUBB; tubulin beta', '4e-17']
    '''
    # Read KEGG data if available
    keggBrites = []
    if os.path.exists(keggBrDetailsFile) == True:
        handle = open(keggBrDetailsFile, 'r')
        for line in handle:
            if line[0] != '#':
                items = line.strip().split('\t')
                #print items
                # Transcript id, brite term, kterm, annotation, category, detailed category, evalue
                category, subCategory = "", ""
                if len(items) < 7:
                    category, subCategory = items[5], "|-|-|-"
                else:
                    category, subCategory = items[5], items[6]
                    if subCategory == "": subCategory = "|-|-|-"
                keggBrites.append((items[0], items[1], items[2], items[3], category, subCategory, items[-1]))
    dCategories = {} # First level
    if len(keggBrites) > 0:
        for item in keggBrites:
            try:
                dCategories[item[4]].append((item[0], item[2], item[1], item[5]))
            except KeyError:
                dCategories[item[4]] = []
                dCategories[item[4]].append((item[0], item[2], item[1], item[5]))
            #print "%s\t%s\t%s" %("KEGGBR", item[0].replace('%', '_'), '\t'.join(item[1:]))

    dThird = {}
    for first in dCategories.iterkeys(): # Go through all the hits
        #print dCategories[first]
        for item in dCategories[first]:
            ktId, trId = item[1], item[0]
            try:
                dThird[first]
            except KeyError:
                dThird[first] = {}
            second = item[-1].split('|')[1].strip("</b>")
            try:
                dThird[first][second]
            except KeyError:
                dThird[first][second] = {}
            third = '-'
            try:
                third = item[-1].split('|')[2].strip("</b>")
            except IndexError:
                pass
            try:
                dThird[first][second][third]
                #dThird[first][second][third].append((ktId, trId))
            except KeyError:
                dThird[first][second][third] = {}
                #dThird[first][second][third] = []
                #dThird[first][second][third].append((ktId, trId))
            fourth = '-'
            try:
                fourth = item[-1].split('|')[3].strip("</b>")
            except IndexError:
                pass
            try:
                dThird[first][second][third][fourth].append((ktId, trId))
            except KeyError:
                dThird[first][second][third][fourth] = []
                dThird[first][second][third][fourth].append((ktId, trId))

    handlePrint.write("\n##########################################################\n")
    handlePrint.write("## Level1\tKtermCnt\tTranscriptCnt\n")
    handleFullPrint.write("\n##########################################################\n")
    handleFullPrint.write("## Level1\tKtermCnt\tTranscriptCnt\n")
    printLevel(1, dThird, handlePrint, handleFullPrint)

    handlePrint.write("\n##########################################################\n")
    handlePrint.write("## Level1\tLevel2\tKtermCnt\tTranscriptCnt\n")
    handleFullPrint.write("\n##########################################################\n")
    handleFullPrint.write("## Level1\tLevel2\tKtermCnt\tTranscriptCnt\n")
    printLevel(2, dThird, handlePrint, handleFullPrint)

    handlePrint.write("\n##########################################################\n")
    handlePrint.write("## Level1\tLevel2\tLevel3\tKtermCnt\tTranscriptCnt\n")
    handleFullPrint.write("\n##########################################################\n")
    handleFullPrint.write("## Level1\tLevel2\tLevel3\tKtermCnt\tTranscriptCnt\n")
    printLevel(3, dThird, handlePrint, handleFullPrint)

    handlePrint.write("\n##########################################################\n")
    handlePrint.write("## Level1\tLevel2\tLevel3\tLevel4\tKtermCnt\tTranscriptCnt\n")
    handleFullPrint.write("\n##########################################################\n")
    handleFullPrint.write("## Level1\tLevel2\tLevel3\tLevel4\tKtermCnt\tTranscriptCnt\n")
    printLevel(4, dThird, handlePrint, handleFullPrint)


#################################################
def readKeggBriteTrees(dataPath):
    '''
    '''
    '''
    #briteFiles = ["ko00001.keg", "ko00002.keg", "ko00194.keg", "ko00199.keg", \ 
    briteFiles = ["ko00194.keg", "ko00199.keg", \
                  "ko00535.keg", "ko01000.keg", "ko01001.keg", "ko01002.keg", \
                  "ko01003.keg", "ko01004.keg", "ko01005.keg", "ko01006.keg", \
                  "ko01007.keg", "ko01008.keg", "ko01020.keg", "ko02000.keg", \
                  "ko02001.keg", "ko02022.keg", "ko02035.keg", "ko02042.keg", \
                  "ko02044.keg", "ko03000.keg", "ko03009.keg", "ko03011.keg", \
                  "ko03012.keg", "ko03016.keg", "ko03021.keg", "ko03032.keg", \
                  "ko03036.keg", "ko03041.keg", "ko03051.keg", "ko03110.keg", \
                  "ko03310.keg", "ko03400.keg", "ko04030.keg", "ko04031.keg", \
                  "ko04040.keg", "ko04050.keg", "ko04052.keg", "ko04090.keg", \
                  "ko04091.keg", "ko04121.keg", "ko04131.keg", "ko04515.keg", \
                  "ko04516.keg", "ko04812.keg"]
    #briteFiles = ["ko01020.keg"]
    '''
    '''
    briteFiles = ["ko00001.keg", "ko00002.keg", "ko00003.keg", "ko00194.keg", "ko00199.keg", \
                  "ko00535.keg", "ko00536.keg", "ko01000.keg", "ko01001.keg", "ko01002.keg", \
                  "ko01003.keg", "ko01004.keg", "ko01005.keg", "ko01006.keg", "ko01007.keg", \
                  "ko01008.keg", "ko01020.keg", "ko02000.keg", "ko02001.keg", "ko02022.keg", \
                  "ko02035.keg", "ko02042.keg", "ko02044.keg", "ko03000.keg", "ko03009.keg", \
                  "ko03011.keg", "ko03012.keg", "ko03016.keg", "ko03021.keg", "ko03032.keg", \
                  "ko03036.keg", "ko03041.keg", "ko03051.keg", "ko03100.keg", "ko03110.keg", \
                  "ko03310.keg", "ko03400.keg", "ko04030.keg", "ko04031.keg", "ko04040.keg", \
                  "ko04050.keg", "ko04052.keg", "ko04090.keg", "ko04091.keg", "ko04121.keg", \
                  "ko04131.keg", "ko04147.keg", "ko04515.keg", "ko04516.keg", "ko04812.keg"]
    '''
    briteFiles = ["ko00001.keg", "ko00002.keg", "ko00003.keg", "ko00194.keg", "ko00199.keg", \
                  "ko00535.keg", "ko00536.keg", "ko01000.keg", "ko01001.keg", "ko01002.keg", \
                  "ko01003.keg", "ko01004.keg", "ko01005.keg", "ko01006.keg", "ko01007.keg", \
                  "ko01008.keg", "ko01009.keg", "ko01020.keg", "ko02000.keg", "ko02001.keg", \
                  "ko02022.keg", "ko02035.keg", "ko02042.keg", "ko02044.keg", "ko03000.keg", \
                  "ko03009.keg", "ko03011.keg", "ko03012.keg", "ko03016.keg", "ko03021.keg", \
                  "ko03029.keg", "ko03032.keg", "ko03036.keg", "ko03041.keg", "ko03051.keg", \
                  "ko03100.keg", "ko03110.keg", "ko03310.keg", "ko03400.keg", "ko04030.keg", \
                  "ko04031.keg", "ko04040.keg", "ko04050.keg", "ko04052.keg", "ko04090.keg", \
                  "ko04091.keg", "ko04121.keg", "ko04131.keg", "ko04147.keg", "ko04515.keg", \
                  "ko04516.keg", "ko04812.keg"]


    # Create trees
    keggTrees = []
    comments = {}
    for filename in briteFiles:
        key = "br:" + filename[:-4]
        keggTree = KeggTree(key)
        comments[key] = {}
        keggTree.parse(dataPath + filename)
        keggTrees.append(keggTree)
        #if filename == "ko01009.keg":
        #    print keggTree

    # Create comment lines for each K term
    for i in xrange(len(keggTrees)):
        stack = []
        keggTree = keggTrees[i]
        #print keggTree
        prevLevel = -1
        for node in keggTree.depthFirst():
            stack.append(node)
            myDiff = prevLevel - node.level
            prevLevel = node.level
            if myDiff >= 0:
                #print "Diff: %d" %myDiff
                for i in xrange(myDiff + 2):
                    stack.pop()
                stack.append(node)
                #print len(stack),
                #for item in stack:
                #    print '|' + item.label,
                #print ""
            #print "(%d %d)" %(prevLevel, node.level), 
            #print node.label
            #print len(stack),
            #for item in stack:
            #    print '|' + item.label,
            #print ""
            items = node.label.split(' ')
            #if items[0] == "K05694":
            #    print "Found: K05694"
            if len(items[0]) == 6 and items[0][0] == 'K': # Is a K-term?
                #print "### " + items[0],
                myStack = [node.label for node in stack]
                myStack[-1] = ' '.join(myStack[-1].split(' ')[1:])
                myStack = '|'.join(myStack)
                #print myStack
                comments[keggTree.name][items[0]] = myStack
                #stack.pop()
                #if items[0] == "K05694":
                #    print myStack
    return comments


#################################################
def pathwayKtFreqs(kterms, dPw):
    # Create Kterm and transcription frequencies per each BRITE term (44)
    dKt, myPaths = {}, {}
    for trId in kterms.iterkeys():
        dKt[kterms[trId]] = True
    for kId in dKt.iterkeys():
        try:
            pathIdList = dPw[kId]
            for pathId in pathIdList:
                try:
                    myPaths[pathId] += 1
                except KeyError:
                    myPaths[pathId] = 1
        except KeyError:
            #print kterms[trId]
            #print trId
            pass
    return myPaths


#################################################
def pathwayTrFreqs(kterms, dPw):
    # Create Kterm and transcription frequencies per each BRITE term (44)
    myPaths = {}
    myPathDetails = []
    for trId in kterms.iterkeys():
        try:
            pathIdList = dPw[kterms[trId]]
            for pathId in pathIdList:
                try:
                    myPaths[pathId] += 1
                except KeyError:
                    myPaths[pathId] = 1
                myPathDetails.append((trId, pathId, kterms[trId]))
        except KeyError:
            #print kterms[trId]
            #print trId
            pass
    #for key in sorted(myPaths.iterkeys()):
    #    print "%s: %d" %(key, myPaths[key])
    return myPaths, myPathDetails


#################################################
def briteKtFreqs(kterms, dBr):
    # Create Kterm and transcription frequencies per each BRITE term (44)
    dKt, myBrites = {}, {}
    for trId in kterms.iterkeys():
        dKt[kterms[trId]] = True
    for kId in dKt.iterkeys():
        try:
            briteIdList = dBr[kId]
            for briteId in briteIdList:
                try:
                    myBrites[briteId] += 1
                except KeyError:
                    myBrites[briteId] = 1
        except KeyError:
            print kterms[trId]
            print trId
    #for key in sorted(myBrites.iterkeys()):
    #    print "%s: %d" %(key, myBrites[key])
    return myBrites


#################################################
def briteTrFreqs(kterms, dBr):
    # Create Kterm and transcription frequencies per each BRITE term (44)
    myBrites = {}
    myBriteDetails = []
    for trId in kterms.iterkeys():
        try:
            briteIdList = dBr[kterms[trId]]
            for briteId in briteIdList:
                try:
                    myBrites[briteId] += 1
                except KeyError:
                    myBrites[briteId] = 1
                myBriteDetails.append((trId, briteId, kterms[trId]))
        except KeyError:
            print kterms[trId]
            print trId
    #for key in sorted(myBrites.iterkeys()):
    #    print "%s: %d" %(key, myBrites[key])
    return myBrites, myBriteDetails


#################################################
def calcFreqs(dAnno):
    '''
    '''
    trs, kts, paths = set(), set(), set()
    #for anno in dAnno.iterkeys():
    #    annoList = dAnno[anno]
    #for item in annoList:
    for item in dAnno:
        trs.add(item[0])
        paths.add(item[1])
        kts.add(item[2])
    return len(kts), len(trs)

#################################################
def createHierarchies(myTerms, myAnnos):
    '''
    '''
    level1, level2, level3 = {}, {}, {}
    maxLen = 0
    for item in myTerms:
        #print "###"
        #print item
        #print list(myAnnos[item[1]])[0]
        items = list(myAnnos[item[1]])[0].strip().split(';')
        #print items
        if len(items) > maxLen: maxLen = len(items)
        annoItems = []
        for anno in items:
            annoItems.append(anno.strip())
        #print item
        try:
            level1[annoItems[0]].append((item[0], item[1], item[2]))
        except KeyError:
            level1[annoItems[0]] = []
            level1[annoItems[0]].append((item[0], item[1], item[2]))
        try:
            #print annoItems
            level2['\t'.join(annoItems[0:2])].append((item[0], item[1], item[2]))
        except KeyError:
            level2['\t'.join(annoItems[0:2])] = []
            level2['\t'.join(annoItems[0:2])].append((item[0], item[1], item[2]))
        try:
            level3['\t'.join(annoItems[0:3])].append((item[0], item[1], item[2]))
        except KeyError:
            level3['\t'.join(annoItems[0:3])] = []
            level3['\t'.join(annoItems[0:3])].append((item[0], item[1], item[2]))
        # TrId, BriteId, Kterm, Annotation
    if maxLen > 3:
       print "Number of hierarchies (%d) was bigger than expected 3" %maxLen
    return level1, level2, level3

#################################################
def readAnnos(dataPath):
    '''
    '''
    inClass, inBrite = False, False
    limits = ["DBLINKS", "MODULE", "GENES", "DISEASE", "///", "ENTRY", "BRITE"] # DBLINKS not always defined, genes are next in line
    #limits = ["DBLINKS", "GENES", "///"] # DBLINKS not always defined, genes are next in line
    limitsBrite = ["DBLINKS", "MODULE", "GENES", "DISEASE", "///", "ENTRY"] # DBLINKS not always defined, genes are next in line
    #myKterm, myClass = None, None
    #unClassKterms = []
    myAnnos, myKtermAnnos = {}, {}
    handle = open(dataPath + "ko", 'r')
    geneId, geneStr, annos, anno = None, None, {}, ["", "", ""]
    for line in handle:
        #if line[0] == '\t': print line,
        if len(line) > 7:
            items = line.strip().split()

            if items[-1][0:4] == "[BR:":
                key = items[-1].strip('[]').lower()
                myItems = list(items)
                if myItems[0] == 'BRITE': myItems = myItems[1:]
                try:
                    myAnnos[key].add(' '.join(myItems[:-1]))
                except KeyError:
                    myAnnos[key] = set()
                    myAnnos[key].add(' '.join(myItems[:-1]))

            #print items

            if inBrite == True and items[0] not in limitsBrite:
                #print items
                if items[-1][0] == '[':
                    inBrite = False
                    #print "Brite inactive"
                    anno = ["", "", ""]
                else:
                    key = "path:ko%s" %items[0]
                    #print key
                    if key not in list(annos.iterkeys()):
                        #print "Checking"
                        if line[14] == ' ':
                            pass
                        elif line[13] == ' ':
                            anno[1] = ' '.join(items)
                        elif line[12] == ' ':
                            anno[0] = ' '.join(items)
                    else:
                        anno[2] =' '.join(items[1:])
                        #print "Third level"
                        try:
                            myAnnos[key].add(';'.join(anno))
                        except KeyError:
                            myAnnos[key] = set()
                            myAnnos[key].add(';'.join(anno))
                        if len(anno) != 3: print "%s\t%s" %(key, anno)
                        #print "%s\t%s" %(key, anno)
                    
            if items[0] == "BRITE" and inClass == True:
                #print "In Brite"
                if items[-1][0] == '[':
                    key = items[-1].strip('[]').lower()
                    #key = items[-1].strip('[').strip(']').lower()
                    if key == "br:ko00001": # Metabolism only read in
                        anno = ["", "", ""]
                        inBrite, inClass = True, False
                        #print "Brite active"

            if inClass == True and items[0] not in limits:
                #print items
                myAnno = items[1:]
                key = "path:%s" %items[0]
                try:
                    annos[key].add(' '.join(myAnno))
                except KeyError:
                    annos[key] = set()
                    annos[key].add(' '.join(myAnno))
                #if "path:ko03020" in list(annos.iterkeys()):
                #    print "Found"
                #    print annos
            if items[0] == "PATHWAY":
                inClass, annos = True, {}
                myAnno = items[2:]
                key = "path:%s" %items[1]
                try:
                    annos[key].add(' '.join(myAnno))
                except KeyError:
                    annos[key] = set()
                    annos[key].add(' '.join(myAnno))

            if items[0] == "ENTRY":
                if geneId != None:
                    myKtermAnnos[geneId] = geneStr
                    #print "%s %s" %(geneId, geneStr)
                geneId = items[1]
            if items[0] == "DEFINITION":
                #DEFINITION  alcohol dehydrogenase [EC:1.1.1.1]
                geneStr = ' '.join(items[1:])
                myKtermAnnos[geneId] = geneStr
                geneId = None
            if items[0] == "NAME":
                geneStr = ' '.join(items[1:])
                if geneStr == geneId: geneStr = "Unclassified"

            #if line[0] != '\t': # and line[0] != '\t':
            if items[0] in limitsBrite:
                inBrite = False
                #print "Brite inactive"
                #if items[0] in limits:
                #inClass = False
                #myClass = items[1]
                #unclassified = False
                #if items[1] == "Unclassified;": unclassified = True
                #if myKterm != None:
                #    if unclassified == True:
                #        unClassKterms.append(myKterm)
                #    myKterm = None
    handle.close()
    
    #myAnnos["br:ko01000"] = set(["Enzymes"])
    #myAnnos["br:ko00001"] = set(["Metabolism"])
    #myAnnos["br:ko00002"] = set(["Pathway"])
    myAnnos["br:ko02022"] = set(["Two component regulatory system"])
    myAnnos["br:ko02000"] = set(["Transport system"])
    myAnnos["br:ko00194"] = set(["Photosynthesis proteins"])
    myAnnos["br:ko03041"] = set(["Spliceosome"])
    myAnnos["br:ko03021"] = set(["Transcription machinery (RNA polymerase)"])
    myAnnos["br:ko03051"] = set(["Proteasome"])
    myAnnos["br:ko02044"] = set(["Secretion system (Protein transport system"])
    myAnnos["br:ko03011"] = set(["Ribosome"])
    myAnnos["br:ko03012"] = set(["Translation factors"])
    myAnnos["br:ko03400"] = set(["DNA repair and recombination proteins"])
    myAnnos["br:ko03032"] = set(["DNA replication proteins"])
    myAnnos["br:ko04121"] = set(["Ubiquitin system"])

    for anno in myAnnos.iterkeys():
        if len(myAnnos[anno]) > 1:
            print "More than one annotation per ko-term: %s" %anno
            print myAnnos[anno]
    return myAnnos, myKtermAnnos


#################################################
def main():
    '''
    '''
    opts = options()
    base = Base()

    if opts.dir != ".": createDir(opts.dir)

    rank = int(opts.rank)
    bfile = opts.bfile
    ofile = opts.bfile.split('/')[-1]
    dataPath = os.environ['KEGGDB']
    dataPath = "%s.new" %dataPath # Uses Dec 2013 version
    if dataPath[-1] != '/': dataPath += '/'

    # Read KEGG BRITE trees from .keg files
    kTermBriteComments = readKeggBriteTrees(dataPath)
    #print kTermBriteComments

    # Read in the gene identifiers if available
    geneIds = []
    if opts.grpFile != '':
        handle = open(opts.grpFile)
        for line in handle:
            geneIds.append(line.strip().strip('"'))
        handle.close()

    # K-term dictionaries and reverse dictionaries
    print >> sys.stderr, "## Reading in KEGG files ..."
    dEn, dBr, dPw = {}, {}, {} # K-term is key here
    rEn, rBr, rPw, rGn = {}, {}, {}, {}
    handle = open(dataPath + "ko_enzyme.list", 'r') # Enzymes: one to many
    for line in handle:
        items = line.strip().split()
        dEn[items[0][3:]] = items[1] # EC numbers
        rEn[items[1]] = items[0][3:] # K-terms
    handle.close()

    handle = open(dataPath + "ko_brite.list", 'r') # Brite terms: one to many
    for line in handle:
        items = line.strip().split()
        #if items[0][3:] == "K05694":
        #    print "Found K05694"
        try:
            dBr[items[0][3:]].add("br:" + items[1][3:]) # Brite ko file identifiers: key is kterm
        except KeyError:
            dBr[items[0][3:]] = set() # Brite ko file identifiers: key is kterm
            dBr[items[0][3:]].add("br:" + items[1][3:]) # Brite ko file identifiers: key is kterm
        #if items[0][3:] == "K05694":
        #    print dBr[items[0][3:]]
        rBr["br:" + items[1][3:]] = True # key is brite term
    handle.close()

    #for key in dBr.iterkeys():
    #    for item in dBr[key]:
    #        if item == "br:ko01008": print "Found 8!"        
    #        if item == "br:ko01009":
    #            print "Found 9!"
    #            break

    handle = open(dataPath + "ko_pathway.list", 'r') # Pathways: many to many
    for line in handle:
        items = line.strip().split()
        try:
            dPw[items[0][3:]].append("path:" + items[1][5:]) # Pathway ko numbers
        except KeyError:
            dPw[items[0][3:]] = []
            dPw[items[0][3:]].append("path:" + items[1][5:]) # Pathway ko numbers
        try:
            rPw["path:" + items[1][5:]].append(items[0][3:]) # K-terms
        except KeyError:
            rPw["path:" + items[1][5:]] = []
            rPw["path:" + items[1][5:]].append(items[0][3:]) # K-terms
    handle.close()

    for key in rPw.iterkeys():
        rPw[key] = set(rPw[key])

    handle = open(dataPath + "ko_genes.list", 'r') # Genes: many to one
    for line in handle:
        items = line.strip().split()
        #print items[1]
        rGn[items[1]] = items[0][3:] # K-terms
        try:
            dBr[items[0][3:]].add("br:ko01000") # Brite ko file identifiers: key is kterm
        except KeyError:
            dBr[items[0][3:]] = set() # Brite ko file identifiers: key is kterm
            dBr[items[0][3:]].add("br:ko01000")
        rBr["br:ko01000"] = True # key is brite term
    handle.close()

    # Read BLAST file
    print >> sys.stderr, "## Extracting K-terms from blast file ..."
    handle = base.ropen(bfile)
    prev = None
    cnt = 0
    skip = False
    dHits, kterms, myGeneIds, myEvalue = {}, {}, {}, {}
    for line in handle:
        items = line.strip().split('\t')
        if line[0] != '#' and float(items[10]) <= float(opts.evalue):
            cnt += 1
            trId = items[0]
            if trId != prev or prev == None:
                cnt = 1
                prev = trId
                skip = False
            if (cnt <= rank or rank == 0) and skip == False:
                geneId = items[1].split()[0]
                '''
                try:
                    kterms[trId]
                except KeyError: # Insert first hit only
                    try:
                        kterms[trId] = rGn[geneId]
                        myGeneIds[trId] = geneId
                        myEvalue[trId] = items[10] # Evalue
                    except KeyError:
                        #print geneId
                        pass
                #print line,
                '''
                try:
                    kterm = rGn[geneId]
                    try:
                        dHits[trId].append((float(items[10]), geneId, kterm, items[10]))
                        skip = True
                    except KeyError:
                        dHits[trId] = []
                        dHits[trId].append((float(items[10]), geneId, kterm, items[10]))
                        skip = True
                except KeyError:
                    pass
    base.rclose()

    # Select best K-term hit per transcript for further processing
    for trId in dHits.iterkeys():
        bestHit = sorted(dHits[trId])[0]
        #if len(dHits[trId]) > 1:
        #    print sorted(dHits[trId])
        #    print bestHit
        myGeneIds[trId] = bestHit[1]
        kterms[trId] = bestHit[2]
        myEvalue[trId] = bestHit[3]

    #handle = open("debug.file", 'w')
    #for myGeneId in myGeneIds.iterkeys():
    #    handle.write("%s\n" %myGeneIds[myGeneId])
    #handle.close()

    # Collect found transcripts in BLAST file and reduce to given geneId list
    if opts.grpFile != '':
        trIds = list(kterms.iterkeys())
        myTrs = list(set(geneIds) & set(trIds))
        print "Lengths of gene ids: ",
        print "blast file %d, given set %d, found KOs %d" %(len(trIds), len(geneIds), len(myTrs))
        kterms = dict((key, kterms[key]) for key in myTrs)

    kLen = float(len(set(kterms.values()))) # Number of kterms
    tLen = float(len(list(kterms.iterkeys()))) # Number of transcripts
    print "## Kt (%d), Tr(%d)" %(int(kLen), int(tLen))

    myAnnos, myKtermAnnos = readAnnos(dataPath)

    # Create Kterm and transcription frequencies per each brite term
    prefix = "%s/%s." %(opts.dir, ofile)
    if opts.prefix != "": prefix = "%s%s." %(prefix, opts.prefix) 
    print >> sys.stderr, "## KEGG BRITE terms frequencies to file %s ..." %("%skegg.brite.freqs" %prefix)
    handle = open("%skegg.brite.freqs" %prefix, 'w')
    handle.write("## BriteId\tKtermCnt\tTranscriptCnt\tKterm-%%\tTranscript-%%\tAnnotation\n")
    myBritesKt = briteKtFreqs(kterms, dBr)
    myBritesTr, myBriteTerms = briteTrFreqs(kterms, dBr)

    #for item in myBriteTerms:
    #    if item[1] == "br:ko01008": print "Found 8!"        
    #    if item[1] == "br:ko01009":
    #        print "Found 9!"
    #        break
    #    #print item        

    for key in sorted(myBritesTr.iterkeys()):
        #print myAnnos[key]
        kRatio = float(myBritesKt[key]) / kLen * 100.0
        tRatio = float(myBritesTr[key]) / tLen * 100.0
        try:
            myTxt = "%s\t%d\t%d\t%.2f\t%.2f\t%s" %(key, myBritesKt[key], myBritesTr[key], kRatio, tRatio, list(myAnnos[key])[0])
            #print myTxt
            handle.write("%s\n" %myTxt)
        except KeyError:
            print >> sys.stderr, "KeyError: %s not found!" %key
    handle.close()

    # Write detailed Kterms and brite terms per transcription
    print >> sys.stderr, "## KEGG BRITE terms defails to file %s ..." %("%skegg.brite.details" %prefix)
    handle = open("%skegg.brite.details" %prefix, 'w')
    handle.write("## TranscriptId\tBriteId\tKtermId\tKtermAnno\tGeneId\tAnnotation\n")
    for item in myBriteTerms:
        # TrId, BriteId, Kterm, Kterm annotation, Original GeneId, Annotation
        keggInfo = ""
        try:
            keggInfo = kTermBriteComments[item[1]][item[2]]
        except KeyError:
            pass
        try:
            myTxt = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" %(item[0], item[1], item[2], myKtermAnnos[item[2]], myGeneIds[item[0]], list(myAnnos[item[1]])[0], keggInfo, myEvalue[item[0]])
        except KeyError:
            print >> sys.stderr, "KeyError: %s not found!" %item[1]
        #print myTxt
        handle.write("%s\n" %myTxt)
    handle.close()

    # Create Kterm and transcription frequencies per each pathway term
    print >> sys.stderr, "## KEGG PATHWAY terms frequencies to file %s ..." %("%skegg.pathway.freqs" %prefix)
    handle = open("%skegg.pathway.freqs" %prefix, 'w')
    handle.write("## PathwayId\tFullCnt\tKtermCnt\tTranscriptCnt\tKterm-%%\tTranscript-%%\tAnnotation\n")
    myPathsKt = pathwayKtFreqs(kterms, dPw)
    myPathsTr, myPathTerms = pathwayTrFreqs(kterms, dPw)
    for key in sorted(myPathsTr.iterkeys()):
        kRatio = float(myPathsKt[key]) / kLen * 100.0
        tRatio = float(myPathsTr[key]) / tLen * 100.0
        myTxt = "%s\t%d\t%d\t%d\t%.2f\t%.2f\t%s" %(key, len(rPw[key]), myPathsKt[key], myPathsTr[key], kRatio, tRatio, list(myAnnos[key])[0])
        #print myTxt
        handle.write("%s\n" %myTxt)
    handle.close()

    # Write detailed Kterms and pathway terms per transcription
    print >> sys.stderr, "## KEGG PATHWAY terms defails to file %s ..." %("%skegg.pathway.details" %prefix)
    handle = open("%skegg.pathway.details" %prefix, 'w')
    handle.write("## TranscriptId\tPathwayId\tFullCnt\tKtermId\tKtermAnno\tGeneId\tAnnotation\n")
    for item in myPathTerms:
        # TrId, PathId, FullCnt, Kterm, Kterm Annotation, original gene id Annotation
        myTxt = "%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t" %(item[0], item[1], len(rPw[item[1]]), item[2], myKtermAnnos[item[2]], myGeneIds[item[0]], list(myAnnos[item[1]])[0], myEvalue[item[0]])
        #print myTxt
        handle.write("%s\n" %myTxt)
    handle.close()

    #############################################################################
    # Write hierarchical files
    print >> sys.stderr, "## KEGG PATHWAY hierarchy defails to file %s ..." %("%skegg.pathway.hierarchy" %prefix)
    level1, level2, level3 = createHierarchies(myPathTerms, myAnnos)
    list1 = []
    for key in level1.iterkeys():
        kLen, tLen = calcFreqs(level1[key])
        list1.append((key, kLen, tLen))
        #list1.append((key, len(level1[key])))
    list1 = sorted(list1, key = lambda s: s[1], reverse = True)
    #for item in list1:
    #    print "%s\t%d" %(item[0], item[1])
    list2 = []
    for key in level2.iterkeys():
        anno = key.split('\t')
        kLen, tLen = calcFreqs(level2[key])
        #print anno
        list2.append((anno[0], anno[1], kLen, tLen))
        #list2.append((anno[0], anno[1], len(level2[key])))
    list2 = sorted(list2, key = lambda s: s[2], reverse = True)
    list2 = sorted(list2, key = lambda s: s[0])
    #for item in list2:
    #    print "%s\t%s\t%d" %(item[0], item[1], item[2])
    list3 = []
    for key in level3.iterkeys():
        anno = key.split('\t')
        kLen, tLen = calcFreqs(level3[key])
        list3.append((anno[0], anno[1], anno[2], kLen, tLen))
        #list3.append((anno[0], anno[1], anno[2], len(level3[key])))
    list3 = sorted(list3, key = lambda s: s[3], reverse = True)
    list3 = sorted(list3, key = lambda s: s[1])
    list3 = sorted(list3, key = lambda s: s[0])
    #for item in list3:
    #    print "%s\t%s\t%s\t%d" %(item[0], item[1], item[2], item[3])

    handle = open("%skegg.pathway.hierarchy" %prefix, 'w')
    handle.write("## Level1\tKtermCnt\tTranscriptCnt\n")
    for item in list1:
        handle.write("%s\t%d\t%d\n" %(item[0], item[1], item[2]))
    handle.write("\n## Level1\tLevel2\tKtermCnt\tTranscriptCnt\n")
    for item in list2:
        handle.write("%s\t%s\t%d\t%d\n" %(item[0], item[1], item[2], item[3]))
    handle.write("\n## Level1\tLevel2\tLevel3\tKtermCnt\tTranscriptCnt\n")
    for item in list3:
        handle.write("%s\t%s\t%s\t%d\t%d\n" %(item[0], item[1], item[2], item[3], item[4]))
    handle.close()

    print >> sys.stderr, "## KEGG BRITE hierarchy defails to file %s ..." %("%skegg.brite.hierarchy" %prefix)
    handle = open("%skegg.brite.hierarchy" %prefix, 'w')
    handleFull = open("%skegg.brite.hierarchy.with.gene.ids" %prefix, 'w')
    printBrTerms("%skegg.brite.details" %prefix, handle, handleFull)
    handleFull.close()
    handle.close()

    #print myBritesTr
    #print myBritesKt
    #print len(myBritesTr)
    #print len(sorted(kterms.iterkeys()))
    #print sum(myBritesTr.values())
    #print sum(myBritesKt.values())


#################################################
if __name__ == "__main__":
    main()
