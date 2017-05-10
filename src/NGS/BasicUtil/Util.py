# -*- coding: UTF-8 -*-
import random,re
import pickle,os,configparser
import src.NGS.BasicUtil.DBManager as dbm

'''
Created on 2013-6-30

@author: rui
'''


currentpath=os.path.realpath(__file__)
currentpath[:currentpath.find("life/src")]+"life/com/config.properties"
cfparser = configparser.ConfigParser()
cfparser.read(currentpath[:currentpath.find("life/src")]+"life/com/config.properties")

ip=cfparser.get("mysqldatabase","ip")
print(currentpath,ip)#currentpath[:currentpath.find("life/src")]+"life/com/config.properties")
username=cfparser.get("mysqldatabase","username")
password=cfparser.get("mysqldatabase","password")
webdbname=cfparser.get("mysqldatabase","webdbname")
genomeinfodbname=cfparser.get("mysqldatabase","genomeinfodbname")
pekingduckchromtable=cfparser.get("mysqldatabase","pekingduckchromtable")
ghostdbname=cfparser.get("mysqldatabase","ghostdbname")
vcfdbname=cfparser.get("mysqldatabase","vcfdbname")
TranscriptGenetable=cfparser.get("mysqldatabase","TranscriptGenetable")
outgroupVCFBAMconfig_beijingref=cfparser.get("mysqldatabase","outgroupVCFBAMconfig_beijingref")
pathtoPython=cfparser.get("mysqldatabase", "pathtoPython")
beijingreffa=cfparser.get("mysqldatabase","beijingreffa")


def random_str(randomlength=8):
    a = list(string.ascii_letters)
    random.shuffle(a)
    return ''.join(a[:randomlength])
def alinmultPopSnpPos(vcfMaplist,jointmode="i"):
    """input:
    two map fomart like this [chrNo:[(pos,REF,ALT,INFO,FORMAT,sample,...),(pos,REF,ALT,INFO,FORMAT,sample,...),,,,,],{chrNo:[]},,,,,,]
    output:
    one map like this {chrNo:[(pos,REF,ALT,(INFO,FORMAT,sample,...),(INFO,FORMAT,sample,...)),(,,,(),()),,,,,],chrNo:[],,,}
                                            from pop1                        from pop2
    """
    multipleVcfMap={}
    if len(vcfMaplist)==1 or jointmode=="o" or jointmode=="l":

        for currentChrom in vcfMaplist[0].keys():
            multipleVcfMap[currentChrom]=[]
            for SNPrec in vcfMaplist[0][currentChrom]:
                posInPop1 = SNPrec[0]#;print(posInPop1,file=open("testpos_8rep.txt0",'a'))
                RefInPop1 = SNPrec[1]
                AltInPop1 = SNPrec[2]
                multipleVcfMap[currentChrom].append([posInPop1,RefInPop1,AltInPop1,SNPrec[3:]])
        if len(vcfMaplist)==1:
            return copy.deepcopy(multipleVcfMap)
        vcfMap_obj_idx=0
        for vcfMap in vcfMaplist[1:]:
            vcfMap_obj_idx+=1
            for currentChrom in vcfMap:
                for SNPrec in vcfMap[currentChrom]:
                    posInPop1 = SNPrec[0]#;print(posInPop1,file=open("testpos_8rep.txt"+str(vcfMap_obj_idx),'a'))
                    RefInPop1 = SNPrec[1]
                    AltInPop1 = SNPrec[2]
                    low=0;high=len(multipleVcfMap[currentChrom])-1
                    while low<=high:
                        mid = (low + high)>>1
                        if multipleVcfMap[currentChrom][mid][0]<posInPop1:
                            low=mid+1
                        elif multipleVcfMap[currentChrom][mid][0]>posInPop1:
                            high=mid-1
                        else:
                            if AltInPop1 == multipleVcfMap[currentChrom][mid][2]:#same alt alle
                                fillNoneNum=vcfMap_obj_idx-(len(multipleVcfMap[currentChrom][mid])-3)
                                for i in range(fillNoneNum):
                                    multipleVcfMap[currentChrom][mid].append(None)
                                multipleVcfMap[currentChrom][mid].append(SNPrec[3:])
                            break
                    else:
                        if jointmode=="o":
                            insertelem=[posInPop1,RefInPop1,AltInPop1]
                            for i in range(0,vcfMap_obj_idx):
                                insertelem.append(None)
                            insertelem.append(SNPrec[3:])
                            multipleVcfMap[currentChrom].insert(low,insertelem)
    #list(multipleVcfMap.keys())[0]==currentChrom
    #when a pos only exist in the former several pops,but not exist in the rear several pops,the loop block under are neccessary
        for REC_idx in range(0,len(multipleVcfMap[list(multipleVcfMap.keys())[0]])):
            #
            for i in range(len(vcfMaplist)+3-len(multipleVcfMap[list(multipleVcfMap.keys())[0]][REC_idx])):
                multipleVcfMap[list(multipleVcfMap.keys())[0]][REC_idx].append(None)

        return copy.deepcopy(multipleVcfMap)
    
    for currentChrom in vcfMaplist[0].keys():
#             self.FstMapByChrom[currentChrom] = []
        multipleVcfMap[currentChrom]=[]
        for SNPrec in vcfMaplist[0][currentChrom]:
            posInPop1 = SNPrec[0]
            RefInPop1 = SNPrec[1]
            AltInPop1 = SNPrec[2]
            elementToAppend=[posInPop1,RefInPop1,AltInPop1,SNPrec[3:]]
            if len(vcfMaplist)==1:
                multipleVcfMap[currentChrom].append(elementToAppend)
                continue
            for vcfMap_obj_idx in range(1,len(vcfMaplist[:])):
                vcfMap_obj=vcfMaplist[vcfMap_obj_idx]
                if currentChrom not in vcfMap_obj or len(vcfMap_obj[currentChrom])==0:
                    print("alinmultPopSnpPos",currentChrom,"didn't find in vcfMap2")
#                     if jointmode=="i":
                    break
#                     elif jointmode=="o":
#                         if vcfMap_obj_idx!=len(vcfMaplist)-1:
#                             elementToAppend.append(None)
#                         else:
#                             elementToAppend.append(None)
#                             multipleVcfMap[currentChrom].append(elementToAppend)
                low = 0
                high = len(vcfMap_obj[currentChrom]) - 1
                
                if re.search(r"[A-Za-z]+,[A-Za-z]+", AltInPop1) != None:  # multiple allels
                    continue
#                dp4 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", SNPrec[3])
#                 print(dp4.group(0))
                
                while low <= high:
                    mid = (low + high)>>1
                    if vcfMap_obj[currentChrom][mid][0]<posInPop1:
                        low=mid+1
                    elif vcfMap_obj[currentChrom][mid][0]>posInPop1:
                        high=mid-1
                    else:
                        if AltInPop1 == vcfMap_obj[currentChrom][mid][2]:#same alt alle
                            if vcfMap_obj_idx!=len(vcfMaplist)-1:
                                elementToAppend.append(vcfMap_obj[currentChrom][mid][3:])
                            elif vcfMap_obj_idx==len(vcfMaplist)-1:
                                elementToAppend.append(vcfMap_obj[currentChrom][mid][3:])
                                multipleVcfMap[currentChrom].append(elementToAppend)
#                         elif jointmode=="i":
#                         print("skip the different allele rec",currentChrom,posInPop1,AltInPop1,vcfMap_obj[currentChrom][mid][2])
#                             print(currentChrom,posInPop1,AltInPop1,vcfMap_obj[currentChrom][mid][2],"different alt allele,should skip this rec,but i have no time to improve this now")
#                         elif innerjoin_outjoin=="o":
#                             if vcfMap_obj_idx!=len(vcfMaplist)-1:
#                                 elementToAppend.append(None)
#                             elif vcfMap_obj_idx==len(vcfMaplist)-1:
#                                 elementToAppend.append(None)
#                                 multipleVcfMap[currentChrom].append(elementToAppend)
                        break
                else:
                    if jointmode=="i":
                        #ignore the rec
                        break
#                     elif jointmode=="o":
#                         if vcfMap_obj_idx!=len(vcfMaplist)-1:
#                             elementToAppend.append(None)
#                         elif vcfMap_obj_idx==len(vcfMaplist)-1:
#                             elementToAppend.append(None)
#                             multipleVcfMap[currentChrom].append(elementToAppend)                              
#                     print("snp not found in vcfMap2",SNPrec)
#                     self.doubleVcfMap[currentChrom].append(SNPrec+)
    return multipleVcfMap

class Window():
    def __init__(self):
        super().__init__()
        self.winValueL = []  # [(startPos,lastPos,value),(),,,,,,]
    def forPhastConsFormat(self, L, L_End_Pos, windowWidth, Caculator, winStart=0):
        """
        without overlap
        L=[startpos,endpos,value]
        """
        self.winValueL = []
        currentIdx = 0
        while currentIdx != len(L):
            if L[currentIdx][0] >= winStart and L[currentIdx][1] <= (winStart + windowWidth):
#                 print(L[currentIdx][1] - L[currentIdx][0])
                Caculator.process(L[currentIdx], L[currentIdx][1] - L[currentIdx][0])
                if L[currentIdx][1] == (winStart + windowWidth):
                    value = Caculator.getResult()
                    self.winValueL.append((winStart, winStart + windowWidth, value))
                    winStart += windowWidth
                
            elif L[currentIdx][0] > winStart and L[currentIdx][0] < (winStart + windowWidth) and L[currentIdx][1] > (winStart + windowWidth):
                print("2")
                frontPartPosNum = winStart + windowWidth - L[currentIdx][0]
                rearPartPosNum = L[currentIdx][1] - (winStart + windowWidth)
                Caculator.process(L[currentIdx], frontPartPosNum)
                value = Caculator.getResult()
                self.winValueL.append((winStart, winStart + windowWidth, value))
                winStart += windowWidth
                Caculator.process(L[currentIdx], rearPartPosNum)
            elif L[currentIdx][0] <= winStart and L[currentIdx][1] > winStart and L[currentIdx][1] < (winStart + windowWidth):
                print("3")
                rearPartPosNum = L[currentIdx][1] - (winStart + windowWidth)
                Caculator.process(L[currentIdx], frontPartPosNum)
            elif (winStart + windowWidth) <= L[currentIdx][0]:
                print("4")
                while  (winStart + windowWidth) <= L[currentIdx][0]:
                    print("Util", winStart + windowWidth , L[currentIdx][0])
                    self.winValueL.append((winStart, winStart + windowWidth, Caculator.getResult()))
                    winStart += windowWidth
            elif L[currentIdx][1] == winStart:
                self.winValueL.append((winStart, winStart + windowWidth, Caculator.getResult()))
                winStart += windowWidth
            currentIdx += 1
        else:
            self.winValueL.append((winStart, winStart + windowWidth, Caculator.getResult()))
                
    def slidWindowOverlap(self, L, L_End_Pos, windowWidth, slideSize, Caculator,L_Start_Pos=0):
        print("L_End_Pos",L_End_Pos,L_Start_Pos)
        """
        window slide from L_Start_Pos to L_End_Pos
        L = [(pos,p1,p2,p3,A_base_idx),(pos,"a,b","c,d","e,f",0),(pos,"a,b","c,d","e,f",1),....] for D-statistics wihtout "no covered"
        or 
        L = [(pos, REF, ALT, INFO,FORMAT,sampleslist),(pos, REF, ALT, INFO,FORMAT,sampleslist),(),...........] for any score need one vcf,
        or 
        L = [(pos,REF,ALT,(INFO,FORMAT,sampleslist),(INFO,FORMAT,sampleslist)),(pos,REF,ALT,(INFO,FORMAT,sampleslist),(INFO,FORMAT,sampleslist)),(),...........] for any score need one or multiple vcf,for example two vcf's compare,eg. fst, one or multiple vcf caculate hp
        like the two situation upside,return a value
        or
        L = [(pos,samples1dp,samples2dp,samples3dp,,,),(pos,samples1dp,samples2dp,samples3dp,,,),(),(),......]#in this situation ,value formation like this ([sample1_pecentage,sample2_pecentage,,,],[sample1_average_depth,sample2_average_depth,,,])
        
        """
#         del self.winValueL[:]
        self.winValueL = []  # notice here
        nextIdx = -1  # always be -1 if windowWidth == slideSize
        currentIdx = 0
        winStart = L_Start_Pos
        FoundNextIdx = False
        firstComeInWin = True
        notjustforsnp = True
        for findfirstidx in range(len(L)):
            if L[findfirstidx][0]>winStart:
                currentIdx=findfirstidx
                break
        while currentIdx != len(L):
#             print(L[currentIdx][0],L[currentIdx])
            if L[currentIdx][0] > winStart and  L[currentIdx][0] <= (winStart + windowWidth):
#                if notjustforsnp or (len(L[currentIdx][1])==1 and re.search(r'[^a-zA-Z]', L[currentIdx][2]) != None and len(L[currentIdx][2])==1):# it's not a snp? indel or cnv
                if firstComeInWin:
                    startPos = L[currentIdx][0]
                    firstComeInWin = False
                lastPos = L[currentIdx][0]
                Caculator.process(L[currentIdx])
                if FoundNextIdx == False and L[currentIdx][0] > (winStart + slideSize):  # always go to |currentIdx+=1|
                    nextIdx = currentIdx
                    FoundNextIdx = True
            else:
                noofsnps, value = Caculator.getResult()
                try:
                    self.winValueL.append((startPos, lastPos, noofsnps, value))
#                     print(startPos, lastPos, noofsnps, value)
                except:
                    print("no snp in first win", len(L), currentIdx, value, L[currentIdx])
                    self.winValueL.append((0, 0, noofsnps, value))
                    winStart += slideSize
                    continue
#                 self.winValueL.append((0, 0, value))
                winStart += slideSize
                firstComeInWin = True
                
                FoundNextIdx = False
                if nextIdx == -1:
                    if slideSize >= windowWidth:
                        while not (L[currentIdx][0] > winStart and  L[currentIdx][0] <= (winStart + windowWidth)) and L[currentIdx][0] > winStart + windowWidth:
                            winStart += slideSize
                            noofsnps, value = Caculator.getResult()
                            self.winValueL.append((0, 0, noofsnps, value))
                        if L[currentIdx][0] < winStart:
                            while currentIdx != len(L):
                                if L[currentIdx][0] > winStart and L[currentIdx][0] <= (winStart + windowWidth):
                                    break
                                elif L[currentIdx][0] < winStart:
                                    winStart += slideSize
                                    noofsnps, value = Caculator.getResult()
                                    self.winValueL.append((0, 0, noofsnps, value))
                                currentIdx += 1
#                             self.winValueL.append((0,0,'NA'))
#                             winStart += slideSize
                    continue  # go to |if L[currentIdx][0] > winStart and L[currentIdx][0] < (winStart + windowWidth):| in upside block
                else:
                    currentIdx = nextIdx
                    nextIdx = -1
                    continue
                
            currentIdx += 1
        else:
            noofsnps, value = Caculator.getResult()
            try :
                self.winValueL.append((startPos, lastPos, noofsnps, value))
#                 print(startPos, lastPos, noofsnps, value)
            except UnboundLocalError:
                self.winValueL.append((0, 0, noofsnps, value))
#             if nextIdx!=-1:
#                 currentIdx = nextIdx
#                 nextIdx = -1
#                 while currentIdx != len(L):
#                     lastPos = L[currentIdx][0]
#                     Caculator.process(L[currentIdx])
#                     currentIdx += 1
#                 else:
#                     noofsnps, value = Caculator.getResult()
#                     try:
#                         self.winValueL.append((startPos, lastPos, noofsnps, value))
#                     except:
#                         self.winValueL.append((0, 0, noofsnps, value))
#            
        
        n = int((L_End_Pos - (len(self.winValueL) * slideSize + windowWidth)) / slideSize) + 1
        for i in range(n):
            noofsnps, value = Caculator.getResult()
            self.winValueL.append((0, 0, noofsnps, value))