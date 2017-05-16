# -*- coding: UTF-8 -*-
'''
Created on 2013-7-2

@author: rui
'''

import re, copy, math, numpy, time,pysam
from itertools import combinations
from src.NGS.BasicUtil import VCFutil




class Caculator():
    def __init__(self):
        #every Caculator which need two or more vcf have the follow two variable
        self.vcfnamelist=[]
        self.vcfnameKEY_vcfobj_pyBAMfilesVALUE={}
    def process(self, T):
        pass
    def getResult(self):
        pass

class Caculate_df(Caculator):
    def __init__(self,pop1idxlist,pop2idxlist,filehanlder):
        super().__init__()
        self.pop1idxlist=pop1idxlist
        self.pop2idxlist=pop2idxlist
        self.COUNTED = 0# nooffixdiff
        self.COUNTEDadditional=[0,[0,0]]#noofheterozygosity (pop1recs,pop2recs)
        self.unsufficentfixediff=0
        self.pop1_indvds=None#=6#when pop1 is none at a pos,and no depth information
        self.pop2_indvds=None#=6
        self.filepoiner=filehanlder
        self.currentchrID=None
    def process(self,T):
        #no matter the present of the multiple allele or not,it's still work correct.
        pop1reffixed=0;pop1altfixed=0;pop2reffixed=0;pop2altfixed=0
        pop1het=0;pop2het=0
        pop1recs=0;pop2recs=0
        GT_idx = (re.split(":", T[3 + 0][1])).index("GT")  # gatk GT:AD:DP:GQ:PL
        for pop1idx in self.pop1idxlist:
            sample=T[3+0][2][pop1idx]
            if len(re.split(":", sample)) == 1:  # ./.
                continue
            pop1recs+=1
            if re.split(":", sample)[GT_idx]=="0/0":
                pop1reffixed+=1
            elif re.split(":", sample)[GT_idx]=="1/1":
                pop1altfixed+=1
            elif re.split(":", sample)[GT_idx]=="0/1":
                pop1het+=1
        for pop2idx in self.pop2idxlist:
            sample=T[3+0][2][pop2idx]
            if len(re.split(":", sample)) == 1:  # ./.
                continue
            pop2recs+=1
            if re.split(":", sample)[GT_idx]=="0/0":
                pop2reffixed+=1
            elif re.split(":", sample)[GT_idx]=="1/1":
                pop2altfixed+=1
            elif re.split(":", sample)[GT_idx]=="0/1":
                pop2het+=1
#                                  pop1 fixed as alt ,pop2 fixed as ref                                                      pop1 fixed as ref ,pop2 fixed as alt
#         print(self.pop2_indvds,pop2reffixed+pop2altfixed,pop2reffixed*pop2altfixed,self.pop1_indvds,pop1altfixed+pop1reffixed,pop1reffixed*pop1altfixed)
#         print(pop1reffixed,"==0",pop1altfixed,"==",self.pop1_indvds,pop2reffixed,"==",self.pop2_indvds,pop2altfixed,"==0",(pop1reffixed==0 and pop1altfixed==self.pop1_indvds and pop2reffixed==self.pop2_indvds and pop2altfixed==0 ),pop1reffixed,"==",self.pop1_indvds,pop1altfixed,"==0",pop2reffixed,"==0",pop2altfixed,"==",self.pop2_indvds,(pop1reffixed==self.pop1_indvds and pop1altfixed==0 and pop2reffixed==0 and pop2altfixed==self.pop2_indvds))
        if (pop1reffixed==0 and pop1altfixed==self.pop1_indvds and pop2reffixed==self.pop2_indvds and pop2altfixed==0 ) or (pop1reffixed==self.pop1_indvds and pop1altfixed==0 and pop2reffixed==0 and pop2altfixed==self.pop2_indvds):
            self.COUNTED+=1
            print(self.currentchrID,T,file=self.filepoiner)
#             print("bingle",file=open("tttttt.txt",'a'))
        elif (pop1reffixed*pop1altfixed!=0 or pop1het!=0) or (pop2reffixed*pop2altfixed!=0 or pop2het!=0 ) :
            self.COUNTEDadditional[0]+=1
        elif (pop1reffixed*pop1altfixed==0 and pop1het==0) and (pop2reffixed*pop2altfixed==0 and pop2het==0 ) :
            self.COUNTEDadditional[1][0]+=pop1recs/self.pop1_indvds;self.COUNTEDadditional[1][1]+=pop2recs/self.pop2_indvds;self.unsufficentfixediff+=1
    def getResult(self):
        if self.unsufficentfixediff!=0:
            pop1unsufficentfixed=self.COUNTEDadditional[1][0]/self.unsufficentfixediff
            pop2unsufficentfixed=self.COUNTEDadditional[1][1]/self.unsufficentfixediff
        else:
            pop1unsufficentfixed=0;pop2unsufficentfixed=0
        noofhet=self.COUNTEDadditional[0]
        nooffixediff=self.COUNTED
        self.COUNTED=0
        self.COUNTEDadditional=[0,[0,0]]
        self.unsufficentfixediff=0
        return [noofhet,(pop1unsufficentfixed,pop2unsufficentfixed)],nooffixediff #self.COUNTEDadditional,self.COUNTED
class Caculate_Hp_master_slave(Caculator):
    def __init__(self, listOftargetpopvcfconfig,outfileprewithpath, minsnps=10,depth=10):
        super().__init__()
        self.minsnps = minsnps
        self.depth=depth
        self.SeqMethodlist=[]
        self.vcfnamelist=[]
        self.vcfnameKEY_vcfobj_pyBAMfilesVALUE={}
        self.outputname=outfileprewithpath
        self.listOfpopvcfRecsmapByAChr=[]
        for vcfconfigfilename in listOftargetpopvcfconfig[:]:
            self.listOfpopvcfRecsmapByAChr.append({})
            vcfconfig=open(vcfconfigfilename,"r")
            for line in vcfconfig:
                vcffilename_obj=re.search(r"vcffilename=(.*)",line.strip())
                if vcffilename_obj!=None:
                    vcfname=vcffilename_obj.group(1).strip()
                    self.vcfnamelist.append(vcfname)
                    self.outputname+=("_"+re.split(r"\.",re.search(r"[^/]*$",vcfname).group(0))[0])[:3]
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(VCFutil.VCF_Data(vcfname))
                elif line.split():
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.Samfile(line.strip(),'rb'))
            vcfconfig.close()
            if re.search(r"indvd[^/]+",vcfname)!=None:
                self.SeqMethodlist.append("indvd")
    
            elif re.search(r"pool[^/]+",vcfname)!=None:
                self.SeqMethodlist.append("pool")
    
            else:
                print("vcfname must with 'pool' or 'indvd'")
                exit(-1) 
        print(self.SeqMethodlist)
        self.COUNTED = [0] * len(self.SeqMethodlist)
        self.CNMI = [0] * len(self.SeqMethodlist)
        self.CNMA = [0] * len(self.SeqMethodlist)
        self.sum_mean_2pq = 0  
    def process(self, T, seqerrorrate=0.008, mode=1):
        if len(T[1]) != len(T[2]) or len(T[2])!=1 or len(T[2])!=1:
            return
        for MethodToSeq_idx in range(len(self.SeqMethodlist)):
            MethofToSeq = self.SeqMethodlist[MethodToSeq_idx]
            if T[3 + MethodToSeq_idx] == None:
                continue
            if MethofToSeq == "pool":
                refdep = 0;altalleledep = 0
                AD_idx = (re.split(":", T[3 + MethodToSeq_idx][1])).index("AD")  # gatk GT:AD:DP:GQ:PL
                for sample in T[3 + MethodToSeq_idx][2]:
                    if len(re.split(":", sample)) == 1:  # ./.
                        continue
                    AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                    try :
                        refdep += int(AD_depth[0])*0.7
                        altalleledep += int(AD_depth[1])*0.7
                    except ValueError:
                        print(sample, end="|")
            elif MethofToSeq == "indvd":
                AF = float(re.search(r"AF=([\d\.e-]+);", T[3 + MethodToSeq_idx][0]).group(1))
                AN = int(re.search(r"AN=(\d+);", T[3 + MethodToSeq_idx][0]).group(1))
                AC = int(re.search(r"AC=(\d+);", T[3 + MethodToSeq_idx][0]).group(1))
                refdep = AN - AC
                altalleledep = AC
            if refdep <= seqerrorrate * (refdep + altalleledep):  # skip fixed as altallele ,ie refdep == 0
                continue
            if refdep + altalleledep < self.depth:
                continue
            self.COUNTED[MethodToSeq_idx] += 1
            if refdep < altalleledep:
                self.CNMI[MethodToSeq_idx] += refdep
                self.CNMA[MethodToSeq_idx] += altalleledep
            else:
                self.CNMA[MethodToSeq_idx] += refdep
                self.CNMI[MethodToSeq_idx] += altalleledep
    def getResult(self):
        HETEROZY = ['NA'] * len(self.SeqMethodlist)
        for MethodToSeq_idx in range(len(self.SeqMethodlist)):
            try:
                HETEROZY[MethodToSeq_idx] = self.CNMA[MethodToSeq_idx] * self.CNMI[MethodToSeq_idx] * 2 / ((self.CNMA[MethodToSeq_idx] + self.CNMI[MethodToSeq_idx]) ** 2)
            except ZeroDivisionError:
                # print("the Heterozigosity value of currentwindow is dividsion by zero,so set it to be NA")
                HETEROZY[MethodToSeq_idx] = 'NA'
        het_count = 0;het_sum = 0;pop_idx=0
        for pop_idx in range(len(HETEROZY)-1,-1,-1) :
            if HETEROZY[pop_idx] != 'NA' and self.COUNTED[pop_idx]>=self.minsnps:
                het_count += 1
                het_sum += HETEROZY[pop_idx]
            else:
                self.COUNTED.pop(pop_idx)
#                 print(HETEROZY[pop_idx],end="\t")
#         print()
        if self.COUNTED==[]:
            noofsnpcount=0
        else:
            noofsnpcount= min(self.COUNTED)
        if het_count == 0 :
            HETEROZY_toreturn = 'NA'
        else:
            HETEROZY_toreturn = het_sum / het_count
        self.COUNTED = [0] * len(self.SeqMethodlist)
        self.CNMA = [0] * len(self.SeqMethodlist)
        self.CNMI = [0] * len(self.SeqMethodlist)
        return noofsnpcount, HETEROZY_toreturn

class Caculate_S_ObsExp_difference(Caculator):
    def __init__(self,mindepthtojudefixed,listOftargetpopvcfconfig,listOfrefpopvcffileconfig,dbvariantstoolstojudgeancestral,toplevelTablejudgeancestralname,outfileprewithpath):
        super().__init__()
        self.dbvariantstoolstojudgeancestral=dbvariantstoolstojudgeancestral
        self.topleveltablejudgeancestralname=toplevelTablejudgeancestralname
        self.MethodToSeqpoplist=[]
        self.mindepthtojudefixed=20
        self.flankseqfafile=None
        self.N_of_targetpop=len(listOftargetpopvcfconfig)
        self.N_of_refpop=len(listOfrefpopvcffileconfig)
        self.listOfpopvcfRecsmapByAChr=[]
        self.vcfnamelist=[]
        self.vcfnameKEY_vcfobj_pyBAMfilesVALUE={}
        self.outputname=outfileprewithpath
        for vcfconfigfilename in listOftargetpopvcfconfig[:]+listOfrefpopvcffileconfig[:]:
            self.listOfpopvcfRecsmapByAChr.append({})
            vcfconfig=open(vcfconfigfilename,"r")
            for line in vcfconfig:
                vcffilename_obj=re.search(r"vcffilename=(.*)",line.strip())
                if vcffilename_obj!=None:
                    vcfname=vcffilename_obj.group(1).strip()
                    self.vcfnamelist.append(vcfname)
                    self.outputname+=("_"+re.split(r"\.",re.search(r"[^/]*$",vcfname).group(0))[0])[:3]
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(VCFutil.VCF_Data(vcfname))
                elif line.split():
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.Samfile(line.strip(),'rb'))
            vcfconfig.close()
            if re.search(r"indvd[^/]+",vcfname)!=None:
                self.MethodToSeqpoplist.append("indvd")
    
            elif re.search(r"pool[^/]+",vcfname)!=None:
                self.MethodToSeqpoplist.append("pool")
    
            else:
                print("vcfname must with 'pool' or 'indvd'")
                exit(-1)   

        self.currentchrID=None
        self.COUNT=0
        self.obsseq=[]
        self.CEXP=0
        self.CfixedDerived=0
        self.freq_xaxisKEY_yaxisVALUERelation=None
        self.minsnps=10
        self.alignedSNP_absentinfo={}#{chrNo:[[pos1,REF,ALT,(),(),(),,],[pos1,REF,ALT,(),(),(),,],[],,]}
        self.dynamicIU_toptable_obj=None
    def __del__(self):
        self.flankseqfafile.close()
    def process(self,T):
        """T=[pos,ref,alt,pop1,pop2,.....,popn]
        in ordered as the self.vcfnamelise
        """
        if len(T[1]) != len(T[2]) or len(T[2])!=1  or len(T[2])!=1:
            return
        snp=self.dbvariantstoolstojudgeancestral.operateDB("select","select * from "+self.topleveltablejudgeancestralname+" where chrID='"+self.currentchrID+"' and snp_pos='"+str(T[0])+"'")
        if not snp or snp[0][9]==None or snp[0][5]==None:#needed info of SNPs are absent, snp[0][9] and snp[0][5] are dependent on the fellow code segment
#             print(self.currentchrID,T,"snp not find,skip")
#             print("append in to alignedSNP_absentinfo",snp,T)
            self.alignedSNP_absentinfo[self.currentchrID].append(T)
            return
        else:
            A_base_idx=100
            fanyadepthlist=re.split(r",",snp[0][9])
            if len(fanyadepthlist)==2 and int(fanyadepthlist[1]) >=self.mindepthtojudefixed and fanyadepthlist[0].strip()=="0":
                A_base_idx=1
            elif len(fanyadepthlist)==2 and int(fanyadepthlist[0])>=self.mindepthtojudefixed and fanyadepthlist[1].strip()=="0":
                A_base_idx=0
            else:
#                 print("skip snp",snp[0][1],snp[0][7],snp[0][9],snp[0][11],snp[0][13])
                return

        ancestrallcontext=snp[0][5].strip()[0].upper()+snp[0][3+A_base_idx].strip().upper()+snp[0][5].strip()[2].upper()
        if "CG" in ancestrallcontext or "GC" in ancestrallcontext:
#             print("skip CG site",ancestrallcontext)
            return
        ##########x-axis
        countedAF=0;target_DAF_sum=0
        for tpopidx in range(3,self.N_of_targetpop+3):
            if T[tpopidx]==None:
                if len(self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[tpopidx-3]])==1:
#                     print("skip this pos",T)
                    continue
                else:
#                     depth_linelist=self.depthobjlist[tpopidx-3].getdepthByPos_optimized(self.currentchrID,T[0])
                    sum_depth=0
                    for samfile in self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[tpopidx-3]][1:]:
                        ACGTdep=samfile.count_coverage(self.currentchrID,T[0]-1,T[0])
                        for dep in ACGTdep:
                            sum_depth+=dep[0]
#                     for idx in self.species_idx_list[tpopidx-3][:]:
#                         sum_depth+=int(depth_linelist[idx])
                    if sum_depth>self.mindepthtojudefixed:
                        AF=0
                    else:
#                         print(sum_depth,"low coverage skip")
                        continue
            else:
                if self.MethodToSeqpoplist[tpopidx-3]=="indvd":
                    AF=float(re.search(r"AF=([\d\.e-]+);", T[tpopidx][0]).group(1))
                    AN = float(re.search(r"AN=([\d]+);", T[tpopidx][0]).group(1))
                    if AN<5:
                        continue
                elif self.MethodToSeqpoplist[tpopidx-3]=="pool":
                    refdep = 0;altalleledep = 0
                    AD_idx = (re.split(":", T[tpopidx][1])).index("AD")
                    for sample in T[tpopidx][2]:
                        if len(re.split(":", sample)) == 1:  # ./.
                            continue
                        AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                        try :
                            refdep += int(AD_depth[0])
                            altalleledep += int(AD_depth[1])
                        except ValueError:
                            print(sample, end="|")
                    if (refdep==altalleledep and altalleledep==0) or altalleledep+ refdep<10:
                        continue
                    AF=altalleledep/(altalleledep+refdep)
            if A_base_idx==0:
                DAF=1-AF
            elif A_base_idx==1:
                DAF=AF
            target_DAF_sum+=DAF;countedAF+=1
        if  countedAF==0 :#or target_DAF_sum/countedAF==0:
#             print("skip this snp,because it fiexd as ancestral or no covered in this pos in target pops",T,snp)
            return
        target_DAF=target_DAF_sum/countedAF
        #########y-axis
        countedAF=0;rer_DAF_sum=0
        for rpopidx in range(3+self.N_of_targetpop,self.N_of_refpop+self.N_of_targetpop+3):
            if T[rpopidx]==None:
                if len(self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[rpopidx-3]])==1:
#                     print("skip this snp",T)
                    continue
                else:
#                     depth_linelist=self.depthobjlist[rpopidx-3-self.N_of_targetpop].getdepthByPos_optimized(self.currentchrID,T[0])
                    sum_depth=0
                    for samfile in self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[rpopidx-3]][1:]:
                        ACGTdep=samfile.count_coverage(self.currentchrID,T[0]-1,T[0])
                        for dep in ACGTdep:
                            sum_depth+=dep[0]
#                     for idx in self.species_idx_list[rpopidx-3-self.N_of_targetpop][:]:
#                         sum_depth+=int(depth_linelist[idx])
                    if sum_depth>self.mindepthtojudefixed:
                        AF=0
                    else:
#                         print(sum_depth,"low depth skip")
                        continue
            else:
                if self.MethodToSeqpoplist[rpopidx-3]=="indvd":
                    AF=float(re.search(r"AF=([\d\.]+);", T[rpopidx][0]).group(1))
                elif self.MethodToSeqpoplist[rpopidx-3]=="pool":
                    refdep = 0;altalleledep = 0
                    AD_idx = (re.split(":", T[rpopidx][1])).index("AD")
                    for sample in T[rpopidx][2]:
                        if len(re.split(":",sample))==1:
                            continue
                        AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                        try :
                            refdep += int(AD_depth[0])
                            altalleledep += int(AD_depth[1])
                        except ValueError:
                            print(sample, end="|")
                    if refdep==altalleledep and altalleledep==0:
                        continue
                    AF=altalleledep/(altalleledep+refdep)
                if A_base_idx==0:
                    DAF=1-AF
                elif A_base_idx==1:
                    DAF=AF
                rer_DAF_sum+=DAF;countedAF+=1
        if  countedAF==0 or rer_DAF_sum==0:
#             print("skip this snp,because it  no covered in this pos in ref pops",T,snp)
            return
        for a,b in sorted(self.freq_xaxisKEY_yaxisVALUERelation.keys()):
            if target_DAF>a and target_DAF<=b:
                self.CEXP+=self.freq_xaxisKEY_yaxisVALUERelation[(a,b)]
                break
        self.obsseq.append(rer_DAF_sum/countedAF)
        self.COUNT+=1
        if rer_DAF_sum/countedAF==1:
            self.CfixedDerived+=1
    def getResult(self):
#         for chrom in self.alignedSNP_absentinfo.keys():
        if self.alignedSNP_absentinfo[self.currentchrID]!=[]:
            self.dynamicIU_toptable_obj.insertorUpdatetopleveltable(self.alignedSNP_absentinfo,self.flankseqfafile,50)
            No_Of_snpT=len(self.alignedSNP_absentinfo[self.currentchrID])
            for whatever in range(No_Of_snpT):
                snpT=self.alignedSNP_absentinfo[self.currentchrID].pop(0)
#                     print(snpT,len(self.alignedSNP_absentinfo[chrom]))
                print(snpT[0:3],"process",len(self.alignedSNP_absentinfo[self.currentchrID]),self.CEXP,self.COUNT)
                self.process(snpT)
            print("length of snpT",len(self.alignedSNP_absentinfo[self.currentchrID]))      
        S1="NA"
        S2="NA"
        try:
            S1=math.log(numpy.sum(self.obsseq)/self.CEXP)
            S2=0#(numpy.sum(self.obsseq)-self.CEXP)/numpy.std(self.obsseq,ddof=1)
        except:
            S1="NA"
            S2="NA"
        noofsnp=self.COUNT
        self.COUNT=0
        self.CEXP=0
        self.obsseq=[]
        self.CfixedDerived=0

                
        if S1=="NA" and S2=="NA" or noofsnp<self.minsnps:
            return noofsnp,"NA"
        return noofsnp,[S1,S2]
class Caculate_pairFst(Caculator):
    def __init__(self,mindepthtojudefixed,listOftargetpopvcfconfig,listOfrefpopvcffileconfig,outfileprewithpath,minsnps=10):
        super().__init__()
        self.minsnps=minsnps
        self.considerfixdiffinfst=False
        self.MethodToSeqpoplist=[]
        self.mindepthtojudefixed=mindepthtojudefixed
        self.N_of_targetpop=len(listOftargetpopvcfconfig)
        self.N_of_refpop=len(listOfrefpopvcffileconfig)
        self.outputname=outfileprewithpath
        self.vcfnamelist=[]#target ref
        self.listOfpopvcfRecsmapByAChr=[]
        self.vcfnameKEY_vcfobj_pyBAMfilesVALUE={}
        for vcfconfigfilename in listOftargetpopvcfconfig[:]+listOfrefpopvcffileconfig[:]:
            self.listOfpopvcfRecsmapByAChr.append({})
            vcfconfig=open(vcfconfigfilename,"r")
            for line in vcfconfig:
                vcffilename_obj=re.search(r"vcffilename=(.*)",line.strip())
                if vcffilename_obj!=None:
                    vcfname=vcffilename_obj.group(1).strip()
                    self.vcfnamelist.append(vcfname)
                    self.outputname+=("_"+re.split(r"\.",re.search(r"[^/]*$",vcfname).group(0))[0])[:3]
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(VCFutil.VCF_Data(vcfname))
                elif line.split():
                    self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.Samfile(line.strip(),'rb'))
            vcfconfig.close()
            if re.search(r"indvd[^/]+",vcfname)!=None:
                self.MethodToSeqpoplist.append("indvd")
    
            elif re.search(r"pool[^/]+",vcfname)!=None:
                self.MethodToSeqpoplist.append("pool")
    
            else:
                print("vcfname must with 'pool' or 'indvd'")
                exit(-1)   
        self.currentchrID=None
        self.COUNT=[(0,0)]*(self.N_of_refpop*self.N_of_targetpop)
        self.CNK=[0]*(self.N_of_refpop*self.N_of_targetpop)
        self.CDK=[0]*(self.N_of_refpop*self.N_of_targetpop)
        self.CfixedDerived=0
        self.freq_xaxisKEY_yaxisVALUERelation=None
        self.minsnps=10
    def process(self,T):
        """T=[pos,ref,alt,pop1,pop2,.....,popn]
        T is in the order as self.vcfnamelist
        """
        if len(T[1]) != len(T[2]) or len(T[2])!=1  or len(T[2])!=1:
            return
#         snp=self.dbvariantstoolstojudgeancestral.operateDB("select","select * from "+self.topleveltablejudgeancestralname+" where chrID='"+self.currentchrID+"' and snp_pos='"+str(T[0])+"'")
        for TT_idx in range(3,3+self.N_of_targetpop):
            TT=T[TT_idx]
            refdep_1=0;altalleledep_1=0
            if self.MethodToSeqpoplist[TT_idx-3]=="pool":
                if TT==None:
                    if len(self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[TT_idx-3]])==1:
                        print("skip this pos",T)
                        continue
                    else:
                        sum_depth=0
                        for samfile in self.vcfnameKEY_vcfobj_pyBAMfilesVALUE[self.vcfnamelist[TT_idx-3]][1:]:
                            ACGTdep=samfile.count_coverage(self.currentchrID,T[0]-1,T[0])
                            for dep in ACGTdep:
                                sum_depth+=dep[0]
                        if sum_depth>self.mindepthtojudefixed:
                            refdep_1=sum_depth;altalleledep_1=0
                        else:
                            continue
                else:
                    AD_idx_1 = (re.split(":", TT[1])).index("AD")  # gatk GT:AD:DP:GQ:PL
            for RT_idx in range(3+self.N_of_targetpop,len(T)):
                RT=T[RT_idx]
                refdep_2=0;altalleledep_2=0
                
        ##########x-axis
        countedAF=0;target_DAF_sum=0
        for tpopidx in range(3,self.N_of_targetpop+3):
            if T[tpopidx]==None:
                if self.depthobjlist==[]:
                    print("skip this pos",T)
                    continue
                else:
                    depth_linelist=self.depthobjlist[tpopidx-3].getdepthByPos_optimized(self.currentchrID,T[0])
                    sum_depth=0
                    for idx in self.species_idx_list[tpopidx-3][:]:
                        sum_depth+=int(depth_linelist[idx])
                    if sum_depth>self.mindepthtojudefixed:
                        AF=0
                    else:
                        continue
            else:
                if self.MethodToSeqpoplist[tpopidx-3]=="indvd":
                    AF=float(re.search(r"AF=([\d\.e-]+);", T[tpopidx][0]).group(1))
                    AN = float(re.search(r"AN=([\d]+);", T[tpopidx][0]).group(1))
                    if AN<5:
                        continue
                elif self.MethodToSeqpoplist[tpopidx-3]=="pool":
                    refdep = 0;altalleledep = 0
                    AD_idx = (re.split(":", T[tpopidx][1])).index("AD")
                    for sample in T[tpopidx][2]:
                        if len(re.split(":", sample)) == 1:  # ./.
                            continue
                        AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                        try :
                            refdep += int(AD_depth[0])
                            altalleledep += int(AD_depth[1])
                        except ValueError:
                            print(sample, end="|")
                    if (refdep==altalleledep and altalleledep==0) or altalleledep+ refdep<10:
                        continue
                    AF=altalleledep/(altalleledep+refdep)
            if A_base_idx==0:
                DAF=1-AF
            elif A_base_idx==1:
                DAF=AF
            target_DAF_sum+=DAF;countedAF+=1
        if  countedAF==0 :#or target_DAF_sum/countedAF==0:
#             print("skip this snp,because it fiexd as ancestral or no covered in this pos in target pops",T,snp)
            return
        target_DAF=target_DAF_sum/countedAF
        #########y-axis
        countedAF=0;rer_DAF_sum=0
        for rpopidx in range(3+self.N_of_targetpop,self.N_of_refpop+self.N_of_targetpop+3):
            if T[rpopidx]==None:
                if self.depthobjlist==[]:
                    print("skip this snp",T)
                    continue
                else:
                    depth_linelist=self.depthobjlist[rpopidx-3-self.N_of_targetpop].getdepthByPos_optimized(self.currentchrID,T[0])
                    sum_depth=0
                    for idx in self.species_idx_list[rpopidx-3-self.N_of_targetpop][:]:
                        sum_depth+=int(depth_linelist[idx])
                    if sum_depth>self.mindepthtojudefixed:
                        AF=0
                    else:
                        continue
            else:
                if self.MethodToSeqpoplist[rpopidx-3-self.N_of_targetpop]=="indvd":
                    AF=float(re.search(r"AF=([\d\.]+);", T[rpopidx][0]).group(1))
                elif self.MethodToSeqpoplist[rpopidx-3-self.N_of_targetpop]=="pool":
                    refdep = 0;altalleledep = 0
                    AD_idx = (re.split(":", T[rpopidx][1])).index("AD")
                    for sample in T[rpopidx][2]:
                        if len(re.split(":",sample))==1:
                            continue
                        AD_depth = re.split(",", re.split(":", sample)[AD_idx])
                        try :
                            refdep += int(AD_depth[0])
                            altalleledep += int(AD_depth[1])
                        except ValueError:
                            print(sample, end="|")
                    if refdep==altalleledep and altalleledep==0:
                        continue
                    AF=altalleledep/(altalleledep+refdep)
                if A_base_idx==0:
                    DAF=1-AF
                elif A_base_idx==1:
                    DAF=AF
                rer_DAF_sum+=DAF;countedAF+=1
        if  countedAF==0 or rer_DAF_sum==0:
#             print("skip this snp,because it  no covered in this pos in ref pops",T,snp)
            return
        for a,b in sorted(self.freq_xaxisKEY_yaxisVALUERelation.keys()):
            if target_DAF>a and target_DAF<=b:
                self.CEXP+=self.freq_xaxisKEY_yaxisVALUERelation[(a,b)]
                break
        self.obsseq.append(rer_DAF_sum/countedAF)
        self.COUNT+=1
        if rer_DAF_sum/countedAF==1:
            self.CfixedDerived+=1
    def getResult(self):
        S1="NA"
        S2="NA"
        try:
            S1=math.log(numpy.sum(self.obsseq)/self.CEXP)
            S2=0#(numpy.sum(self.obsseq)-self.CEXP)/numpy.std(self.obsseq,ddof=1)
        except:
            S1="NA"
            S2="NA"
        noofsnp=self.COUNT
        self.COUNT=0
        self.CEXP=0
        self.obsseq=[]
        self.CfixedDerived=0
        if S1=="NA" and S2=="NA" or noofsnp<self.minsnps:
            return noofsnp,"NA"
        return noofsnp,[S1,S2]

