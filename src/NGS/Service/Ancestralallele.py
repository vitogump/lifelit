'''
Created on 2014-11-30

@author: liurui
'''
import copy, time, pysam
import os, numpy, sys, re
import pickle


from NGS.BasicUtil import Util
import NGS.BasicUtil.DBManager as dbm


mindepthforJudgefixref=10
class dynamicInsertUpdateAncestralContext():#here use the fast edition
    def __init__(self,dbvariantstools,refseqfa,toplevelsnptablename="mspsgjlksy10pop_toplevel_pekingduckref_new",insertauthority=True,changetableauthority=True):
        self.refseqfahandler=open(refseqfa,'r')
        self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE={}#{vcfname:[vcfobj,pysam.samfile1,pysam.samfile2,,,,],,,,}
        self.mapofSNPrecsForeachVCFpop_mapBYchrom={}
        try:
            self.refseqfafasteridx=pickle.load(open(refseqfa+".myfasteridx",'rb'))
        except IOError:
            Util.generateFasterRefIndex(refseqfa, refseqfa+".myfasteridx")
            self.refseqfafasteridx=pickle.load(open(refseqfa+".myfasteridx",'rb'))
#         self.vcfnameKEY_vcfobjVALUE={}#{vcfname:vcfobj,,,,}
        self.currentchrLen=None
        self.toplevelsnptablename=toplevelsnptablename
        self.dbvariantstools=dbvariantstools
        self.toplevelsnptable_titlelist=[a[0].strip() for a in self.dbvariantstools.operateDB("select", "select column_name  from information_schema.columns where table_schema='" + Util.vcfdbname + "' and table_name='" + toplevelsnptablename + "'")]
        outgroupBAMconfigs=re.split(r";",Util.outgroupVCFBAMconfig_beijingref)
        for configfile in outgroupBAMconfigs:
            fp=open(configfile,'r')
            for line in fp:
                vcffilename_obj=re.search(r"vcffilename=(.*)",line.strip())
                if vcffilename_obj!=None:
                    vcfname=vcffilename_obj.group(1).strip()
                    archicpop_colname=re.search(r'[^/]*$',vcfname).group(0)
                    archicpop_colname=re.sub(r"[^\w^\d]","_",archicpop_colname)
                    print("dynamicInsertUpdateAncestralContext",vcfname,archicpop_colname)
                    for outgroupname_ALT in self.toplevelsnptable_titlelist[6::2]:
                        if archicpop_colname+"_alt"==outgroupname_ALT:
                            break
                    else:
                        print(archicpop_colname+"_alt does not exist in topleveltable" )
                        if changetableauthority:
                            self.dbvariant.operateDB("callproc", "mysql_sp_add_column", data=(Util.vcfdbname, toplevelsnptablename, archicpop_colname+"_alt", "char(128)", "default null"))
                            self.dbvariant.operateDB("callproc", "mysql_sp_add_column", data=(Util.vcfdbname, toplevelsnptablename, archicpop_colname+"_dep", "char(128)", "default null"))
                        else:
                            pass
                    #always a
                    self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
                    self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(VCFutil.VCF_Data(vcffilename_obj.group(1).strip()))
                elif line.split():
                    print("dynamicInsertUpdateAncestralContext",line.strip())
                    self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.Samfile(line.strip(),'rb'))
                    
            fp.close()
    def __del__(self):
        for vcfname in self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE.keys():
            for pysamobj in self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname][1:]:
                pysamobj.close()
    def insertorUpdatetopleveltable(self,PopSnpAligned,flankseqfafile,SNP_flanklen):
        """
            first element of a snp is the snp position
            PopSnpAligned are required to be ordered
        """
        for chrom in PopSnpAligned.keys():
            #pre process for context
            RefSeqMap=Util.getRefSeqBypos_faster(self.refseqfahandler, self.refseqfafasteridx, chrom, PopSnpAligned[chrom][0][0]-SNP_flanklen, PopSnpAligned[chrom][-1][0]+SNP_flanklen, self.currentchrLen)
            #pre process for context end
            SNPrec_of_one_chrom_invcf=[]
            insertsql_statement_list=[];updatesql_statement_list=[]
            insertsql_data_list=[];updatesql_date_list=[]
            for snp in PopSnpAligned[chrom]:
                snp_pos=snp[0]
                snp_recINtopleveltable=self.dbvariantstools.operateDB("select","select * from "+self.toplevelsnptablename+" where chrID='"+chrom+"' and snp_pos="+str(snp_pos))
                if not snp_recINtopleveltable or snp_recINtopleveltable==0:#empty
                    insertsql_statement="insert into "+self.toplevelsnptablename + " (chrID,snp_pos,snpID,ref_base,alt_base,context,"
                    insertsql_date=[chrom,snp_pos,".",snp[1],snp[2]]
                    updatesql_statement=None
                else:
                    """
                    snp_recINtopleveltable [('KB743962.1', 60011, '.', 'A', 'T', 'TAA', None, None, 'T', '39,0', None, None, 'T', '21,0')]

                    """
                    print("snp_recINtopleveltable",snp_recINtopleveltable)
                    print(snp)
                    insertsql_statement=None
                    updatesql_date=[]
                    updatevcflist=[]
                    alt_idx=6
                    for outgroupname_ALT in snp_recINtopleveltable[0][6::2]:
                        if None==outgroupname_ALT:
                            updatevcflist.append(self.toplevelsnptable_titlelist[alt_idx][:-4])#know which field need to be filled
                        alt_idx+=2
                        
                    updatesql_statement="update "+ self.toplevelsnptablename+" set "
                for vcfname in self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE.keys():
                    archicpop_colname=re.search(r'[^/]*$',vcfname).group(0)
                    archicpop_colname=re.sub(r"[^\w^\d]","_",archicpop_colname)
                    if insertsql_statement==None and archicpop_colname in updatevcflist:
                        updatesql_statement+=(archicpop_colname+"_alt =%s,"+archicpop_colname+"_dep =%s,")
                    elif updatesql_statement!=None and insertsql_statement==None:
                        print("warning! "+archicpop_colname+" not in the ",updatevcflist)
                        print("updatesql_statement is ",updatesql_statement)
                        continue# this vcf is filled
                    elif updatesql_statement==None and insertsql_statement!=None:
                        insertsql_statement+=(archicpop_colname+"_alt,"+archicpop_colname+"_dep,")
                    
                    if chrom in  self.mapofSNPrecsForeachVCFpop_mapBYchrom:
                        pass
                    else:
                        self.mapofSNPrecsForeachVCFpop_mapBYchrom={}
                        self.mapofSNPrecsForeachVCFpop_mapBYchrom[chrom]={}
                        for vcfname_tmp in self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE.keys():
                            self.mapofSNPrecsForeachVCFpop_mapBYchrom[chrom][vcfname_tmp]=self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname_tmp][0].getVcfListByChrom(chrom,MQfilter=None)
                    SNPrec_of_one_chrom_invcf=self.mapofSNPrecsForeachVCFpop_mapBYchrom[chrom][vcfname]
                    """state: updatesql_statement= update toplevelsnptablename set processed_vcfname1_alt=%s,processed_vcfname1_dep=%s,
                            insertsql_statement=insert into toplevelsnptablename (chrID,snp_pos,snpID,ref_base,alt_base,context,processed_vcfname_m_alt,processed_vcfname_m_dep,
                    """
                    low=0;ALT=snp[2];high=len(SNPrec_of_one_chrom_invcf)-1
                    while low<=high:
                        mid=(low+high)>>1
                        if SNPrec_of_one_chrom_invcf[mid][0]<snp_pos:
                            low=mid+1
                        elif SNPrec_of_one_chrom_invcf[mid][0]>snp_pos:
                            high=mid-1
                        else:#find the pos
                            print("here is a  bug still not finished",mid)
                            pos, REF, archicpop_ALT, INFO,FORMAT,samples = SNPrec_of_one_chrom_invcf[mid]
                            dp4=re.search(r"DP4=(\d*)", INFO)
                            AF=re.search(r"AF=([\d\.e-]+)", INFO).group(1)
                            refdep=0;altalleledep=0
                            if dp4!=None and AF!=None: #samtools
                                refdep = int(dp4.group(1))*(1- float(AF))
                                altalleledep = int(dp4.group(1)) * float(AF)
                            else:#vcf from indvd 
                                if "AD" not in FORMAT:
                                    print("AD not in ",FORMAT )
                                    continue
                                AD_idx=(re.split(":",FORMAT)).index("AD")#gatk GT:AD:DP:GQ:PL 
                                for sample in samples:
                                    if len(re.split(r":",sample))==1:# ./.
                                        continue
                                    AD_depth=re.split(r",",re.split(":",sample)[AD_idx])
                                    try:
                                        refdep+=int(AD_depth[0])
                                        altalleledep+=int(AD_depth[1])
                                    except ValueError:
                                        print("Ancestralallele.fillAncestral except ValueError",sample,end="")
                            popsdata_alt=archicpop_ALT
                            popsdata_dep=str(refdep)+","+str(altalleledep)
                            break
                    else:
                        sum_depth=0
                        popsdata_alt=ALT
                        for samfile in self.outgroupVcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname][1:]:
                            ACGTdep=samfile.count_coverage(chrom,snp_pos-1,snp_pos)
                            for dep in ACGTdep:
                                sum_depth+=dep[0]
                        if sum_depth<4:
                            popsdata_dep="no covered"
                        else:
                            popsdata_dep=str(sum_depth) + ",0"
                    #continnue sql_statement
                    """state: updatesql_date=[],insertsql_date=[chrom,snp_pos,".",snp[1],snp[2]]
                    updatesql_statement= update toplevelsnptablename set processed_vcfname1_alt=%s,processed_vcfname1_dep=%s,
                    insertsql_statement=insert into toplevelsnptablename (chrID,snp_pos,snpID,ref_base,alt_base,context,processed_vcfname_m_alt,processed_vcfname_m_dep,
                    """
                    if updatesql_statement!=None and insertsql_statement==None:
                        updatesql_date+=[popsdata_alt,popsdata_dep]
                        print("updatesql_date",updatesql_date)
                    elif  updatesql_statement==None and insertsql_statement!=None:
                        insertsql_date+=[popsdata_alt,popsdata_dep]
                    else:
                        print("error: updatesql_statement!=None or insertsql_statement!=None")
                #for one snp
                #process context 
                if (updatesql_statement!=None and  snp_recINtopleveltable[0][5]==None) or insertsql_statement!=None:
                    snpID=chrom+"_"+str(snp[0])
                    if snp_pos + SNP_flanklen<=RefSeqMap[chrom][0]+len(RefSeqMap[chrom])-1 and snp_pos- SNP_flanklen>RefSeqMap[chrom][0]:
                        snpflankseq="".join(RefSeqMap[chrom][(snp_pos-SNP_flanklen-RefSeqMap[chrom][0]):(snp_pos+SNP_flanklen-RefSeqMap[chrom][0]+1)])
                        if insertsql_statement!=None:
                            insertsql_date.insert(5, snpflankseq[SNP_flanklen-1:SNP_flanklen+2])                        
                        elif snp_recINtopleveltable[0][5]==None:
                            updatesql_statement+="context=%s"
                            updatesql_date.append(snpflankseq[SNP_flanklen-1:SNP_flanklen+2])
                            print("updatesql_date",updatesql_date)
                        currentsnpID=chrom+"_"+str(snp_pos)+snpflankseq[SNP_flanklen]+":"+snp[1]
                        snpflankseq=snpflankseq[0:SNP_flanklen]+'N'+snpflankseq[SNP_flanklen+1:]

                    elif snp_pos <=RefSeqMap[chrom][0]+len(RefSeqMap[chrom])-1 and snp_pos-SNP_flanklen>RefSeqMap[chrom][0]:
                        snpflankseq="".join(RefSeqMap[chrom][(snp_pos-SNP_flanklen-RefSeqMap[chrom][0]):(snp_pos-RefSeqMap[chrom][0]+1)])
                        if insertsql_statement!=None:
                            insertsql_date.insert(5,snpflankseq[SNP_flanklen-1:SNP_flanklen+1]+"N")                        
                        elif snp_recINtopleveltable[0][5]==None:
                            updatesql_statement+="context=%s"
                            updatesql_date.append(snpflankseq[SNP_flanklen-1:SNP_flanklen+1]+"N")
                            print("updatesql_date",updatesql_date)
                        currentsnpID=chrom+"_"+str(snp_pos)+snpflankseq[SNP_flanklen]+":"+snp[1]
                        snpflankseq=snpflankseq[0:SNP_flanklen]+'N'
                        
                    elif snp_pos-SNP_flanklen<=RefSeqMap[chrom][0] and snp_pos + SNP_flanklen<=RefSeqMap[chrom][0]+len(RefSeqMap[chrom])-1:
                        snpflankseq="".join(RefSeqMap[chrom][(snp_pos - RefSeqMap[chrom][0]):(snp_pos + SNP_flanklen - RefSeqMap[chrom][0] + 1)])
                        if insertsql_statement!=None:
                            insertsql_date.insert(5,"N"+snpflankseq[0:2])                        
                        elif snp_recINtopleveltable[0][5]==None:
                            updatesql_statement+="context=%s"
                            updatesql_date.append("N"+snpflankseq[0:2])
                            print("updatesql_date",updatesql_date)
                        currentsnpID=chrom+"_"+str(snp_pos)+snpflankseq[0]+":"+snp[1]
                        snpflankseq = 'N'+snpflankseq[1:SNP_flanklen+1]
                        
                    else:
                        print(len(RefSeqMap[chrom]))
                        currentsnpID=">error"
                        snpflankseq="error"
                        print("Error! what's wrong with the func insertorUpdatetopleveltable?")
                    print(">" + currentsnpID + "\n" + snpflankseq, end='\n', file=flankseqfafile)
                #context process end ,append into multiple statements    
                if updatesql_statement!=None:
                    if snp_recINtopleveltable[0][5]!=None:
                        updatesql_statement=updatesql_statement[:-1]
                    updatesql_statement+=" where snp_pos="+str(snp_pos)+" and chrID='"+chrom+"'"
                    updatesql_statement_list.append(updatesql_statement)
                    print("updatesql_date_list",updatesql_date_list)
                    updatesql_date_list.append(tuple(updatesql_date))
                    print("updatesql_date_list",updatesql_date_list)
                elif insertsql_statement!=None:
                    insertsql_statement=insertsql_statement[:-1]+")VALUES ("+"%s,"*(len(insertsql_date))
                    insertsql_statement=insertsql_statement[:-1]+")"
                    insertsql_statement_list.append(insertsql_statement)
                    insertsql_data_list.append(tuple(insertsql_date))
#                     self.dbvariantstools.operateDB("",insertsql_statement,data=tuple(insertsql_date))
            if insertsql_statement_list !=[]:
#                 print(insertsql_statement_list)
                self.dbvariantstools.operateDB("insert",*insertsql_statement_list,data=insertsql_data_list)
#                 snp=self.dbvariantstools.operateDB("select","select * from "+self.toplevelsnptablename+" where chrID='"+insertsql_data_list[0][0]+"' and snp_pos='"+str(insertsql_data_list[0][1])+"'")
#                 print(snp)
            elif updatesql_statement_list !=[] and len(updatesql_date_list)==len(updatesql_statement_list):
                print(updatesql_statement_list,updatesql_date_list)
                self.dbvariantstools.operateDB("update",*updatesql_statement_list,data=updatesql_date_list)
                
#             print("insertsql_data_list",insertsql_data_list)
#             print("updatesql_date_list",updatesql_date_list)
#             snp=self.dbvariantstools.operateDB("select","select * from "+self.topleveltablejudgeancestralname+" where chrID='"+insertsql_data_list[0][0]+"' and snp_pos='"+str(insertsql_data_list[0][1])+"'")
               
                