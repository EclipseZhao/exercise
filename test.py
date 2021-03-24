#!/usr/bin/python
#coding: utf-8

import os
import re
import sys
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description = 'annotate the inheritance mode of trio family', formatter_class = RawTextHelpFormatter)
##parser.add_argument('-M', '--model', metavar = 'String', required = True, help = 'Dominant(D)|Recessive(R)|Compound heterozygous(C).')
parser.add_argument('-p', '--ped', metavar = 'File', required = True, help = '#FamilyID, SampleID and parient must be inclued.')
parser.add_argument('-i', '--in', metavar = 'File', required = True, help = 'The merged file use to do inheritance analysis with .xls(.gz) format.')
#parser.add_argument('-', '--parent',metavar = 'Float', required = True, help = 'the frequency to filter the variation site, default=0.05',default=0.05)
parser.add_argument('--denovo',metavar = 'String', required = True, help = 'the denovo file, if there are more than one denovo vcf file, then split wiht ",", for example denovo1.vcf,denovo2.vcf,denovo3.vcf....')
parser.add_argument('-t', '--type', metavar = 'String', required = True, help = '', choices=['trio3', 'other'], default='trio3')
parser.add_argument('-o', '--out', metavar = 'File', default  = './', help = 'The output file ')
argv = vars(parser.parse_args())

inf = argv['in'].strip()
pedf = argv['ped'].strip()
outf = argv['out'].strip()
denovof=argv['denovo'].strip()
##f = argv['freq'].strip()
t=argv['type'].strip() 

def getpos(name,title,name2=''):
    if name in title:
        pos = title.index(name)
    elif name2 != '' and name2 in title:
        pos = title.index(name2)
    else:
        if name2 != '':
           exit('Warning: %s and %s not in title.' %(name,name2))
        else:
           exit('Warning: %s not in title.' %name)
    return pos


##return list gt=[type1,type2,genotype,DV,DP,vr,lowGQ,gt]
##the item "lowGQ" in gt list only exists when the format has the lowGQ tag
def get_gt(line,fm_list,sample):
    ad_pos=getpos('AD',fm_list)
    dp_pos=getpos('DP',fm_list)
    gt_list=sample.split(':')
    if gt_list[0] != "./.":
        if re.match(r'\d+',sample):
            sample_format = re.search('(\d+)[\/|\|](\d+)',sample)
            gt=[sample_format.group(1),sample_format.group(2)]
            if sample_format.group(1) == sample_format.group(2) and sample_format.group(1) != '0':
                gt.append('hom')
            elif sample_format.group(1) == sample_format.group(2) and sample_format.group(1) == '0':
                gt.append('ref')
            elif sample_format.group(1) != sample_format.group(2) and \
                (sample_format.group(1) == '0' or sample_format.group(2) == '0'):
                gt.append('het')
            else:
                gt.append('chet')  ##compound heterozygous
            if 'ref' not in gt:
                if int(gt[0]) <2 and int(gt[1]) < 2:
                    dv=int(gt_list[ad_pos].split(',')[1])
                    dp=int(gt_list[dp_pos])
                    if dv==0 or dp==0:
                        vr=0
                        print 'check DV:\t'+line
                    else:
                        vr=float(dv)/float(dp)
                else:
                    v_dep=gt_list[ad_pos].split(',')
                    dp=int(gt_list[dp_pos])
                    dv=0
                    for i in range(len(v_dep)):
                        if i ==0:
                            continue
                        else:
                            dv+=int(v_dep[i])
                    if dv==0 or dp ==0:
                        vr=0
                        print 'check DV:\t'+line
                    else:
                        vr=float(dv)/float(dp)
            else:
                dv=0
                dp=gt_list[2]
                vr=0
            gt.append(dv)
            gt.append(dp)
            gt.append(vr)
            if 'lowGQ' in gt_list:
                gt.append('lowGQ')
            gt.append(gt_list[0])
            return gt
    else:
        return "missing genotype"

denovo_d={}

if denovof != 'no':
    if len(denovof.split(',')) ==1:
        with open(denovof,'r') as denovof:
            for line in denovof:
                if line.startswith('#'):
                    continue
                else:
                    if 'hiConfDeNovo' in line:
                        dat_l=line.strip().split('\t')
                        k=dat_l[0]+'_'+dat_l[1]+'_'+dat_l[3]+'_'+dat_l[4]
                        denovo_d[k]=line
    else:
        denovol=denovof.split(',')
        for i in range(len(denovol)):
            with open(denovol[i],'r') as tmpf:
                for line in tmpf:
                    if line.startswith('#'):
                        continue
                    else:
                        if 'hiConfDeNovo' in line:
                            dat_l=line.strip().split('\t')
                            k=dat_l[0]+'_'+dat_l[1]+'_'+dat_l[3]+'_'+dat_l[4]
                            denovo_d[k]=line



if t=="trio3":
    with open(pedf,'r') as pedf:
        patient_l=[]
        normal_l=[]
        p_sum=0
        n_sum=0
        for line in pedf:
            d_ped=line.strip().split('\t')
            fid=d_ped[0]
            if d_ped[2] != '0' and d_ped[3] != '0':
                kid=d_ped[1]
                fa=d_ped[2]
                ma=d_ped[3]
                if d_ped[4]=='1':
                    kid_gender='M'
                elif d_ped[4]=='2':
                    kid_gender="F"
                else:
                    print "The gender in ped file is mistaken, please check"
            if d_ped[5] =="2":
                patient_l.append(d_ped[1])
                p_sum+=1
            else:
                normal_l.append(d_ped[1])
                n_sum+=1
        dnkid='hiConfDeNovo='+kid
        print dnkid
    ##outf= outd+'/Trio_'+fid+"_inheritance_mode.xls"
    #outf= outd+'/Trio_'+fa+"_"+ma+"_"+kid+"_inheritance_mode.xls"
    print kid
    with open(inf,'r') as inf, open(outf,'w') as outf:
        p_pos=[]
        n_pos=[]
        for line in inf:
            dat_l=line.strip().split('\t')
            p_n=0
            n_n=0
            ##mod_l list will have 3 items:1 it's a patient's site or not;2 inherit from;3 Hereditary type 
            if 'Ref' in dat_l and 'Alt' in dat_l:
                outf.write('\t'.join(dat_l))
                outf.write('\tDenovo\tCompound_heterozygosity\tInheritance_Mode\n')
                kid_pos=getpos(kid,dat_l)
                ma_pos=getpos(ma,dat_l)
                fa_pos=getpos(fa,dat_l)
                chr_pos=getpos('Chr',dat_l)
                start_pos=getpos('Start',dat_l)
                ref_pos=getpos('Ref',dat_l)
                alt_pos=getpos('Alt',dat_l)
                ##gene_pos=getpos('Gene.refGeneWithVer',dat_l)
                info_pos=getpos('INFO',dat_l)
                fm_pos=getpos('FORMAT',dat_l)
                for i in patient_l:
                    p_pos.append(getpos(i,dat_l))
                for j in normal_l:
                    n_pos.append(getpos(j,dat_l))
            else:
                fm_list=dat_l[fm_pos].split(':')
                kid_gt=get_gt(line,fm_list,dat_l[kid_pos])
                fa_gt=get_gt(line,fm_list,dat_l[fa_pos])
                ma_gt=get_gt(line,fm_list,dat_l[ma_pos])
                info_list=dat_l[info_pos].split(':')
                ##if the site belong to the patient's site and normal sample don't have the same varitation
                for i in p_pos:
                    if 'ref' not in get_gt(line,fm_list,dat_l[i]) and "missing genotype" not in get_gt(line,fm_list,dat_l[i]):
                        p_n+=1
                for j in n_pos:
                    if 'ref' not in get_gt(line,fm_list,dat_l[j]) and "missing genotype" not in get_gt(line,fm_list,dat_l[j]):
                        n_n+=1
                ##determine the hereditary mode
                ke=dat_l[chr_pos]+'_'+dat_l[start_pos]+'_'+dat_l[ref_pos]+'_'+dat_l[alt_pos]
                if ke in denovo_d:
                    dat_l.append('Y')
                else:
                    dat_l.append('N')
                if dat_l[chr_pos] != 'ChrX' and dat_l[chr_pos] !="ChrY":
                    if kid_gt[2]=='chet' or fa_gt[2]=='chet' or ma_gt[2]=='chet':
                        dat_l.append('Y')
                    else:
                        dat_l.append('N')
                    if kid_gt[2]=='hom':
                        if kid_gt[0] in fa_gt and kid_gt[0] in ma_gt:
                            if 'hom' not in fa_gt and 'hom' not in ma_gt:
                                if fa not in patient_l and ma not in patient_l and kid in patient_l:
                                    dat_l.append("Recessive")
                                else:
                                    dat_l.append("UNKNOWN")
                            elif 'hom' in fa_gt and 'hom' not in ma_gt:
                                if fa in patient_l and kid in patient_l:
                                    dat_l.append("Recessive")
                                else:
                                    dat_l.append("UNKNOWN")
                            elif 'hom' not in fa_gt and 'hom' in ma_gt:
                                if ma in patient_l and kid in patient_l:
                                    dat_l.append("Recessive")
                                else:
                                    dat_l.append("UNKNOWN")
                            elif 'hom' in fa_gt and 'hom' in ma_gt:
                                if fa in patient_l and kid in patient_l and ma in patient_l:
                                    dat_l.append("Recessive")
                                else:
                                    dat_l.append("UNKNOWN")
                            else:
                                dat_l.append("UNKNOWN")
                        else:
                            dat_l.append("UNKNOWN")
                        #elif kid_gt[0] in fa_gt and kid_gt[0] not in ma_gt:
                        #    mod_l.append("Paternal")
                        #    mod_l.append(".")
                        #elif kid_gt[0] not in fa_gt and kid_gt[1] in ma_gt:
                        #    mod_l.append("Maternal")
                        #    mod_l.append(".")
                        #else:
                        #    mod_l.append(".")
                        #    mod_l.append("UNKNOWN")
                        #    print 'unknown:\t'+ dat_l[kid_pos]+'\t'+dat_l[fa_pos]+'\t'+dat_l[ma_pos]+'\t'+line
                        #mod_line=";".join(mod_l[1:])
                    elif kid_gt[2]=='het':
                        if kid_gt[1] in fa_gt and kid_gt[1] in ma_gt:
                            if fa in patient_l and kid in patient_l and ma in patient_l:
                                dat_l.append("Dominant")
                            else:
                                dat_l.append("UNKNOWN")
                        elif kid_gt[1] in fa_gt and kid_gt[1] not in ma_gt:
                            if fa in patient_l and kid in patient_l:
                                dat_l.append("Dominant")
                            else:
                                dat_l.append("UNKNOWN")
                        elif kid_gt[1] in ma_gt and kid_gt[1] not in fa_gt:
                            if ma in patient_l and kid in patient_l:
                                dat_l.append("Dominant")
                            else:
                                dat_l.append("UNKNOWN")
                        elif kid_gt[1] in fa_gt and kid_gt[1] not in ma_gt:
                            if ma in patient_l and kid in patient_l:
                                dat_l.append("Dominant")
                            else:
                                dat_l.append("UNKNOWN")
                        else:
                            dat_l.append("UNKNOWN")
                    elif kid_gt[2]=='chet':
                        if (kid_gt[0] in fa_gt and kid_gt[1] in ma_gt) or (kid_gt[0] in ma_gt and kid_gt[1] in fa_gt):
                            dat_l.append("Compound heterozygosity")
                        else:
                            dat_l.append("UNKNOWN")
                    else:
                        dat_l.append("UNKNOWN")
                elif dat_l[chr_pos] == "ChrX":
                    if kid_gender=="M":
                        if 'chet' in ma_gt:
                            dat_l.append('Y')
                        else:
                            dat_l.append('N')
                        if kid_gt[0] in ma_gt:
                            if kid in patient_l and ((ma_gt[2]=='het' and ma not in patient_l) or ((ma_gt[2] == 'hom' or ma_gt[2]=='chet') and ma in patient_l)) and pa not in patient_l:
                                dat_l.append('Recessive_X')
                            else:
                                dat_l.append("UNKNOWN")
                        else:
                            dat_l.append("UNKNOWN")
                    elif kid_gender=="F":
                        if kid_gt[2]=='hom':
                            if kid_gt[0] in fa_gt and kid_gt[0] in ma_gt:
                                if kid in patient_l and fa in patient_l and ma in patient_l and ('hom' in ma_gt or 'chet' in ma_gt):
                                    dat_l.append("Recessive_X")
                                elif kid in patient_l and fa in patient_l and ma not in Paternal and 'hom' not in ma_gt and 'chet'  not in ma_gt:
                                    dat_l.append("Recessive_X")
                                else:
                                    dat_l.append("UNKNOWN")
                            else:
                                dat_l.append("UNKNOWN")
                        elif kid_gt[2]=='het':
                            if kid_gt[1] in fa_gt and kid_gt[1] in ma_gt:
                                if kid in patient_l and fa in patient_l and ma in patient_l:
                                    dat_l.append("Dominant_X")
                                else:
                                    dat_l.append('UNKNOWN')
                            elif kid_gt[1] in fa_gt and kid_gt[1] not in ma_gt:
                                if kid in patient_l and fa in patient_l and ma not in patient_l:
                                    mod_l.append("Dominant_X")
                                else:
                                    dat_l.append('UNKNOWN')
                            elif kid_gt[1] not in fa_gt and kid_gt[1] in ma_gt:
                                if kid in patient_l and fa not in patient_l and ma in patient_l:
                                    dat_l.append("Dominant_X")
                                else:
                                    dat_l.append('UNKNOWN')
                            else:
                                dat_l.append("UNKNOWN")
                        elif kid_gt[2]=='chet':
                            if (kid_gt[0] in fa_gt and kid_gt[1] in ma_gt) or (kid_gt[0] in ma_gt and kid_gt[1] in fa_gt):
                                dat_l.append("Compound heterozygosity")
                            else:
                                dat_l.append("UNKNOWN")
                        else:
                            dat_l.append("UNKNOWN")
                    else:
                        print "Genda mistake"
                elif dat_l[chr_pos] =="ChrY":
                    dat_l.append('N')
                    if kid_gender == "M":
                        if 'hom' in kid_gt and  kid_gt[0] in fa_gt:
                            if kid in patient_l and fa in patient_l and ma not in patient_l:
                                dat_l.append('Y_link')
                            else:
                                mod_l.append("UNKNOWN")
                        else:
                            dat_l.append("UNKNOWN")
                    elif kid_gender=="F":
                        dat_l.append("UNKNOWN")
                else:
                    dat_l.append("UNKNOWN")
                if dat_l[-1] != 'UNKNOWN' or dat_l[-2] != 'N' or dat_l[-3] !='N':
                    outf.write('\t'.join(dat_l))
                    outf.write('\n')
elif t=="other":
    patient_d={}
    with open(pedf,'r') as pedf:
        patient_l=[]
        normal_l=[]
        gender_p={}
        gender_n={}
        p_sum=0
        n_sum=0
        for line in pedf:
            d_ped=line.strip().split('\t')
            fid=d_ped[0]
            if d_ped[5] =="2":
                patient_l.append(d_ped[1])
                p_sum+=1
                if d_ped[4] == '1':
                    gender_p[d_ped[1]]='M'
                elif d_ped[4] == '2':
                    gender_p[d_ped[1]]='F'
            else:
                normal_l.append(d_ped[1])
                n_sum+=1
                if d_ped[4] == '1':
                    gender_n[d_ped[1]]='M'
                elif d_ped[4] == '2':
                    gender_n[d_ped[1]]='F'
    ##outf= outd+'/Trio_'+fid+"_inheritance_mode.xls"
    #outf= outd+'/Trio_'+fa+"_"+ma+"_"+kid+"_inheritance_mode.xls"
    with open(inf,'r') as inf, open(outf,'w') as outf:
        p_pos=[]
        n_pos=[]
        for line in inf:
            dat_l=line.strip().split('\t')
            mod_l=[]
            p_n=0
            n_n=0
            ##mod_l list will have 3 items:1 it's a patient's site or not;2 inherit from;3 Hereditary type 
            if 'Ref' in dat_l and 'Alt' in dat_l:
                outf.write('\t'.join(dat_l))
                outf.write('\tDenovo\tCompound_heterozygosity\tInheritance_Mode\n')
                #kid_pos=getpos(kid,dat_l)
                #ma_pos=getpos(ma,dat_l)
                #fa_pos=getpos(fa,dat_l)
                chr_pos=getpos('Chr',dat_l)
                start_pos=getpos('Start',dat_l)
                ref_pos=getpos('Ref',dat_l)
                alt_pos=getpos('Alt',dat_l)
                ##gene_pos=getpos('Gene.refGeneWithVer',dat_l)
                info_pos=getpos('INFO',dat_l)
                fm_pos=getpos('FORMAT',dat_l)
                for i in patient_l:
                    p_pos.append(getpos(i,dat_l))
                for j in normal_l:
                    n_pos.append(getpos(j,dat_l))
            else:
                patient_d={}
                normal_d={}
                fm_list=dat_l[fm_pos].split(':')
                for i in range(len(patient_l)):
                    patient_d[patient_l[i]]=get_gt(line,fm_list,dat_l[p_pos[i]])
                for i in range(len(normal_l)):
                    normal_d[normal_l[i]]=get_gt(line,fm_list,dat_l[n_pos[i]])
                #kid_gt=get_gt(line,fm_list,dat_l[kid_pos])
                #fa_gt=get_gt(line,fm_list,dat_l[fa_pos])
                #ma_gt=get_gt(line,fm_list,dat_l[ma_pos])
                info_list=dat_l[info_pos].split(':')
                ##if the site belong to the patient's site and normal sample don't have the same varitation
                ke=dat_l[chr_pos]+'_'+dat_l[start_pos]+'_'+dat_l[ref_pos]+'_'+dat_l[alt_pos]
                if denovof == 'no':
                    dat_l.append('N')
                elif ke in denovo_d:
                    dat_l.append('Y')
                else:
                    dat_l.append('N')
                if dat_l[chr_pos] !="chrX" and dat_l[chr_pos]!="chrY":
                    het_p=0
                    hom_p=0
                    chet_p=0
                    no_p=0
                    no_n=0
                    het_n=0
                    chet_n=0
                    hom_n=0
                    for i in range(len(patient_l)):
                        if patient_d[patient_l[i]][2]=="het":
                            het_p+=1
                        elif patient_d[patient_l[i]][2]=="chet":
                            chet_p+=1
                        elif patient_d[patient_l[i]][2]=="hom":
                            hom_p+=1
                        else:
                            no_p+=1
                    for i in range(len(normal_l)):
                        if normal_d[normal_l[i]][2]=='het':
                            het_n+=1
                        elif normal_d[normal_l[i]][2]=='chet':
                            chet_n+=1
                        elif normal_d[normal_l[i]][2]=='hom':
                            hom_n+=1
                        else:
                            no_n+=1
                    if chet_p !=0 and chet_n!=0:
                        dat_l.append('Y')
                    elif chet_p !=0 and chet_n==0:
                        dat_l.append('patient_Y')
                    elif chet_n != 0 and chet_p==0:
                        dat_l.append('normal_Y')
                    else:
                        dat_l.append('N')
                    if (het_p+chet_p+hom_p)==p_sum and no_n == n_sum:
                        dat_l.append('Dominant')
                    elif (chet_p+hom_p)==p_sum and (het_n+no_n)==n_sum:
                        dat_l.append('Recessive')
                    else:
                        dat_l.append('UNKNOWN')
                elif dat_l[chr_pos] =="chrX":
                    fhet_p=0
                    mhom_p=0
                    fhom_p=0
                    fchet_p=0
                    no_p=0
                    no_n=0
                    fhet_n=0
                    fchet_n=0
                    mhom_n=0
                    fhom_n=0
                    for i in range(len(patient_l)):
                        if gender_p[patient_l[i]] == 'M' and patient_d[patient_l[i]][2]=="hom":
                            mhom_p+=1
                        elif  gender_p[patient_l[i]] == 'F' and (patient_d[patient_l[i]][2]=="hom" or patient_d[patient_l[i]][2]=="chet"):
                            fhom_p+=1
                            fchet_p+=1
                        elif gender_p[patient_l[i]] == 'F' and patient_d[patient_l[i]][2]=="het":
                            fhet_p+=1
                        else:
                            no_p+=1
                    for i in range(len(normal_l)):
                        if gender_n[normal_l[i]] == 'M' and normal_d[normal_l[i]][2]=="hom":
                            mhom_n+=1
                        elif  gender_n[normal_l[i]] == 'F' and (normal_d[normal_l[i]][2]=="hom" or normal_d[normal_l[i]][2]=="chet"):
                            fhom_n+=1
                            fchet_n+=1
                        elif gender_n[normal_l[i]] == 'F' and normal_d[normal_l[i]][2]=="het":
                            fhet_n+=1
                        else:
                            no_n+=1
                    if fchet_p !=0 and fchet_n!=0:
                        dat_l.append('Y')
                    elif fchet_p !=0 and fchet_n==0:
                        dat_l.append('patient_Y')
                    elif fchet_n != 0 and fchet_p==0:
                        dat_l.append('normal_Y')
                    else:
                        dat_l.append('N')
                    if (mhom_p+fhom_p)==p_sum and (no_n+fhet_n)==n_sum:
                        dat_l.append('Recessive_X')
                    elif (mhom_p+fhom_p+fchet_p+fhet_p)==p_sum and no_n==n_sum:
                        dat_l.append('Dominant_X')
                    else:
                        dat_l.append('UNKNOWN')
                elif dat_l[chr_pos] =="chrY":
                    p=0
                    mn=0
                    fn=0
                    dat_l.append('N')
                    for i in range(len(patient_l)):
                        if gender_p[patient_l[i]] == 'M' and patient_d[patient_l[i]][2]=="hom":
                            p+=1
                    for i in range(len(normal_l)):
                        if gender_n[normal_l[i]] == 'M' and normal_d[normal_l[i]][2] != "hom":
                            mn+=1
                        elif gender_n[normal_l[i]] == 'F':
                            fn+=1
                    if p==p_sum and mn == 0:
                        dat_l.append('Y_link')
                    else:
                        dat_l.append('UNKNOWN')
                if dat_l[-1] != 'UNKNOWN' or dat_l[-2] != 'N' or dat_l[-3] !='N':
                    outf.write('\t'.join(dat_l))
                    outf.write('\n')

