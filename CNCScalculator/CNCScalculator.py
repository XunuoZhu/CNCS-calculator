#!/usr/bin/python
# -*- coding: UTF-8 -*-

from re import match,findall,compile
from math import log,e,inf
from decimal import Decimal as dec
from scipy.stats import chi2,chi2_contingency
from scipy.special import gamma as ga
from progressbar import ProgressBar,ETA,Bar,Percentage
from sys import exit
from pkg_resources import resource_filename as data_path
import pandas as pd
from numpy import array
import copy

def open_file(file:str,mode:str = 'r'):
    '''
    file:dir of file to open\n
    mode:file open mode\n
    open file with warning modification.
    '''
    try:
        RET = open(file,mode)
        return RET
    except:
        exit("Can't open file: %s ."%file)
    
def count_mutation(mut_file:str,Nonsynonymous_site_with_mutation_initial:dict,Synonymous_site_with_mutation_initial:dict,Mut_missense_initial:dict,Mut_silent_initial:dict):
    '''
    mut_file:dir of mutation file\n
    Nonsynonymous_site_with_mutation_initial:dict of initial nonsynonymous site with mutation, key is enst, value is zero\n
    Mut_missense_initial:dict of initial missense mutation, key is enst, value is zero\n
    Mut_silent_initial:dict of initial silent mutation, key is enst, value is zero\n
    Count the number of non-synonymous, synonymous and nonsense mutations for each gene
    '''
    MUT = pd.read_csv(mut_file,sep='\t',header=0)
    Aa_mut_count = {}
    Nucl_mis_mut_count = {}
    Nucl_syn_mut_count = {}
    Nonsynonymous_site_with_mutation = copy.deepcopy(Nonsynonymous_site_with_mutation_initial)
    Synonymous_site_with_mutation = copy.deepcopy(Synonymous_site_with_mutation_initial)
    Mut_missense = copy.deepcopy(Mut_missense_initial)
    Mut_silent = copy.deepcopy(Mut_silent_initial)
    # use package 'pandas to realize the function of count_mutation'
    if 'Hugo_Symbol' in MUT.columns:
        mut_type_pos = MUT.columns.get_loc('Consequence')
        var_type_pos = MUT.columns.get_loc('Variant_Type')
        aa_mut_pos = MUT.columns.get_loc('HGVSp_Short')
        ensg_pos = MUT.columns.get_loc('Gene')
        enst_pos = MUT.columns.get_loc('Feature')
        nucl_mut_pos = MUT.columns.get_loc('HGVSc')
        for index in MUT.index:
            mutation_type = MUT.iloc[index,mut_type_pos]
            var_type = MUT.iloc[index,var_type_pos]
            if (mutation_type == 'missense_variant') and (var_type == 'SNP'):
                enst = MUT.iloc[index,enst_pos]
                try:
                    Mut_missense[enst] += 1
                except KeyError:
                    continue
                pro_mut = MUT.iloc[index,aa_mut_pos].split('.')[1]
                id = '\t'.join([MUT.iloc[index,ensg_pos],enst,pro_mut])
                Aa_mut_count[id] = Aa_mut_count.get(id,0) + 1
                nucl_mut_type = MUT.iloc[index,nucl_mut_pos]
                if findall(r'c.(\d+)[A-Z]',nucl_mut_type):
                    nucl_pos = findall(r'c.(\d+)[A-Z]',nucl_mut_type)[0]
                nucl_mut = enst + ':' + nucl_pos
                if nucl_mut in Nucl_mis_mut_count.keys():
                    Nucl_mis_mut_count[nucl_mut] += 1
                else:
                    Nucl_mis_mut_count[nucl_mut] = 1
                    Nonsynonymous_site_with_mutation[enst] += 1
            elif mutation_type == 'synonymous_variant':
                enst = MUT.iloc[index,enst_pos]
                try:
                    Mut_silent[enst] += 1
                except KeyError:
                    continue
                nucl_mut_type = MUT.iloc[index,nucl_mut_pos]
                if findall(r'c.(\d+)[A-Z]',nucl_mut_type):
                    nucl_pos = findall(r'c.(\d+)[A-Z]',nucl_mut_type)[0]
                nucl_mut = enst + ':' + nucl_pos
                if nucl_mut in Nucl_syn_mut_count.keys():
                    Nucl_syn_mut_count[nucl_mut] += 1
                else:
                    Nucl_syn_mut_count[nucl_mut] = 1
                    Synonymous_site_with_mutation[enst] += 1

    del MUT
    print("Count the number of somatic mutations for each gene! File name: %s\n"%mut_file)
    return Aa_mut_count,Nonsynonymous_site_with_mutation,Synonymous_site_with_mutation,Mut_missense,Mut_silent

def prepare_fas_cds_file(cds_file = None ,fas_file = None):
    '''
    cds_file:dir of cds file, None means use default file\n
    fas_file:dir of fasta file, None means use default file\n
    Read cds file and fasta file, return dicts of gene features, synonymous sites, nonsynonymous sites and so on.
    '''
    ## Load cds file.
    if cds_file is None:
        CDS = open_file(data_path('CNCScalculator.CNCScalculator','data/TCGA_match_GRCh37.cds'))
    else:
        CDS = open_file(cds_file)

    head_cds = CDS.readline()
    Aa_length = {}
    Gene_id = {}
    Synonymous ={}
    Non_synonymous = {}
    Mut_missense_initial = {} ; Mut_silent_initial = {}
    Nonsynonymous_site_with_mutation_initial = {}
    Synonymous_site_with_mutation_initial = {}
    for line in CDS.readlines():
        line = line.strip('\n').strip('\r')
        cds = line.split('\t')
        Aa_length[cds[1]] = cds[5]
        Gene_id[cds[1]] = '\t'.join(cds[0:4])
        if cds[1] not in Synonymous.keys():
            Synonymous[cds[1]] = 0 ; Non_synonymous[cds[1]] = 0
            Mut_missense_initial[cds[1]] = 0 ; Mut_silent_initial[cds[1]] = 0
            Nonsynonymous_site_with_mutation_initial[cds[1]] = 0 ; Synonymous_site_with_mutation_initial[cds[1]] = 0
    CDS.close()
    print("Load cds file successfully!\n")

    def sub_matrix():
        A_freq,C_freq,T_freq,G_freq = {},{},{},{}
        MOD = open_file(data_path('CNCScalculator.CNCScalculator','data/pancancer_mutation_rate.txt')) 
        head = MOD.readline()
        for line in MOD.readlines():
            line = line.strip('\n').strip('\r').split('\t')
            if match(r'\wA\w',line[1]):
                A_freq[line[0]] = eval(line[2])
            elif match(r'\wT\w',line[1]):
                T_freq[line[0]] = eval(line[2])
            elif match(r'\wC\w',line[1]):
                C_freq[line[0]] = eval(line[2])
            elif match(r'\wG\w',line[1]):
                G_freq[line[0]] = eval(line[2])
        MOD.close()
        print("Load mutation rate profile!")
        return {"A":A_freq,"T":T_freq,"C":C_freq,"G":G_freq}

    def code_matrix():
        CODE = open_file(data_path('CNCScalculator.CNCScalculator',"data/genetic_code.txt"))
        head = CODE.readline()
        substitution = []
        for line in CODE.readlines():
            substitution.append(line.strip('\n').strip('\r'))
        CODE.close()
        print("Load genetic code!")
        return(substitution)
    Substitution = code_matrix()
    freq_dict = sub_matrix()

    # Load fasta file.
    print("START READING FAS FILE..")
    if fas_file is None:
        FAS = open_file(data_path('CNCScalculator.CNCScalculator',"data/human_TCGA_match_GRCh37.fas"))
    else:
        FAS = open_file(fas_file)
    
    Seq = FAS.read().split('>')[1:]
    Gene_name = {};Sequence = {}

    print("THERE ARE %s SEQ TO ANALYSE"%(len(Seq)))
    bar = ProgressBar(widgets=[Percentage(),'',Bar('#'),'',ETA(),''])
    for k in bar(range(len(Seq))):
        fas = Seq[k].split('\n')
        gene_name = fas[0]
        gene = gene_name.split('|')
        enst = gene[1]
        Sequence[enst] = ''.join(fas[1:])
        Gene_name[enst] = gene_name
        sequence = Sequence[enst]
        seq_length = len(sequence)

        ## Calculate mutation rate of each site
        for i in range(3,seq_length-3):
            fragment = sequence[i-1:i+2]
            if not fragment:exit("No fragment.")
            if 'N'in fragment: print("{}\t{}\t{}\t".format(enst,i,fragment));continue
            start = i - i%3
            codon = sequence[start:start+3]
            if not codon:exit("No codon.")
            trans_table = ''.maketrans('ATCG','TAGC')
            
            def mutation_cal(base):
                fragment_rc = fragment[::-1].translate(trans_table)
                base_rc = base.translate(trans_table)
                if fragment in freq_dict[base].keys():
                    mutation = freq_dict[base][fragment]
                elif fragment_rc in freq_dict[base_rc].keys():
                    mutation = freq_dict[base_rc][fragment_rc]
                else:
                    mutation = 0
                return mutation

            mutation_T,mutation_A,mutation_G,mutation_C = map(mutation_cal,['T','A','G','C'])
            mutation = sum([mutation_T,mutation_A,mutation_G,mutation_C])
            if mutation == 0: exit("erro:mutation == 0\n{}\t{}\t{}".format(enst,i,fragment))
            prop_T,prop_C,prop_A,prop_G = map(lambda x:x/mutation,[mutation_T,mutation_C,mutation_A,mutation_G])
                
            for sub in Substitution:
                sub_type = sub.split('\t')
                if not sub_type[0]: exit("no sub_type!")
                if sub_type[0] == codon:
                    j = i%3
                    t = 4*j+1
                    c = 4*j+2
                    a = 4*j+3
                    g = 4*j+4
                    for stype,prop in [(sub_type[t],prop_T),(sub_type[c],prop_C),(sub_type[a],prop_A),(sub_type[g],prop_G)]:
                        if stype and prop:
                            if stype == 'synonymous':Synonymous[enst] += prop
                            elif stype in ('non-synonymous','stop-loss'):Non_synonymous[enst] += prop
                    #print("gene: {}\nsynonymous: {}\nNon-synonymous: {}".format(gene[0],Synonymous[enst],Non_synonymous[enst]))
                    break

    FAS.close()
    print("Load gene sequence!")
    print("Calculate non-synonymous and synonymous sites!")
    return Gene_id,Aa_length,Gene_name,Synonymous,Non_synonymous,Nonsynonymous_site_with_mutation_initial,Synonymous_site_with_mutation_initial,Mut_missense_initial,Mut_silent_initial

def cncs(mut_file:str,outputfile:str,fas_cds:list=[],cds_file = None, fas_file = None):
    '''
    Calculate CN/CS value for each gene.\n
    mut_file: mutation filepath\n
    outputfile: output filepath\n
    fas_cds: the output of prepare_fas_cds_file function\n
    cds_file: cds file path, if fas_cds is not empty list, this parameter will be ignored\n
    fas_file: fasta file path, if fas_cds is not empty list, this parameter will be ignored\n
    Calculate the CN/CS of each gene.
    Print result to a file. And return a dic whose keys is enst id and values is CN/CS.
    '''
    if fas_cds == []:
        Gene_id,Aa_length,Gene_name,L_S,L_N,Nonsynonymous_site_with_mutation_initial,Synonymous_site_with_mutation_initial,Mut_missense_initial,Mut_silent_initial = prepare_fas_cds_file(cds_file,fas_file)
    else:
        Gene_id,Aa_length,Gene_name,L_S,L_N,Nonsynonymous_site_with_mutation_initial,Synonymous_site_with_mutation_initial,Mut_missense_initial,Mut_silent_initial = fas_cds

    # Count the number of non-synonymous, synonymous and nonsense mutations for each gene
    Aa_mut_count,S_N,S_S,Mut_N,Mut_S = count_mutation(mut_file,Nonsynonymous_site_with_mutation_initial,Synonymous_site_with_mutation_initial,Mut_missense_initial,Mut_silent_initial)

    OUT = open_file(outputfile,mode='w')
    OUT.write("Gene_id\tnonsynonymous_site\tsynonymous_site\tnonsynonymous_count\tsynonymous_count\tnonsynonymous_site_with_mutation\tsynonymous_site_with_mutation\tCN\tCS\tCN/CS\tp_value\n")
    Cncs = {}
    Cncs_p = {}
    for enst in Gene_name.keys():
        gene_name = Gene_name[enst]
        nonsynonymous_site = L_N[enst]
        synonymous_site = L_S[enst]
        nonsynonymous_count = Mut_N[enst]
        synonymous_count = Mut_S[enst]
        nonsynonymou_site_with_mutation = S_N[enst]
        synonymous_site_with_mutation = S_S[enst]
        ## chi2_contingency
        try :
            chi2_con_array = array([[nonsynonymous_site,synonymous_site],[nonsynonymous_count,synonymous_count]])
            chi2_value,p_value,df,expected = chi2_contingency(chi2_con_array)        
        except:
            p_value = 'NA'

        ## cncs calculate
        if synonymous_site == 0 or nonsynonymous_site == 0:
            cn_cs = 0
        elif synonymous_count == 0:
            cs = (synonymous_count + 0.5)/(synonymous_site + 0.5)
            cn = nonsynonymous_count/nonsynonymous_site
            cn_cs = cn/cs
        else:
            cs = synonymous_count/synonymous_site
            cn = nonsynonymous_count/nonsynonymous_site
            cn_cs = cn/cs
        Cncs[enst] = cn_cs
        Cncs_p[enst] = p_value
        nonsynonymous_site = "{:.2f}".format(nonsynonymous_site)
        synonymous_site = "{:.2f}".format(synonymous_site)
        temp = '\t'.join(map(str,[gene_name,nonsynonymous_site,synonymous_site,nonsynonymous_count,synonymous_count,nonsynonymou_site_with_mutation,synonymous_site_with_mutation,cn,cs,cn_cs,p_value])) + '\n'
        OUT.write(temp)
    OUT.close()
    print("Finish calculate CN/CS for each gene and print result!")
    return Cncs,Cncs_p,Aa_mut_count,L_N,L_S,S_N,S_S,Mut_N,Mut_S

def CNCSPipe(mut_file:str,alpha:float,out_cncs:str,out_gene:str,out_site:str,out_two_component_cncs:str,out_H_test:str,fas_cds:list=[],cds_file = None, fas_file = None):
    '''
    Run the pipeline for CN/CS-calculator.\n
    mut_file:input\n
    alpha:flaot or you can set it as inf\n
    out_cncs:file dir.write cncs result to it\n
    out_gene:file dir.write gene with at least a driver site\n
    out_site:file dir.write driver site\n
    out_two_component_cncs:file dir.write cncs ratio of passenger sites and driver sites\n
    out_H_test:file dir.write H_test result to it.\n
    fas_cds: the output of prepare_fas_cds_file function\n
    cds_file: cds file path, if fas_cds is not empty list, this parameter will be ignored\n
    fas_file: fasta file path, if fas_cds is not empty list, this parameter will be ignored\n
    '''
    def gamma(n):
        return ga(n+1)
    ## calculate likelihood ratio test value
    def LRT(max_hit,eta,m1,m0,mean):
        global e
        E = dec(e)
        lg = lambda x: dec(log(dec(x)))
        mean = dec(mean)
        m0 = dec(m0)
        m1 = dec(m1)
        eta = dec(eta)
        if max_hit > 50: max_hit = 50
        fun_single = dec(1)
        fun_two = dec(1)
        for i in range(max_hit+1):
            fun_single += lg((mean ** i)/((E ** mean) * dec(gamma(i))))
            fun_two += lg((1-eta)*(m0**i)/((E ** m0) * dec(gamma(i))) + eta * (m1**i)/((E**m1) * dec(gamma(i))))
        return float(2*(fun_two-fun_single))

    if fas_cds == []:
        Gene_id,Aa_length,Gene_name,Synonymous,Non_synonymous,Nonsynonymous_site_with_mutation_initial,Synonymous_site_with_mutation_initial,Mut_missense_initial,Mut_silent_initial = prepare_fas_cds_file(cds_file, fas_file)
    else:
        Gene_id,Aa_length,Gene_name,Synonymous,Non_synonymous,Nonsynonymous_site_with_mutation_initial,Synonymous_site_with_mutation_initial,Mut_missense_initial,Mut_silent_initial = fas_cds
    
    print("START CAL CNCS.")
    Cncs,Cncs_p,Aa_mut_count,L_N,L_S,S_N,S_S,Mut_N,Mut_S = cncs(mut_file,out_cncs,fas_cds = [Gene_id,Aa_length,Gene_name,Synonymous,Non_synonymous,Nonsynonymous_site_with_mutation_initial,Synonymous_site_with_mutation_initial,Mut_missense_initial,Mut_silent_initial])
    print("FINISHED CNCS CAL")

    OUT1 = open_file(out_gene,mode='w')
    OUT2 = open_file(out_site,mode='w')
    OUT3 = open_file(out_two_component_cncs,mode='w')
    OUT4 = open_file(out_H_test,mode='w')
    OUT1.write("Gene_id\tTranscript_id\tGene_name\tChr\tMax_hit\tCount\tProtein_length\tDistribution\tMut_sites\tf0\tmean\tvariance\tm0\tm1\teta\tLRT_p-value\tQ(z)\n")
    OUT2.write("Gene_id\tTranscript_id\tGene_name\tChr\tProtein_position\tz\tQ(z)\tProtein_mutation\tM(z)\tOmega(z)\n")
    OUT3.write("Gene_id\tTranscript_id\tGene_name\tChr\tOmega\tOmega_p\tOmegadri\tOmegapass\tOmegapass_p\n")
    OUT4.write("Gene_id\tTranscript_id\tGene_name\tChr\tCN/CS\t1-H\tHtest_p_value\n")

    Count = dict()
    Site_mut = []
    Site_mut_list = dict()
    Max = dict()
    Hit = dict()

    Omega0 = {}
    Omega1 = {}
    Omega0_p = {}
    ## Read mutation records
    for key in Aa_mut_count.keys():
        mut = key.split('\t')
        gene = mut[1]
        aa_mut = mut[2]
        if gene in Gene_id.keys():
            pattern = compile(r'[a-zA-Z](\d+)[a-zA-z]')
            pos = pattern.findall(aa_mut)
            if pos:
                pos = pos[0]
            count = Aa_mut_count[key]
            site_mut = gene + ':' + pos
            if not pos: print(key)
            Site_mut_list[site_mut] = Site_mut_list.get(site_mut,'') + aa_mut + ':' + str(count) + ';'
            if site_mut not in Count.keys():
                    Count[site_mut] = count
                    Site_mut.append(site_mut)
            else:
                Count[site_mut] = Count.get(site_mut,0) + count
    print("Read mutation records!")

    ## Sort site mutation
    for site_mut in Site_mut:
        mut_count = {}
        mutation = Site_mut_list[site_mut].split(';')
        for mutation_ele in mutation:
            if not mutation_ele: continue
            aa_mut = mutation_ele.split(':')[0]
            mut_count[aa_mut] = mut_count.get(aa_mut,'') + mutation_ele.split(':')[1]

        Site_mut_list[site_mut] = ''
        for aa_mut , _ in sorted(mut_count.items(),key=lambda x: x[1],reverse=True):
            Site_mut_list[site_mut] += aa_mut + ':' + mut_count[aa_mut] + ';'


    gene_list = []
    ## Mutation profile for each gene
    for site_mut in Site_mut:
        gene = site_mut.split(':')[0]
        pos = site_mut.split(':')[1]
        if gene not in Count.keys():
            Count[gene] = Count[site_mut]
            Max[gene] = Count[site_mut]
            Hit[gene] = str(Count[site_mut]) + ':' + pos + ';' 
            gene_list.append(gene)
        else:
            Count[gene] += Count[site_mut]
            Hit[gene] += str(Count[site_mut]) + ':' +pos + ';'
            if (Count[site_mut] > Max[gene]) : Max[gene] = Count[site_mut] 
    print("Count the mutation profile for each gene!")
    
    
    # sort mutation distribution
    for gene in gene_list:
        site_hit = dict()
        site_hit_list = Hit[gene].split(';')
        for hit in site_hit_list:
            if not hit : continue
            pos = hit.split(':')[1]
            site_hit[pos] = int(hit.split(':')[0])
        Hit[gene] = ""
        for pos,_ in sorted(site_hit.items(),key=lambda x: (-x[1],x[0])):
            Hit[gene] += str(site_hit[pos]) + ':' + pos + ';'

    ## Calculate paramaters
    for gene in gene_list:
        id = Gene_id[gene]
        count_int = Count[gene]
        count_trim = int(count_int)
        aa = eval(Aa_length[gene])
        maxhit = Max[gene]
        dis = Hit[gene]
        dis_list = dis.split(';')[:-1]
        aa_trim = int(aa)
        z1,z0 = 0,0
        
        for i in range(len(dis_list)):
            dis1 = int(dis_list[i].split(':')[0])
            if dis1 > 100:
                count_trim -= dis1
                aa_trim -= 1

        h,q,m = [0] * len(dis_list),[0] * len(dis_list),[0] * len(dis_list)
        Omega = [0] * len(dis_list)
        hit = 0

        if isinstance(aa,int):
            mean = count_trim / aa_trim
            num = len(dis_list) 
            nohit = aa - num
            if nohit == 0:
                continue
            else:
                f0 = nohit / aa
                m0 = -log(f0)
                var = nohit * (mean **2)
                for i in range(len(dis_list)):
                    dis1 = int(dis_list[i].split(':')[0])
                    if dis1 < 100: var += (dis1 - mean) **2 
                var = var / (aa_trim - 1)
                if alpha == 'inf':
                    m1 = mean + (var - mean)/(mean - m0)
                else:
                    m1 = (mean + (m0/alpha) + (var - mean)/(mean - m0))/(1+1/alpha)
                eta = (mean - m0) / (m1 - m0)
                if m1 < 0 or eta < 0 or eta > 1:
                    continue

                ################################################
                # Calculate omega0 & omega1
                omega = Cncs[gene]
                omega_p = Cncs_p[gene]
                # parameters for omega0 and omega1
                mut_n,mut_s,ln,ls,sn,ss = Mut_N[gene],Mut_S[gene],L_N[gene],L_S[gene],S_N[gene],S_S[gene]

                # calculate omega0 (SN/SS)/(LN/LS)
                if ls == 0 or ln == 0:
                    omega0 = 'NA'
                elif ss == 0:
                    omega0 = (sn/ln)/((ss+0.5)/(ls+0.5))
                else:
                    omega0 = (sn/ln)/(ss/ls)
                # calculate omega0_p
                try:
                    chi2_con_array_omega0 = array([[ln,ls], [sn, ss]])
                    chi2_value_omega0, p_value_omega0, df_omega0, expected_omega0 = chi2_contingency(chi2_con_array_omega0)
                except:
                    p_value_omega0 = 'NA'
                # store
                Omega0[gene] = omega0
                Omega0_p[gene] = p_value_omega0

                # calculate omega1
                try:
                    omega1_numerator = ((mut_n/mut_s) - (1-eta) * (sn/ss))/eta
                    omega1 = omega1_numerator / (ln/ls)
                    if omega1 < 0:
                        omega1 = 'NA'
                except ZeroDivisionError:
                    omega1 = 'NA'
                # store
                Omega1[gene] = omega1

                if (omega0 != 'NA') and (omega1 != 'NA'):
                    temp_two_comp = [id,omega,omega_p,omega1,omega0,p_value_omega0]
                    s_two_comp = '\t'.join(list(map(str,temp_two_comp))) + '\n'
                    OUT3.write(s_two_comp)
                    
                ################################################
                #H test
                H = (var - mean) / (var + mean * (mean - 1))
                H0=1-H
                # statistical test for H
                try:
                    chi2_con_array_H = array([[mut_n,mut_s], [ln*(1-H), ls]])
                    chi2_value_H, Htest_p_value, df_H, expected_H = chi2_contingency(chi2_con_array_H)
                except:
                    Htest_p_value = 'NA'
                # mut_n_hat = (mut_n + mut_s) * ((ln*(1-H))/(ls+ln*(1-H)))
                # mut_s_hat = (mut_n + mut_s) * (ls/(ls+ln*(1-H)))
                # ln0_hat = ln - mut_n_hat
                # ls0_hat = ls - mut_s_hat
                # ln0 = ln - mut_n
                # ls0 = ls - mut_s
                # # caculate chi-square value
                # Htest_chi_square = ((mut_n - mut_n_hat) ** 2 / mut_n_hat) + ((mut_s - mut_s_hat) ** 2 / mut_s_hat) + ((ln0 - ln0_hat) ** 2 / ln0_hat) + ((ls0 - ls0_hat) ** 2 / ls0_hat)
                # Htest_p_value = chi2.sf(Htest_chi_square,3)
                temp_H = [id,omega,H0,Htest_p_value]
                s_H = '\t'.join(list(map(str,temp_H)))+'\n'
                OUT4.write(s_H)

            ################################################
            # candris site & gene calculation for gene that has at least one site with mutation for three or more times
            if maxhit > 2:
                qz = ''
                lrt = LRT(maxhit,eta,m1,m0,mean)
                p = chi2.sf(lrt,1)
                for i in range(len(dis_list)):
                    dis1,pos = dis_list[i].split(':')[0:2]
                    dis0 = dis_list[i-1].split(':')[0] 
                    z = 50 if int(dis1) > 50 else int(dis1)
                    if alpha == 'inf':
                        ratio = log(eta/(1-eta)) + z*log(m1/m0)-(m1-m0)
                        q[i] = 1 - 1/float(1+dec(e)**dec(ratio))
                        try:
                            Omega[i] = (1-q[i])*omega0+q[i]*omega1
                        except TypeError:
                            Omega[i] = 'NA'
                    else:
                        h[i] = (gamma(z+alpha-1)/gamma(alpha-1)) * e**(z*(log(m1)-log(m0))+alpha*log(alpha)+m0-(z+alpha)*log(m1+alpha))
                        q[i] = eta * h[i]/(1+eta*(h[i]-1))
                        m[i] = (1-q[i])*m0+q[i]*((z+alpha)/(m1+alpha))*m1
                        try:
                            Omega[i] = (1-q[i])*omega0 + q[i]*((z+alpha)/(m1+alpha))*omega1
                        except TypeError:
                            Omega[i] = 'NA'
                    site_mut = gene + ':' + pos
                    OUT2.write('\t'.join([id,pos,dis1,str(q[i]),Site_mut_list[site_mut],str(m[i]),str(Omega[i])]) + '\n')

                    ## The distribution of driver sites in each gene
                    if len(dis_list) == 1:
                        hit = 1
                        qz += '{0}:{1}:{2};'.format(q[i],dis1,hit)
                    elif i == 0:
                        hit = 1
                    elif (i+1 < len(dis_list)) and (dis1 == dis0):
                        hit += 1
                    elif (i+1 < len(dis_list)) and (dis1 != dis0):
                        qz += '{0}:{1}:{2};'.format(q[i-1],dis0,hit)
                        hit = 1
                    elif (i+1 == len(dis_list)) and (dis1 == dis0):
                        hit += 1
                        qz += '{0}:{1}:{2};'.format(q[i],dis1,hit)
                    elif (i+1 == len(dis_list)) and (dis1 != dis0):
                        qz += '{0}:{1}:{2};'.format(q[i-1],dis0,hit)
                        hit = 1
                        qz += '{0}:{1}:{2};'.format(q[i],dis1,hit)

                if alpha == 'inf':
                    q0 = 1 - 1/(1 + (eta/(1-eta))* e**(m0-m1))
                else:
                    h0 = e**(alpha*log(alpha)+m0-(alpha)*log(m1+alpha))
                    q0 = eta*h0/(1+eta*(h0-1))

                qz += str(q0) + ":0:" + str(nohit) + ';'          
                temp_can = [id,maxhit,count_int,aa,dis,num,f0,mean,var,m0,m1,eta,p,qz]
                s_can = '\t'.join(list(map(str,temp_can))) + '\n'
                OUT1.write(s_can)
    print("Finish calculate paramates!")
    OUT1.close()
    OUT2.close()
    OUT3.close()
    OUT4.close()

def H_test(mut_file:str,alpha:float,out_cncs:str,out_H_test:str,fas_cds:list=[],cds_file = None, fas_file = None):
    '''
    Run H_test\n
    mut_file:input\n
    alpha:flaot or you can set it as inf\n
    out_cncs:file dir.write cncs result to it\n
    out_H_test:file dir.write H_test result to it\n
    fas_cds:the output of prepare_fas_cds_file function\n
    cds_file: cds file path, if fas_cds is not empty list, this parameter will be ignored\n
    fas_file: fasta file path, if fas_cds is not empty list, this parameter will be ignored\n
    '''
    if fas_cds == []:
        Gene_id,Aa_length,Gene_name,Synonymous,Non_synonymous,Nonsynonymous_site_with_mutation_initial,Synonymous_site_with_mutation_initial,Mut_missense_initial,Mut_silent_initial = prepare_fas_cds_file(cds_file, fas_file)
    else:
        Gene_id,Aa_length,Gene_name,Synonymous,Non_synonymous,Nonsynonymous_site_with_mutation_initial,Synonymous_site_with_mutation_initial,Mut_missense_initial,Mut_silent_initial = fas_cds
    
    print("START CAL CNCS.")
    Cncs,Cncs_p,Aa_mut_count,L_N,L_S,S_N,S_S,Mut_N,Mut_S = cncs(mut_file,out_cncs,fas_cds = [Gene_id,Aa_length,Gene_name,Synonymous,Non_synonymous,Nonsynonymous_site_with_mutation_initial,Synonymous_site_with_mutation_initial,Mut_missense_initial,Mut_silent_initial])
    print("FINISHED CNCS CAL")

    OUT1 = open_file(out_H_test,mode='w')
    OUT1.write("Gene_id\tTranscript_id\tGene_name\tChr\tCN/CS\t1-H\tHtest_p_value\n")

    Count = dict()
    Site_mut = []
    Site_mut_list = dict()
    Max = dict()
    Hit = dict()

    ## Read mutation records
    for key in Aa_mut_count.keys():
        mut = key.split('\t')
        gene = mut[1]
        aa_mut = mut[2]
        if gene in Gene_id.keys():
            pattern = compile(r'[a-zA-Z](\d+)[a-zA-z]')
            pos = pattern.findall(aa_mut)
            if pos:
                pos = pos[0]
            count = Aa_mut_count[key]
            site_mut = gene + ':' + pos
            if not pos: print(key)
            Site_mut_list[site_mut] = Site_mut_list.get(site_mut,'') + aa_mut + ':' + str(count) + ';'
            if site_mut not in Count.keys():
                    Count[site_mut] = count
                    Site_mut.append(site_mut)
            else:
                Count[site_mut] = Count.get(site_mut,0) + count
    print("Read mutation records!")

    ## Sort site mutation
    for site_mut in Site_mut:
        mut_count = {}
        mutation = Site_mut_list[site_mut].split(';')
        for mutation_ele in mutation:
            if not mutation_ele: continue
            aa_mut = mutation_ele.split(':')[0]
            mut_count[aa_mut] = mut_count.get(aa_mut,'') + mutation_ele.split(':')[1]

        Site_mut_list[site_mut] = ''
        for aa_mut , _ in sorted(mut_count.items(),key=lambda x: x[1],reverse=True):
            Site_mut_list[site_mut] += aa_mut + ':' + mut_count[aa_mut] + ';'


    gene_list = []
    ## Mutation profile for each gene
    for site_mut in Site_mut:
        gene = site_mut.split(':')[0]
        pos = site_mut.split(':')[1]
        if gene not in Count.keys():
            Count[gene] = Count[site_mut]
            Max[gene] = Count[site_mut]
            Hit[gene] = str(Count[site_mut]) + ':' + pos + ';' 
            gene_list.append(gene)
        else:
            Count[gene] += Count[site_mut]
            Hit[gene] += str(Count[site_mut]) + ':' +pos + ';'
            if (Count[site_mut] > Max[gene]) : Max[gene] = Count[site_mut] 
    print("Count the mutation profile for each gene!")
    
    
    # sort mutation distribution
    for gene in gene_list:
        site_hit = dict()
        site_hit_list = Hit[gene].split(';')
        for hit in site_hit_list:
            if not hit : continue
            pos = hit.split(':')[1]
            site_hit[pos] = int(hit.split(':')[0])
        Hit[gene] = ""
        for pos,_ in sorted(site_hit.items(),key=lambda x: (-x[1],x[0])):
            Hit[gene] += str(site_hit[pos]) + ':' + pos + ';'

    ## Calculate paramaters
    for gene in gene_list:
        id = Gene_id[gene]
        count_int = Count[gene]
        count_trim = int(count_int)
        aa = eval(Aa_length[gene])
        maxhit = Max[gene]
        dis = Hit[gene]
        dis_list = dis.split(';')[:-1]
        aa_trim = int(aa)
        z1,z0 = 0,0
        
        for i in range(len(dis_list)):
            dis1 = int(dis_list[i].split(':')[0])
            if dis1 > 100:
                count_trim -= dis1
                aa_trim -= 1

        h,q,m = [0] * len(dis_list),[0] * len(dis_list),[0] * len(dis_list)
        Omega = [0] * len(dis_list)
        hit = 0

        if isinstance(aa,int):
            mean = count_trim / aa_trim
            num = len(dis_list) 
            nohit = aa - num
            if nohit == 0:
                continue
            else:
                f0 = nohit / aa
                m0 = -log(f0)
                var = nohit * (mean **2)
                for i in range(len(dis_list)):
                    dis1 = int(dis_list[i].split(':')[0])
                    if dis1 < 100: var += (dis1 - mean) **2 
                var = var / (aa_trim - 1)
                if alpha == 'inf':
                    m1 = mean + (var - mean)/(mean - m0)
                else:
                    m1 = (mean + (m0/alpha) + (var - mean)/(mean - m0))/(1+1/alpha)
                eta = (mean - m0) / (m1 - m0)
                if m1 < 0 or eta < 0 or eta > 1:
                    continue

                # used parameters
                omega = Cncs[gene]
                mut_n,mut_s,ln,ls,sn,ss = Mut_N[gene],Mut_S[gene],L_N[gene],L_S[gene],S_N[gene],S_S[gene]
                # Htest
                H = (var - mean) / (var + mean * (mean - 1))
                H0=1-H
                # statistical test for H
                try:
                    chi2_con_array_H = array([[mut_n,mut_s], [ln*(1-H), ls]])
                    chi2_value_H, Htest_p_value, df_H, expected_H = chi2_contingency(chi2_con_array_H)
                except:
                    Htest_p_value = 'NA'
                
                temp_H = [id,omega,H0,Htest_p_value]
                s_H = '\t'.join(list(map(str,temp_H)))+'\n'
                OUT1.write(s_H)
    print("Finish calculate paramates!")
    OUT1.close()

def two_component_cncs(mut_file:str,alpha:float,out_cncs:str,out_two_component_cncs:str,fas_cds:list=[],cds_file = None, fas_file = None):
    '''
    Run two_component_cncs.\n
    mut_file:input\n
    alpha:flaot or you can set it as inf\n
    out_cncs:file dir.write cncs result to it\n
    out_two_component_cncs:file dir.write cncs ratio of passenger sites and driver sites\n
    fas_cds:the output of prepare_fas_cds_file function\n
    cds_file: cds file path, if fas_cds is not empty list, this parameter will be ignored\n
    fas_file: fasta file path, if fas_cds is not empty list, this parameter will be ignored\n
    '''
    if fas_cds == []:
        Gene_id,Aa_length,Gene_name,Synonymous,Non_synonymous,Nonsynonymous_site_with_mutation_initial,Synonymous_site_with_mutation_initial,Mut_missense_initial,Mut_silent_initial = prepare_fas_cds_file(cds_file, fas_file)
    else:
        Gene_id,Aa_length,Gene_name,Synonymous,Non_synonymous,Nonsynonymous_site_with_mutation_initial,Synonymous_site_with_mutation_initial,Mut_missense_initial,Mut_silent_initial = fas_cds
    
    print("START CAL CNCS.")
    Cncs,Cncs_p,Aa_mut_count,L_N,L_S,S_N,S_S,Mut_N,Mut_S = cncs(mut_file,out_cncs,fas_cds = [Gene_id,Aa_length,Gene_name,Synonymous,Non_synonymous,Nonsynonymous_site_with_mutation_initial,Synonymous_site_with_mutation_initial,Mut_missense_initial,Mut_silent_initial])
    print("FINISHED CNCS CAL")

    OUT1 = open_file(out_two_component_cncs,mode='w')
    OUT1.write("Gene_id\tTranscript_id\tGene_name\tChr\tOmega\tOmega_p\tOmegadri\tOmegapass\tOmegapass_p\n")

    Count = dict()
    Site_mut = []
    Site_mut_list = dict()
    Max = dict()
    Hit = dict()

    Omega0 = {}
    Omega1 = {}
    Omega0_p = {}
    ## Read mutation records
    for key in Aa_mut_count.keys():
        mut = key.split('\t')
        gene = mut[1]
        aa_mut = mut[2]
        if gene in Gene_id.keys():
            pattern = compile(r'[a-zA-Z](\d+)[a-zA-z]')
            pos = pattern.findall(aa_mut)
            if pos:
                pos = pos[0]
            count = Aa_mut_count[key]
            site_mut = gene + ':' + pos
            if not pos: print(key)
            Site_mut_list[site_mut] = Site_mut_list.get(site_mut,'') + aa_mut + ':' + str(count) + ';'
            if site_mut not in Count.keys():
                    Count[site_mut] = count
                    Site_mut.append(site_mut)
            else:
                Count[site_mut] = Count.get(site_mut,0) + count
    print("Read mutation records!")

    ## Sort site mutation
    for site_mut in Site_mut:
        mut_count = {}
        mutation = Site_mut_list[site_mut].split(';')
        for mutation_ele in mutation:
            if not mutation_ele: continue
            aa_mut = mutation_ele.split(':')[0]
            mut_count[aa_mut] = mut_count.get(aa_mut,'') + mutation_ele.split(':')[1]

        Site_mut_list[site_mut] = ''
        for aa_mut , _ in sorted(mut_count.items(),key=lambda x: x[1],reverse=True):
            Site_mut_list[site_mut] += aa_mut + ':' + mut_count[aa_mut] + ';'


    gene_list = []
    ## Mutation profile for each gene
    for site_mut in Site_mut:
        gene = site_mut.split(':')[0]
        pos = site_mut.split(':')[1]
        if gene not in Count.keys():
            Count[gene] = Count[site_mut]
            Max[gene] = Count[site_mut]
            Hit[gene] = str(Count[site_mut]) + ':' + pos + ';' 
            gene_list.append(gene)
        else:
            Count[gene] += Count[site_mut]
            Hit[gene] += str(Count[site_mut]) + ':' +pos + ';'
            if (Count[site_mut] > Max[gene]) : Max[gene] = Count[site_mut] 
    print("Count the mutation profile for each gene!")
    
    
    # sort mutation distribution
    for gene in gene_list:
        site_hit = dict()
        site_hit_list = Hit[gene].split(';')
        for hit in site_hit_list:
            if not hit : continue
            pos = hit.split(':')[1]
            site_hit[pos] = int(hit.split(':')[0])
        Hit[gene] = ""
        for pos,_ in sorted(site_hit.items(),key=lambda x: (-x[1],x[0])):
            Hit[gene] += str(site_hit[pos]) + ':' + pos + ';'

    ## Calculate paramaters
    for gene in gene_list:
        id = Gene_id[gene]
        count_int = Count[gene]
        count_trim = int(count_int)
        aa = eval(Aa_length[gene])
        maxhit = Max[gene]
        dis = Hit[gene]
        dis_list = dis.split(';')[:-1]
        aa_trim = int(aa)
        z1,z0 = 0,0
        
        for i in range(len(dis_list)):
            dis1 = int(dis_list[i].split(':')[0])
            if dis1 > 100:
                count_trim -= dis1
                aa_trim -= 1

        h,q,m = [0] * len(dis_list),[0] * len(dis_list),[0] * len(dis_list)
        Omega = [0] * len(dis_list)
        hit = 0

        if isinstance(aa,int):
            mean = count_trim / aa_trim
            num = len(dis_list) 
            nohit = aa - num
            if nohit == 0:
                continue
            else:
                f0 = nohit / aa
                m0 = -log(f0)
                var = nohit * (mean **2)
                for i in range(len(dis_list)):
                    dis1 = int(dis_list[i].split(':')[0])
                    if dis1 < 100: var += (dis1 - mean) **2 
                var = var / (aa_trim - 1)
                if alpha == 'inf':
                    m1 = mean + (var - mean)/(mean - m0)
                else:
                    m1 = (mean + (m0/alpha) + (var - mean)/(mean - m0))/(1+1/alpha)
                eta = (mean - m0) / (m1 - m0)
                if m1 < 0 or eta < 0 or eta > 1:
                    continue

                ################################################
                # Calculate omega0 & omega1
                omega = Cncs[gene]
                omega_p = Cncs_p[gene]
                # parameters for omega0 and omega1
                mut_n,mut_s,ln,ls,sn,ss = Mut_N[gene],Mut_S[gene],L_N[gene],L_S[gene],S_N[gene],S_S[gene]
                
                # calculate omega0 (SN/SS)/(LN/LS)
                if ls == 0 or ln == 0:
                    omega0 = 'NA'
                elif ss == 0:
                    omega0 = (sn/ln)/((ss+0.5)/(ls+0.5))
                else:
                    omega0 = (sn/ln)/(ss/ls)
                # calculate omega0_p
                try:
                    chi2_con_array_omega0 = array([[ln,ls], [sn, ss]])
                    chi2_value_omega0, p_value_omega0, df_omega0, expected_omega0 = chi2_contingency(chi2_con_array_omega0)
                except:
                    p_value_omega0 = 'NA'
                # store
                Omega0[gene] = omega0
                Omega0_p[gene] = p_value_omega0

                # calculate omega1
                try:
                    omega1_numerator = ((mut_n/mut_s) - (1-eta) * (sn/ss))/eta
                    omega1 = omega1_numerator / (ln/ls)
                except ZeroDivisionError:
                    omega1 = 'NA'
                # store
                Omega1[gene] = omega1
                
                if (omega0 != 'NA') and (omega1 != 'NA'):
                    temp_two_comp = [id,omega,omega_p,omega1,omega0,p_value_omega0]
                    s_two_comp = '\t'.join(list(map(str,temp_two_comp))) + '\n'
                    OUT1.write(s_two_comp)
    print("Finish calculate paramates!")
    OUT1.close()
