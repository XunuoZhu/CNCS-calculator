# Description

CN/CS-calculator is a site-specific selection pressure analysis method, which consists of three integral modules: CN/CS-Pos, CN/CS-H and CN/CS-Neg. CN/CS-calculator leverages the information of mutation frequency, equipping it with the ability to decompose genes into distinct driver and passenger components, thereby capturing selective pressures that differ across various gene sites. This innovative strategy enables a comprehensive analysis of selection pressures across the cancer genome.

Version: 2.0; 

Date: 2024.1.25.

# Contact

Xun Gu		xgu@iastate.edu

Zhou Zhan	zhanzhou@zju.edu.cn

# Analysis steps

1. Choose specific version of Ensembl reference data (eg. TCGA(75_GRCh37.p13))
2. Data preprocessing for input (eg. collect mutation data from TCGA(.maf))
3. Site-specific selection pressure prediction steps as follows:
   1. choose calculation model (collect mutation profiles from TCGA data for mutation-profile based model)
   2. choose specific version of sequence data (.fas)
   3. load genetic code file
   4. perform CNCSPipe
   5. output files in your path

# Usage

Files required(information of protein coding genes, mutation rate profile file, protein coding sequences, genetic code) is included in the package. You can change the files if necessary.
	
## Install

```shell
tar -zxvf CNCScalculator-2.0.tar.gz
cd CNCScalculator
python setup.py install
```

## Example:

```python
from CNCScalculator import CNCSPipe,cncs,H_test,two_component_cncs
## Help information
print(help(CNCSPipe))
print(help(cncs))
print(help(H_test))
print(help(two_component_cncs))
## CNCSPipe(mut_file:str,alpha:float,out_cncs:str,out_gene:str,out_site:str,out_two_component_cncs:str,out_H_test:str)
CNCSPipe("mutation_data.maf",1.51,"cncs.txt","TMM_gene.txt","TMM_site.txt","two_compo_cncs.txt","H_test.txt")

## cncs(mut_file:str,outputfile:str,arg:list = [])
cncs("mutation_data.maf","cncs.txt")

## H_test(mut_file:str,alpha:float,out_cncs:str,out_H_test:str)
H_test("mutation_data.maf",1.51,"cncs.txt","H_test.txt")

## two_component_cncs(mut_file:str,alpha:float,out_cncs:str,out_two_component_cncs:str)
two_component_cncs("mutation_data.maf",1.51,"cncs.txt","two_compo_cncs.txt")
```

## Input File

Mini Sample:

```
Hugo_Symbol	Consequence	Variant_Type	HGVSp_Short	Gene	Feature	HGVSc
TACC2	missense_variant	SNP	p.T38M	ENSG00000138162	ENST00000369005	c.113C>T
JAKMIP3	synonymous_variant	SNP	p.D723D	ENSG00000188385	ENST00000298622	c.2169C>T
PANX3	missense_variant	SNP	p.R296Q	ENSG00000154143	ENST00000284288	c.887G>A
SPI1	missense_variant	SNP	p.P127T	ENSG00000066336	ENST00000227163	c.379C>A
NAALAD2	missense_variant	SNP	p.R65C	ENSG00000077616	ENST00000534061	c.193C>T
FAT3	synonymous_variant	SNP	p.P3444P	ENSG00000165323	ENST00000298047	c.10332G>A
MTERFD3	missense_variant	SNP	p.L213S	ENSG00000120832	ENST00000552029	c.638T>C
BTBD11	missense_variant	SNP	p.E770K	ENSG00000151136	ENST00000280758	c.2308G>A
NOS1	upstream_gene_variant	SNP	.	ENSG00000089250	ENST00000338101	.
```

## Output File

Result of cncs.txt :

```markdown
## Title
Gene_id	nonsynonymous_site	synonymous_site	nonsynonymous_count	synonymous_count	nonsynonymous_site_with_mutation	CN	CS	CN/CS p_value
## Explanation
Gene_id: "ensg|enst|gene|chr|length_of_coding_region"
p_value: the p value of chi-continency of categorical variables:(nonsynonymous_site,synonymous_site) and (nonsynonymous_count,synonymous_count)
```

Result of TMM_gene.txt:

```markdown
## Title
Gene_id	Transcript_id	Gene_name	Chr	Max_hit	Count	Protein_length	Distribution	Mut_sites	f0	mean	variance	m0	m1	eta	LRT_p-value	Q(z)
## Explanation
f0: nohit in aa length
mean: mean of z(count in aa length)
eta: the probability of any site being a driver
```

Result of TMM_site.txt:

```markdown
## Title
Gene_id	Transcript_id	Gene_name	Chr	Protein_position	z	Q(z)	Protein_mutation	M(z)	Omega(z)
## Explanation
z: the number of somatic missense mutations
Q(z): the posterior probability of an amino acid site being a cancer-driver when the number (z) of somatic mutations is observed at this site
M(z): the posterior mean of recurrent mutations at this site
Omega(z): the posterior mean of CN/CS ratio at this site
## mini sample
ENSG00000088256	ENST00000078429	GNA11	19	209	34	1.0	Q209L:34;	33.06147756474409	8900.207358728772
ENSG00000088256	ENST00000078429	GNA11	19	166	1	0.03295076422042006	R166H:1;	0.07036662907632019	18.94281913709293
ENSG00000088256	ENST00000078429	GNA11	19	183	1	0.03295076422042006	R183C:1;	0.07036662907632019	18.94281913709293
ENSG00000115524	ENST00000335508	SF3B1	2	625	14	1.0	R625H:8;R625C:6;	13.020953605053787	12411.247750114364
ENSG00000115524	ENST00000335508	SF3B1	2	666	2	0.9718624368788241	K666T:2;	2.53100158441306	2412.4874930734936
```

Result of H_test.txt:

```markdown
## Title
Gene_id	Transcript_id	Gene_name	Chr	CN/CS	H_value	H_test
## Explanation
H_value: a relative measure of rate variation among sites
H_test: CN/CS>1-H for major and minor driver genes
```

Result of two_compo_cncs.txt:

```markdown
## Title
Gene_id	Transcript_id	Gene_name	Chr	Omega	Omega_p	Omega1	Omega0	Omega0_p
## Explanation
Omega: CN/CS ratio of this gene
Omega_p: the p value of chi-continency assert whether null hypothesis of CN/CS=1 can be statistically rejected.
Omega1: CN/CS ratio of driver sites of this gene
Omega0: CN/CS ratio of passenger sites of this gene
Omega0_p: the p value of chi-continency assert whether null hypothesis of Omega0=1 can be statistically rejected

##  Build-in Data File 

```
CNCScalculator-2.0
│	README.md
│	MANIFEST.in
│  setup.py    
│	...
└───CNCScalculator
│   │   __init__.py
│   │   CNCScalculator.py
│   └───data
│       │   genetic_code.txt
│       │   human_TCGA_match_GRCh37.fas
│       │   pancancer_mutation_rate.txt
│       │   TCGA_match_GRCh37.cds
└───CNCScalculator.egg-info
    │	....
```

**genetic_code.txt**                       symonymous or non-synonymous in codon mutation

**human_TCGA_match_GRCh37.fas**            fasta file of human TCGS(GRCh37)

**pancancer_mutation_rate.txt**	          Pancancer_rate of Ref mutate to Alt	

**TCGA_match_GRCh37.cds**	                cds information of genes