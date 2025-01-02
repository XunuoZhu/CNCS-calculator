# *C<sub>N</sub>/C<sub>S</sub>*-calculator

*C<sub>N</sub>/C<sub>S</sub>*-calculator is a site-specific selection pressure analysis method designed to capture selective pressures acting across different gene sites. The rationale behind *C<sub>N</sub>/C<sub>S</sub>*-calculator lies in the recognition that not all variants within a gene contribute equally to oncogenic processes, with some variants driving cancer progression while others are constrained by functional requirements.

*C<sub>N</sub>/C<sub>S</sub>*-calculator consists of three modules. With the input cancer somatic mutation information, the first module of *C<sub>N</sub>/C<sub>S</sub>*-calculator computes the selection pressure at the gene level, where *C<sub>N</sub>/C<sub>S</sub>* \> 1 denotes gene under positive selection while 1-H \< *C<sub>N</sub>/C<sub>S</sub>* \< 1 denotes that only several sites are under weak positive selection. Then, for each positively selected gene, the mutational landscape of different amino acid sites is inscribed using a two-component mixture module, which models the proportion of driver component ($\eta$) and passenger component ($1 - \eta$) in a gene. Finally, for the site components computed by the second module, the third module calculates site component-specific selection pressure *Ω<sub>driver</sub>* and *Ω<sub>passenger</sub>*, respectively.

## Install

```shell
tar -zxvf CNCScalculator-2.0.tar.gz
cd CNCScalculator
python setup.py install
```

## Analysis Steps

1. Choose specific version of Ensembl reference data (eg. TCGA(75_GRCh37.p13))
2. Data preprocessing for input (eg. collect mutation data from TCGA(.maf))
3. Site-specific selection pressure prediction steps as follows:
   1. choose calculation model (collect mutation profiles from TCGA data for mutation-profile based model)
   2. choose specific version of sequence data (.fas)
   3. load genetic code file
   4. perform CNCSPipe
   5. output files in your path

## Usage

Files required (information of protein coding genes, mutation rate profile file, protein coding sequences, genetic code) is included in the package. You can change the files if necessary.

### Example

```python
from CNCScalculator import CNCSPipe,cncs,H_test,two_component_cncs
## Help information
print(help(CNCSPipe))
print(help(cncs))
print(help(H_test))
print(help(two_component_cncs))
## CNCSPipe(mut_file:str,alpha:float,out_cncs:str,out_gene:str,out_site:str,out_two_component_cncs:str,out_H_test:str)
CNCSPipe("Mutation_file.maf",1.51,"cncs.txt","gene.txt","site.txt","two_compo_cncs.txt","H_test.txt")

## cncs(mut_file:str,outputfile:str,arg:list = [])
cncs("Mutation_file.maf","cncs.txt")

## H_test(mut_file:str,alpha:float,out_cncs:str,out_H_test:str)
H_test("Mutation_file.maf",1.51,"cncs.txt","H_test.txt")

## two_component_cncs(mut_file:str,alpha:float,out_cncs:str,out_two_component_cncs:str)
two_component_cncs("Mutation_file.maf",1.51,"cncs.txt","two_compo_cncs.txt")
```

### Input File

Mini Sample:

| **Hugo_Symbol** | **Consequence**      | **Variant_Type**      | **HGVSp_Short** | **Gene**          | **Feature**       | **HGVSc**    |
|-----------------|----------------------|-----------------------|-----------------|-------------------|-------------------|--------------|
| TACC2          | missense_variant     | SNP                   | p.T38M          | ENSG00000138162   | ENST00000369005   | c.113C>T     |
| JAKMIP3        | synonymous_variant   | SNP                   | p.D723D         | ENSG00000188385   | ENST00000298622   | c.2169C>T    |
| PANX3          | missense_variant     | SNP                   | p.R296Q         | ENSG00000154143   | ENST00000284288   | c.887G>A     |
| SPI1           | missense_variant     | SNP                   | p.P127T         | ENSG00000066336   | ENST00000227163   | c.379C>A     |
| NAALAD2        | missense_variant     | SNP                   | p.R65C          | ENSG00000077616   | ENST00000534061   | c.193C>T     |
| FAT3           | synonymous_variant   | SNP                   | p.P3444P        | ENSG00000165323   | ENST00000298047   | c.10332G>A   |
| MTERFD3        | missense_variant     | SNP                   | p.L213S         | ENSG00000120832   | ENST00000552029   | c.638T>C     |
| BTBD11         | missense_variant     | SNP                   | p.E770K         | ENSG00000151136   | ENST00000280758   | c.2308G>A    |
| NOS1           | upstream_gene_variant| SNP                   | .               | ENSG00000089250   | ENST00000338101   | .            |

### Output File

#### Result of ```cncs.txt```

- Title

| **Gene_id** | **nonsynonymous_site** | **synonymous_site** | **nonsynonymous_count** | **synonymous_count** | **nonsynonymous_site_with_mutation** | **CN** | **CS** | **CN/CS** | **p_value** |
|-------------|------------------------|---------------------|-------------------------|----------------------|---------------------------------------|--------|--------|-----------|------------|

- Explanation
    - CN/CS: the ratio of nonsynonymous mutation to synonymous mutation
    - p_value: the p value of chi-continency assert whether null hypothesis of $C_{N}/C_{S} = 1$ can be statistically rejected

#### Result of ```gene.txt```

- Title

| **Gene_id** | **Transcript_id** | **Gene_name** | **Chr** | **Max_hit** | **Count** | **Protein_length** | **Distribution** | **Mut_sites** | **f0** | **mean** | **variance** | **m0** | **m1** | **eta** | **LRT_p-value** | **Q(z)** |
|-------------|-------------------|---------------|---------|-------------|-----------|-------------------|------------------|--------------|-------|---------|------------|-------|-------|--------|----------------|----------|

- Explanation
    - f0: nohit in aa length
    - mean: mean of $z$
    - variance: variance of $z$
    - m0: the recurrent rate of somatic mutations at a passenger site
    - m1: the recurrent rate of somatic mutations at a driver site
    - eta: $\eta$, the probability of any site being a driver in a gene
    - LRT_p-value: the p value of approximate likelihood ratio test for the null hypothesis that all sites are passengers ($\eta = 0$)
    - Q(z): $Q(z) = P(driver|z)$, the posterior probability of an amino acid site being a cancer-driver when the number ($z$) of somatic mutations is observed at this site

#### Result of ```site.txt```

- Title

| **Gene_id** | **Transcript_id** | **Gene_name** | **Chr** | **Protein_position** | **z** | **Q(z)** | **Protein_mutation** | **M(z)** | **Omega(z)** |
|-------------|-------------------|---------------|---------|----------------------|------|---------|----------------------|---------|--------------|

- Explanation
    - z: the count of somatic missense mutations at a site
    - M(z): the posterior mean of recurrent mutations at the site where $z$ mutations are observed
    - Omega(z): the site-specific selection pressure at the site where $z$ mutations are observed

#### Result of ```H_test.txt```

- Title

| **Gene_id** | **Transcript_id** | **Gene_name** | **Chr** | **CN/CS** | **1-H** | **H_test** | **Htest_p_value** |
|-------------|-------------------|---------------|---------|----------|---------|-----------|-------------------|

- Explanation
    - 1-H: $H$ is a relative measure of evolutionary rate variation among sites
    - H_test: the comparison of $1-H$ and $C_{N}/C_{S}$
    - Htest_p_value: the p value of chi-continency assert whether null hypothesis of $C_{N}/C_{S} = 1-H$ can be statistically rejected

#### Result of ```two_compo_cncs.txt```

- Title

| **Gene_id** | **Transcript_id** | **Gene_name** | **Chr** | **Omega** | **Omega_p** | **Omegadri** | **Omegapass** | **Omegapass_p** |
|-------------|-------------------|---------------|---------|-----------|-------------|------------|------------|--------------|

- Explanation
    - Omega: $C_{N}/C_{S}$ ratio of this gene
    - Omega_p: the p value of chi-continency assert whether null hypothesis of $C_{N}/C_{S} = 1$ can be statistically rejected
    - Omegadri: $\Omega_{dri}$, $C_{N}/C_{S}$ ratio of driver sites of the gene
    - Omegapass: $\Omega_{pass}$, $C_{N}/C_{S}$ ratio of passenger sites of the gene
    - Omegapass_p: the p value of chi-continency assert whether null hypothesis of $\Omega_{pass} = 1$ can be statistically rejected

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

```genetic_code.txt```                       symonymous or non-synonymous in codon mutation

```human_TCGA_match_GRCh37.fas```            fasta file of human genes (GRCh37)

```pancancer_mutation_rate.txt```	          Pancancer rate of Ref mutate to Alt	

```TCGA_match_GRCh37.cds```	                cds information of human genes

## Contact

Xun Gu		xgu@iastate.edu

Zhou Zhan	zhanzhou@zju.edu.cn