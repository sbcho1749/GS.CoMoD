# GS.CoMoD

# disease genes and reference data
dm.genes: disease genes of type 2 diabetes from DisGeNET database

htn.genes: disease genes of hypertension from DisGeNET database

hg.ref: human genes


# functional.genesets.Rdata: a binary file having gobp, gomf, and reactome.gs

  gobp: a list containing gene ontology biological process gene sets

  gomf: a list containing gene ontology molecular function gene sets

  reactome.gs: a list containing reactome database gene sets


# mench.Rdata: a binary file having mench.data, mench.mapping, and mench.disease.genes

  mench.data: comorbidity scores from mench data

  mench.mapping: disease ID mapping result in mench data

  mench.disease.genes: a list of disease genes in mench data


# rubio.Rdata: a binary file having rubio.dat.V1, rubio.dat.V3, rubio.mapping, rubio.disease.genes.V1, and rubio.disease.genes.V3

rubio.dat.V1: comorbidity scores from rubio V1 disease gene data

rubio.dat.V3: comorbidity scores from rubio V3 disease gene data

rubio.mapping: disease ID mapping result in rubio data

rubio.disease.genes.V1: a list of disease genes from DisGeNET database V1

rubio.disease.genes.V3: a list of disease genes from current DisGeNET database. 



# GS.CoMoD.R: source code file for GS.CoMoD analysis (function + analysis examples)
