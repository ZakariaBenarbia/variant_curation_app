-------------------------
COSMIC Cancer Gene Census
-------------------------
 A list of all cancer census genes from the current release in a comma separated table. The census table is exported from https://cancer.sanger.ac.uk/census and the format is the same. [ Cosmic_CancerGeneCensus_v103_GRCh38.tsv.gz ]

  File Description

[column number:label] Heading                                           Description
--------------------------------------------------------------------------------------------------------
[1:A]                GENE_SYMBOL                  The gene name for which the data has been curated in COSMIC. In most cases this is the accepted HGNC identifier.
[2:B]                NAME                         Gene descriptive name.
[3:C]                COSMIC_GENE_ID               A unique COSMIC gene identifier (COSG) is used to identify a gene within the file. This identifier can be used to retrieve additional Gene information from the Cosmic_Genes file. Hallmarks of cancer annotations exist for a proportion of Tier 1 genes, these can be accessed with this identifier from the Cosmic_CancerGeneCensusHallmarksOfCancer file.
[4:D]                CHROMOSOME                   The chromosome location of a given mutation census (1-22, X, Y or MT).
[5:E]                GENOME_START                 The start coordinate of a given mutation census.
[6:F]                GENOME_STOP                  The end coordinate of a given mutation census.
[7:G]                CHR_BAND                     Chromosome (1-22, X, Y or MT), arm (p or q) and cytogenic band.
[8:H]                SOMATIC                      Somatic mutations have been detected (y/n).
[9:I]                GERMLINE                     Germline mutations have been detected (y/n).
[10:J]               TUMOUR_TYPES_SOMATIC         Somatic mutations in the gene are associated with the following diseases (see abbreviations tab for details: https://cancer.sanger.ac.uk/cosmic/help/census#abbrev).
[11:K]               TUMOUR_TYPES_GERMLINE        Germline mutations in the gene are associated with the following diseases (see abbreviations tab for details: https://cancer.sanger.ac.uk/cosmic/help/census#abbrev).
[12:L]               CANCER_SYNDROME              Syndrome associated with germline mutation.
[13:M]               TISSUE_TYPE                  Type of tissue, see abbreviations tab for details: https://cancer.sanger.ac.uk/cosmic/help/census#abbrev.
[14:N]               MOLECULAR_GENETICS           See abbreviations tab for details: https://cancer.sanger.ac.uk/cosmic/help/census#abbrev.
[15:O]               ROLE_IN_CANCER               Role in Cancer: oncogene: hyperactivity of the gene drives the transformation; TSG: loss of gene function drives the transformation. Some genes can play either of these roles depending on cancer type. Fusion: the gene is known to be involved in oncogenic fusions.
[16:P]               MUTATION_TYPES               Types of mutation: See abbreviations tab for details: https://cancer.sanger.ac.uk/cosmic/help/census#abbrev.
[17:Q]               TRANSLOCATION_PARTNER        Gene symbol of fusion partner.
[18:R]               OTHER_GERMLINE_MUT           Other germline mutations not implicated in cancer.
[19:S]               OTHER_SYNDROME               Other non-cancerous syndrome.
[20:T]               TIER                         Indicates to which tier of the Cancer Gene Census the gene belongs (1/2)
[21:U]               SYNONYMS                     Gene alternative names.
