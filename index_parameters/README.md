## Annotation versions

Homo_sapiens
    Fasta: GRCh38
    GTF: Gencode40, "basic" annotations
Mus_musculus
    Fasta: GRCm39
    GTF: GencodeM29, "basic" annotations
Drosophila_melanogaster
    Fasta: BDGP6.28
    GTF: Ensembl99, "basic" annotations


For mixed genomes, chromosomes names in GTF and FASTA files were modified to include a species-specific suffix. E.g, instead of:

Homo_sapiens.GRCh38.Gencode40.gtf
    ##description: evidence-based annotation of the human genome (GRCh38), version 40 (Ensembl 106)
    ##provider: GENCODE
    ##contact: gencode-help@ebi.ac.uk
    ##format: gtf
    ##date: 2022-01-20
    chr1    HAVANA  gene    11869   14409   .       +       .       gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102"; havana_gene "OTTHUMG00000000961.2"; transcript_id "";

we have

Homo_sapiens.GRCh38.Gencode40_sed.gtf
    ##description: evidence-based annotation of the human genome (GRCh38), version 40 (Ensembl 106)
    ##provider: GENCODE
    ##contact: gencode-help@ebi.ac.uk
    ##format: gtf
    ##date: 2022-01-20
    chr1_HS HAVANA  gene    11869   14409   .       +       .       gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102"; havana_gene "OTTHUMG00000000961.2"; transcript_id "";

so that reads specific to each specie could be splitted in downstream analysis. 
We have added a code that does it: src/01_processing/03_codes/create_sed.jl

List of acronyms
    Homo sapiens: HS
    Mus musculus: MM
    Drosophila melanogaster: DM