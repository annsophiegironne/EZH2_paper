Transposable elements annotations were obtained by running RepeatMasker (v4.1.5).
Annotation was performed by combining RepeatMasker with the HMMER tool (v3.3.1)
and the Dfam database (v3.7). RepeatMasker results in fa.out.gff format were converted
to a GTF file using a script provided by the Hammell lab (https://github.com/mhammell-laboratory).

This was the command used:
perl RepeatMasker -species [human/mouse] -gff -pa 32 -u -xm [fasta_genome_file] -dir [output_dir]