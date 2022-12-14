# ggCaller_manuscript
Repository of scripts used in [ggCaller](https://github.com/samhorsfield96/ggCaller) manuscript.

## Generating a simulated pangenome

Data is available in ```data/simulation_pangenome```.

### Population Simulation
To generate a simulated population for pangenome analysis, run:

```python scripts/simulate_full_pangenome.py --gff data/simulated_pangenome/SP_ATCC700669.gff3 --nisolates 100 --n_sim_genes 1000 --pop_size 10e-6 --out sim_pangenome --mutation_rate 1e-14 --gain_rate 1e-12 --loss_rate 1e-12```

This will generate a list of fasta files in the directory ```sim_pangenome```.

### Adding fragmentation
To fragment assemblies, first copy ```sim_pangenome```. 

```cp sim_pangenome sim_pangenome_fragmented```

Then choose a directory ```--ref_dir``` to get empirical read length distributions from (should contain ```.fasta``` files).
Finally, run:

```
python scripts/fragment_fasta.py --ref_dir reference_dir --sample_dir sim_pangenome_fragmented --min_length 10
```

### Adding contaminants
To fragment assemblies, first copy ```sim_pangenome```. 

```cp sim_pangenome sim_pangenome_contaminated```

Then choose a infile ```--infile``` to generate fragments from, and --outfile to append to (this is done IN PLACE).
Finally, run:

```
python scripts/insert_random_genome_fragments.py --infile contaminant.fasta --outfile sim_pangenome_fragmented/genome1.fasta --frag_size 10000
```

### Simulating assemblies
Now simulate paired-end reads from your assemblies with [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm):

```art_illumina -ss HS25 -sam -i sim_pangenome/genome1.fasta -p -l 150 -f 20 -m 200 -s 10 -o genome1_reads```

Assemble reads with [SPADES](https://github.com/ablab/spades):

```spades.py --phred-offset 33 --isolate --threads 1 -1 genome1_reads1.fq -2 genome1_reads2.fq -o genome1_assembly```

### Gene identification and pangenome analysis
Now run [Prokka](https://github.com/tseemann/prokka) on SPADES assemblies, using the gene annotation from ```simulate_full_pangenome.py```:

```prokka --cpus 8 --proteins sim_pangenome_prokka_DB.fa --force --outdir prokka_genome1 --notrna --norrna --prefix genome1 genome1_assembly/scaffolds.fasta```

Run [Panaroo](https://github.com/gtonkinhill/panaroo) on Prokka gene-calls (using moderate clean-mode):

```panaroo -i prokka_gffs/*.gff -o panaroo_moderate_out --clean-mode moderate -t 8```

Run [Roary](https://github.com/sanger-pathogens/Roary) on Prokka gene-calls:

```roary -p 8 -f roary prokka_gffs/*.gff```

Run [PEPPAN](https://github.com/zheminzhou/PEPPAN) on Prokka gene-calls (use PEPPAN_parser to generate gene presence/absence matrix):

```
PEPPAN -t 8 -p PEPPAN_out prokka_gffs/*.gff
PEPPAN_parser -g PEPPAN_out.PEPPAN.gff -s PEPPAN_out
```

Run ggCaller on all assemblies using the gene annotation from ```simulate_full_pangenome.py```:

```ggcaller --refs sim_pangenome_list.txt --threads 8 --annotation sensitive --diamonddb sim_pangenome_prokka_DB.fa --out ggCaller_out --clean-mode moderate```

### Analysing simulated pangenome results

Copy gene presence/absence matrices from respective tools (```gene_presence_absence_roary.csv``` for Roary, Panaroo and ggCaller, ```*.PEPPAN.gene_content.Rtab``` for PEPPAN).
Also copy Prokka gff files for each simulation.

Use ```scripts/compare_simulated_gene_pa.Rmd``` to generate prokka mapping files from gffs, and to compare different pangenome analysis tools. Data is available in ```data/simulated_pangenome```.

This will summary tables detailing false positives, false negatives and correctly-called COGs containing set numbers of errors.

## Comparing real bacterial pangenomes

Gene presence/absence matrices for M. tuberculosis, S. pneumoniae and E. coli used in the ggCaller paper are available in ```data/real_pangenome```

For a chosen dataset, run Prokka, Roary, PEPPAN, Panaroo and ggcaller using the parameters in [Gene identification and pangenome analysis](###Gene identification and pangenome analysis).

Respective gene annotations were provided from:
- [M. tuberculosis](https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3)
- [S. pneumoniae](https://datadryad.org/stash/downloads/file_stream/67467)
- [E. coli](https://microbiology.figshare.com/ndownloader/files/26133194)

## Contig break analysis

Fragment your chosen sequence:

```
python fragment_at_gene.py --CDS data/contig_break/fasta/CDS/CR931662_Streptococcus_pneumoniae_strain_34359_serotype_14_CDS.fa --infile data/contig_break/fasta/CDS/original/CR931662_Streptococcus_pneumoniae_strain_34359_serotype_14.fa
```

Call genes using [GeneMarkS-2](http://exon.gatech.edu/genemark/genemarks2.cgi), or using Prokka. Then use ggCaller or Panaroo as detailed [Gene identification and pangenome analysis](###Gene identification and pangenome analysis).

Analyse gene recall using ```combined_DNA_CDS.fasta``` from Panaroo and ```gene_calls.fasta``` from ggCaller

```
python scripts/gene_recall.py --query data/contig_break/ggc/ggc_group3_fragmented/gene_calls.ffn --seq data/contig_break/ggc/fasta/all_seqs.fasta --CDS data/contig_break/ggc/fasta/all_CDS.fasta --exact
```

This will print precision and recall statistics to the console, and generate ```length_proportions.txt``` which describes the proportions of real genes covered by the respective ORF call.

## Gene end comparison

Fasta, alignment and summary files from previous analysis are available in ```data/gene_end_comparison```.

For a chosen dataset, run Prokka, Panaroo and ggCaller using the parameters in [Gene identification and pangenome analysis](###Gene identification and pangenome analysis).

For panaroo, generate a prokka mapping file out a given gene whilst in a directory containing all prokka gffs:

```grep "<target>" *.gff > prokkamap.txt ```

Pull out all genes of interest from ggCaller, reference (example [here](https://datadryad.org/stash/downloads/file_stream/67467)) or panaroo output directories:

```
python scripts/parse_genes.py --fastafile ggCaller_out/gene_calls.faa --targets CLS02009,CLS02800,CLS00182 --outpref aggc_pspA
python scripts/parse_genes.py --fastafile reference/proteins.faa --targets CLS02009,CLS02800,CLS00182 --outpref ref_pspA
python scripts/parse_genes.py --prokkamap prokkamap.txt --panaroodata panaroo_out/gene_data.csv --outpref panaroo_pspA
```

Append a reference protein sequence to the start of each file. Run mafft on files to generate MSA for all genes and tools

```
cat pspA.faa ggc_pspA.faa > ggc_pspA_cat.faa
mafft ggc_pspA_cat.faa > GGC_pbp1a.aln
```

Analyse alignments of gene ends. All ```.aln``` files should be in the same directory

```
python scripts/gene_end_comparison.py --indir all_alignments --outpref results
```

This will generate a series of summary graphs describing gene start position comparisons and percent identity based on MSAs.

## PanGenome wide association study (PGWAS)

Generate a presence/absence matrix for each isolate in the study. Examples are available in ```data/PGWAS/tetracycline/tetracycline_resistance.txt```
and ```data/PGWAS/erythromycin/erythromycin_resistance.txt```

For this study, S. pneumoniae gene annotations were supplied from [here](https://datadryad.org/stash/downloads/file_stream/67467).

Generate a ```.nwk``` tree using ggCaller (use strict cut-offs for graph cleaning and tree generation) and root at midpoint.

```
ggcaller --refs genome_list.txt --out pyseer_results --clean-mode strict --ignore-pseduogenes --alignment core --aligner def --annotation sensitive --diamonddb CDS_protein_sequences.dmnd --save --core-threshold 1
python scripts/midpoint_root.py --infile pyseer_results/core_tree.nwk --outfile core_tree_NJ_midpoint.nwk
```

Convert the nwk tree to distance matrix using pyseer

```
python pyseer/phylogeny_distance.py --lmm core_tree_NJ_midpoint.nwk > phylogeny_K.tsv
```

Generate unitigs using unitig-caller

```
unitig-caller --call --refs refs.txt --out PGWAS_unitigs
```

Calculate unitig associations with phenotype and determine adjusted signficance cut-off threshold.

```
pyseer --lmm --phenotypes erythromycin_resistance.txt --kmers PGWAS_unitigs.pyseer.gz --similarity phylogeny_K.tsv --output-patterns unitigs_patterns.txt --cpu 8 > unitigs_hits.txt
python pyseer/count_patterns.py unitigs_patterns.txt
cat <(head -1 tet_unitigs.txt) <(awk '$4<THRESHOLD {print $0}' unitigs_hits.txt) > significant_unitigs.txt
```

Map unitigs to a single reference with phandango and annotate

```
phandango_mapper significant_kmers.txt ref.fa phandango.plot
annotate_hits_pyseer significant_unitigs.txt references.txt annotated_unitigs.txt
```

Generate query in correct format, query singificant hits in ggCaller graph

```
awk -F '\t' 'NR>1 && NF=1' significant_unitigs.txt > significant_unitigs_query.txt
ggcaller --graph genome_list.gfa --colours genome_list.bfg_colors --threads 16 --out pyseer_query --data pyseer_results/ggc_data --query significant_unitigs_query.txt --query-id 1.0
```

Analyse hits. Annotation files should match that used [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5708525/#SD1)

```
python scripts/count_annotations.py --fasta pyseer_query/matched_queries.fasta --outpref annotations --annotations annotations_file.xlsx
```

This will generate a series of graphs describing the gene annotations covered by significant unitigs.

## Computational benchmarking

Data is available in ```data/computational_benchmarking```.

Tools were run as detailed in [Gene identification and pangenome analysis](###Gene identification and pangenome analysis).