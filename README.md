# ggCaller_manuscript
Repository of scripts used in ggCaller manuscript.

## Generating a simulated pangenome

### Population Simulation
To generate a simulated population for pangenome analysis, run:

```python simulate_full_pangenome.py --gff data/simulated_pangenome/SP_ATCC700669.gff3 --nisolates 100 --n_sim_genes 1000 --pop_size 10e-6 --out sim_pangenome --mutation_rate 1e-14 --gain_rate 1e-12 --loss_rate 1e-12```

This will generate a list of fasta files in the directory ```sim_pangenome```.

### Adding fragmentation
To fragment assemblies, first copy ```sim_pangenome```. 
```cp sim_pangenome sim_pangenome_fragmented```

Then choose a directory ```--ref_dir``` to get empirical read length distributions from (should contain ```.fasta``` files).
Finally, run:
```
python fragment_fasta.py --ref_dir reference_dir --sample_dir sim_pangenome_fragmented --min_length 10
```

### Adding contaminants
To fragment assemblies, first copy ```sim_pangenome```. 
```cp sim_pangenome sim_pangenome_contaminated```

Then choose a infile ```--infile``` to generate fragments from, and --outfile to append to (this is done IN PLACE).
Finally, run:
```
python insert_random_genome_fragments.py --infile contaminant.fasta --outfile sim_pangenome_fragmented/genome1.fasta --frag_size 10000
```

### Simulating assemblies
Now simulate paired-end reads from your assemblies with ART:

```art_illumina -ss HS25 -sam -i sim_pangenome/genome1.fasta -p -l 150 -f 20 -m 200 -s 10 -o genome1_reads```

Assemble reads with SPADES:

```spades.py --phred-offset 33 --isolate --threads 1 -1 genome1_reads1.fq -2 genome1_reads2.fq -o genome1_assembly```

### Gene identification and pangenome analysis
Now run Prokka on SPADES assemblies, using the gene annotation from ```simulate_full_pangenome.py```:
```prokka --cpus 8 --proteins sim_pangenome_prokka_DB.fa --force --outdir prokka_genome1 --notrna --norrna --prefix genome1 genome1_assembly/scaffolds.fasta```

Run Panaroo on Prokka gene-calls (using moderate clean-mode):
```panaroo -i prokka_gffs/*.gff -o panaroo_moderate_out --clean-mode moderate -t 8```

Run Roary on Prokka gene-calls:
```roary -p 8 -f roary prokka_gffs/*.gff```

Run PEPPAN on Prokka gene-calls (use PEPPAN_parser to generate gene presence/absence matrix):
```
PEPPAN -t 8 -p PEPPAN_out prokka_gffs/*.gff
PEPPAN_parser -g PEPPAN_out.PEPPAN.gff -s PEPPAN_out
```

Run ggCaller on all assemblies using the gene annotation from ```simulate_full_pangenome.py```:
```ggcaller ggcaller --refs sim_pangenome_list.txt --threads 8 --annotation sensitive --diamonddb sim_pangenome_prokka_DB.fa --out ggCaller_out --clean-mode moderate```

### Analysing simulated pangenome results
Copy gene presence/absence matrices from respective tools (```gene_presence_absence_roary.csv``` for Roary, Panaroo and ggCaller, ```*.PEPPAN.gene_content.Rtab``` for PEPPAN).
Also copy Prokka gff files for each simulation.

Use ```scripts/compare_simulated_gene_pa.Rmd``` to generate prokka mapping files from gffs, and to compare different pangenome analysis tools. Data is available in ```data/simulated_pangenome```.

## Comparing real bacterial pangenomes
Gene presence/absence matrices for M. tuberculosis, S. pneumoniae and E. coli used in the ggCaller paper are available in ```data/real_pangenome```

Prokka, Roary, PEPPAN, Panaroo and ggcaller were run using the parameters in the previous section "Gene identification and pangenome analysis".

## Contig break analysis
