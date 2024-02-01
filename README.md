
## Mutational Signature Calculator
An easy to install and run mutational signature calculator.

This library estimates exposures, given a list of existing signatures.

## Installation

On spartan:
```
module load Python/3.6.4-intel-2017.u2
module load cURL/7.60.0-spartan_gcc-6.2.0
```

```
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Inputs
* genome.fa
* sample.vcf
* signatures.txt

## Usage
First calculate counts:
```
python count.py --genome genome.fa --vcf sample.vcf > sample.count
```

Make a signature plot:
```
python plot_counts.py --target sample.png < sample.count
```

Find the most likely combination of signatures:
```
python decompose.py --signatures signatures.txt --counts sample.count > calculation.txt
```

Generate a signature plot for multiple samples
```
python plot_components.py --threshold 0.0 --show_signature --target out.png --order '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' --descriptions 'age related' 'APOBEC' 'double-strand break-repair failure' 'tobacco' '' 'MMRd' 'UV exposure' '' '' 'POLE mutations' 'alkylating agents' '' 'APOBEC' '' 'MMRd' '' '' '' '' 'MMRd' 'unknown aetiology' 'aristolochic acid' 'unknown aetiology' 'aflatoxin exposure' '' 'MMRd' '' 'unknown aetiology' 'tobacco' '' < example/sample.sigs
```

Convert a downloaded signature file from COSMIC to format for this software
```
python mutational_signature/convert.py --conversion sbs32 < ~/sigProfiler_SBS_TCGA_WES_ColoRect-AdenoCa_local_signatures.csv > data/wes-crc.csv
```

Compare sets of signatures for similarity:
```
python mutational_signature/compare_signatures.py --signatures data/signatures_cosmic_v3_sbs.txt --signatures2 data/wes-crc.csv
```

Count available contexts for potential exposure adjustment:
```
python mutational_signature/count_contexts.py --genome genome.fa --bed regions.wgs.bed > opportunity.sbs.wgs.tsv
```

Adjust counts for different captures based on opportunity
```
python mutational_signature/adjust_counts.py --adjust_from opportunity.sbs.cog.tsv --adjust_to opportunity.sbs.exome.tsv --verbose < sample.count
```

Assign artefact probability
```
```

Annotate VCF with context and extended sequence
```
python mutational_signature/annotate_context.py --genome genome.fa --vcf example/sample.vcf.gz --sequence 6 > x.vcf
python mutational_signature/annotate_context_summary.py --vcfs x.vcf --sequence 6 --lengths 3 2 1 > y.tsv
```

### Generate de novo signatures
```
python mutational_signature/generate.py \
  --counts *.counts \
  --components 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 \
  --compare_to data/signatures_cosmic_v3_sbs.txt \
  --compare_to_matches 1 \
  --min_mutation_proportion 0.01 \
  --bootstraps 100 \
  --bootstrap_plot 'bootstraps.png' \
  --component_plot 'component.png'
```

### Backward selection based on cosine error:
```
python mutational_signature/refine_signatures.py \
  --signatures data/signatures_cosmic_v3.1_sbs.txt \
  --verbose < counts_file
```

### Signature distribution panel
```
python mutational_signature/histograms.py --exclude_below 0.05 --exposures *.exposures --target hists.sbs.png --cols 6
```

### Simulation
```
python mutational_signature/simulate.py --verbose --n 1 --sigdefs data/signatures_cosmic_v3.2_sbs.txt --noise 10,5 --exposures SBS1,100,40 --injection SBS88,0.05
```

## Functionality
* combine_counts: for multiple count files, combine into a single tsv
* combine_signatures: for multiple signature files, combine into a single tsv
* compare_signatures: measures cosine similarity between input signatures
* count_indels: counts indels in repeat regions using a an annotated bed file
* count_maf: counts SNVs in context from a maf file
* count: counts SNVs in context from a VCF file
* decompose: calculate signature profile for a single sample
* generate: build a new set of base signatures (incomplete)
* plot_components: plot a signature profile
