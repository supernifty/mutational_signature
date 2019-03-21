
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
python plot_components.py sample.png < sample.count
```

Find the most likely combination of signatures:
```
python decompose.py --signatures signatures.txt --counts sample.count > calculation.txt
```

Generate a signature plot for multiple samples
```
python plot_components.py --threshold 0.0 --show_signature --target out.png --order '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' --descriptions 'age related' 'APOBEC' 'double-strand break-repair failure' 'tobacco' '' 'MMRd' 'UV exposure' '' '' 'POLE mutations' 'alkylating agents' '' 'APOBEC' '' 'MMRd' '' '' '' '' 'MMRd' 'unknown aetiology' 'aristolochic acid' 'unknown aetiology' 'aflatoxin exposure' '' 'MMRd' '' 'unknown aetiology' 'tobacco' '' < example/sample.sigs
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
