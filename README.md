
## Mutational Signature Calculator
An easy to install and run mutational signature calculator.

This library estimates exposures, given a list of existing signatures.

## Installation
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt

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
python plot.py sample.png < sample.count
```

Then find the most likely combination of signatures:
```
python signature.py --signatures signatures.txt --counts sample.count > sample.signature
```
