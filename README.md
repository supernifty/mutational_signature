
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
cd decompose
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
python plot.py sample.png < sample.count
```

Then find the most likely combination of signatures:
```
python decompose.py --signatures signatures.txt --counts sample.count > sample.signature
```
