
# pdb2af v0.0.20

pdb2af uses SIFTS (https://www.ebi.ac.uk/pdbe/docs/sifts/) to download protein structures from the AlphaFold Protein Structure Database that overlap with PDB structures.

## Installation

pdb2af is installed using pip (https://packaging.python.org/en/latest/tutorials/installing-packages/). 

`pip install pdb2af`


## Running pdb2af

pdb2af reads in a text file with PDB IDs. You need to be connected to the internet.

`pdb2af -i pdb_ids.txt` 

Options:

`-i` text file with list of PDB IDs  

`-u` update the SIFTS mapping between PDB IDs and UniProt accession numbers


## Help

Contact Nick Fowler (njfowler.com) for support, queries or suggestions.











