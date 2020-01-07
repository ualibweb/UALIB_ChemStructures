# UALIB_ChemStructures
Chemical substances from The University of Alabama Dissertations and Theses

## What is this?
This repository contains the original indexed chemical substances (non-standardized) 
from The University of Alabama Dissertations and Theses (hereafter, theses). Each 
thesis data folder contains the name, SMILES, and InChI of all synthesized chemical
structures within the thesis along with a permalink to the thesis fulltext or 
catalog link (if not yet available online), Moreover, an SDfile containing the name, 
SMILES, InChI, citation, permalink, and local structure registry ID is included. 
The SDfile is what we submit directly to PubChem. PubChem is likely where you want to download our 
data from as PubChem handles the standardization of the chemical structures:
[The University of Alabama Libraries PubChem Data Source](https://pubchem.ncbi.nlm.nih.gov/source/15645).
However, this GitHub repository is useful for those seeking to download our non-standardized data
and understand our workflow.

## How are we indexing the chemical structures?

1. This is a manually curated effort, no machine extraction. 
We are redrawing and encoding the structures ourselves.

2. We include only structures that can be represented with SMILES. 
Essentially this includes small molecule organic chemistry with some
limited organometallic chemical substances.

3. Chemical substances indexed must have synthetic experimental characterization
details such as NMR, mass spec, elemental analysis, and melting point. 
No judgment is made on the the accuracy of the reported syntheses. 
If the authors claimed they synthesized the substance in an experimental section and
includes characterization data, we indexed it. 

4. The substance name given in the dissertation is used, not a systematic name. 
If no name is given, we generated an IUPAC name.

5. Duplicates are assigned the same local Registry ID (i.e., UALIB-###). We check 
for duplicates via the InChIKey, then update the original PubChem Substance record 
with the additional thesis reference. 

A complete list of all structures indexed is available in the file 
UALIB_Chemical_Structures_REGID.txt

## What is your workflow? 

We're still working on creating an optimal workflow that is highly transparent,
open and reproducible. Here is an overview of our current strategy:

1. Draw chemical structures in [PubChem Sketcher V2.4](https://pubchem.ncbi.nlm.nih.gov/edit3/index.html)
2. Export SMILES
3. Create a tab delimited file with an index number for each compound (i.e., UALIB-###), 
the name, citation, and permalink to the thesis. 
4. Import the tabbed text file into an RDKit Pandas dataframe to calculate the 
InChIKeys, RDKit SMILES, and generate the SDfile.
5. Update the UALIB_Chemical_Structures_REGID.txt file and submit the SDfile to PubChem.
6. Create a record on our [Institutional Repository](https://ir.ua.edu/) with the SDfile 
and InChIkeys associated with the dissertation/thesis. 

## Notes on Copyright and Reuse

Disclaimer: Not legal advice, just our own personal (non-lawyer) thoughtful 
interpretation.

The purpose of The University of Alabama Dissertation and Thesis indexing project 
is to allow greater discovery and credit of the original authors' theses, 
not to claim any ownership of the written thesis content. The thesis authors hold the 
copyright to their own thesis.

We have only extracted scientific facts (i.e., the chemical structures and bibliographic
information) from the theses. Such facts and bibliographic data are not subject 
to U.S. copyright protection: [Compendium of U.S. Copyright Office Practices](https://www.copyright.gov/comp3/).
See specifically section 313.3(A), where examples are listed that are excluded 
from copyright protection, one of which mentions chemical structures:

..."DNA sequences and other genetic, biological, or chemical substances or 
compounds, regardless of whether they are man-made or produced by nature..."

In an effort to fully comply with copyright, Fair Use, and standard scholarly 
practice of reusing a thesis, we did not use any automated chemical structure 
extraction software. Instead, we produced our own chemical structures by using 
the theses as a reference and encoding the structures manually ourselves.

With that said, we are grateful to the thesis authors for their contributions
to science and have endeavored to credit their work respectfully and wherever 
possible by including a citation reference and permalink. The citation references and
permalinks are included on all shared chemical structure data including within
this Git Repository, The University of Alabama Institutional Repository,
and PubChem (Substance Pages).

Our intention is for you to reuse the chemical structures 
however you like. All chemical structure data in this repository is licensed 
with [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/). Please give the 
original authors of the theses credit by citing their work and following 
standard scholarly practice for reuse of the scientific literature, 
particularly if the chemical structure has led you to useful content within the 
thesis (as noted above, each thesis Author holds their own copyright). If you are
reusing a large corpus of the structures as a dataset, then it is appreciated 
if you cite our work as we put the effort into the indexing, reproduction of 
the structures, and compilation. Lastly, any code in this repository is 
licensed under the BSD-2 license.

## References

Much of our inspiration for this project came from the following similar project:

Andrews, D. M.; Broad, L. M.; Edwards, P. J.; Fox, D. N. A.; Gallagher, T.;
Garland, S. L.; Kidd, R.; Sweeney, J. B. The Creation and Characterisation of 
a National Compound Collection: The Royal Society of Chemistry Pilot. Chem. Sci. 2016,
7 (6), 3869–3878. [DOI:10.1039/C6SC00264A](https://doi.org/10.1039/C6SC00264A)





