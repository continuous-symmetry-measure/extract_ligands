# extract_ligands

Extract “main ligands” (HETATM residues) from PDB files and write each ligand into a separate PDB file.

This script scans each input `*.pdb` file, identifies HETATM residues, filters out waters/ions/small molecules/non-standard residues using **hard-coded residue name lists**, and writes each remaining ligand residue as its own file.

Output files are named:

`lig_<pdb-id>_<n>_<ligand-resname>.pdb`

Example:
`lig_1abc_1_INH.pdb`

Where n is ligands counter

---
## What the script considers a “ligand”

A residue is written out only if:

1. It is a **HETATM** residue (i.e., not a standard protein/nucleic residue)
2. It is **not water**
3. Its residue name (3-letter code) is **not** in any of the hard-coded exclusion lists:
   - `ions`
   - `small_mols`
   - `non_standard_res`

These lists are defined at the top of the script and should be edited to match your project.

---

## Requirements

- Python 3.11
- Biopython >= 1.81

## Installation
Install the tool from the source directory
```bash
pip install .
```

## Usage

Run on all PDB files in a directory:
```
extract_ligands.py -d /path/to/pdbs
```
**Optional:** provide a file list

If you only want to process specific files, create a text file with one PDB filename per line, e.g. files.txt:

1abc.pdb
2xyz.pdb
...

Then run:

```

extract_ligands.py -d /path/to/pdbs -l files.txt

``` 
**Output directory**

By default, ligand files are written to a subdirectory called:

<input-directory>/main-ligands/

You can override with -o:

```
extract_ligands.py -d /path/to/pdbs -o /path/to/output_ligands

```
**Output**

For each input pdb filr xxxx.pdb, the script writes one file per ligand residue found:

lig_xxxx_1_YYY.pdb
lig_xxxx_2_ZZZ.pdb
...

Where:

xxxx = input PDB id (filename without .pdb)

1, 2, ... = counter (per input file)

YYY = ligand residue name from the PDB (e.g., 017, INH, NAP)

Notes:

Files that already start with lig_ are skipped.

The counter resets for each input PDB file.

---

## Important limitations and caveats

1. **Hard-coded filtering**

The lists ions, small_mols, and non_standard_res are embedded in the script.

This is project-specific and should be reviewed/edited for new datasets.

2. **Ligands are extracted as residues**

Each output file contains the atoms of a single HETATM residue.

Multi-residue ligands (rare but possible) will be split across multiple files.

3. **Models**

The script iterates over all models, but the output naming does not include the model number.

If the same ligand appears in multiple models, outputs may overwrite or appear duplicated depending on PDB content.

4. **No chemical bonding information**

Output is a coordinate-only PDB. Bond orders/chemistry are not inferred or validated.

5. **Waters**

Waters are excluded using the residue id flag check (residue.id[0] == "W").


### How it works (brief)

Parses PDB with Bio.PDB.PDBParser(PERMISSIVE=1, QUIET=True)

Loops over: model -> chain -> residue

Checks:

is_het(residue) (not standard residue and not water)

is_main_ligand(residue) (not in exclusion lists)

Writes each accepted residue using Bio.PDB.PDBIO with a Select subclass to isolate the chain/residue.

## Version

1.0.0

## Troubleshooting

- **Nothing is written**

    - Confirm your ligands are actually HETATM residues in the PDB.

    - Check whether your ligand residue name is accidentally included in the exclusion lists.

- **Too many outputs / salts included**

    - Add residue codes to the exclusion lists at the top of the script.

- **Missing ligands** 
    - Verify the ligand is not considered a non-standard residue or small molecule by your current lists.

## License

This project is licensed under the MIT License.  
Copyright (c) 2026 Inbal Tuvi-Arad and Yaffa Shalit  
See 'LICENSE' file for details 

## Contact 
For questions or issues related to this script, please contact:

Inbal Tuvi-Arad  
inbaltu@openu.ac.il