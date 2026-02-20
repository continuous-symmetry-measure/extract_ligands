# program to extract the HETATMs from pdb files
# the progam uses hard coded filtering ligand-ids
import os
from Bio.PDB import PDBParser, PDBIO, Select
import argparse


ions = ['AU', 'BR', 'CA', 'CL', 'CO', 'EU', 'FE', 'K', 'KR', 'IOD',
        'MG', 'MN', 'NA', 'NI', 'RB', 'UNX', 'XE', 'Y1', 'YB', 'ZN']

small_mols = ['06B', '1F1', '1X4', '2FX', '3CH', '4DX', 'ACE', 'ACN', 'ACT',
              'ACY', 'AKG', 'ASN', 'AZI', 'B3P', 'BCT', 'BEN', 'BGC', 'BME',
              'BNZ', 'BTB', 'CAD', 'CHT', 'CMT', 'CO3', 'DMS', 'DTD', 'DIQ',
              'EDO', 'EOH', 'FG7', 'FMT', 'GAI', 'GOL', 'HED', 'IIL', 'AIB',
              '7GC', 'HV9', 'HGU', 'HPH', 'IX4', 'MAN', 'MNM', 'MES', 'MAH',
              'MLI', 'MOO', 'MPD', 'MRD', 'MTE', 'MTN', 'N4B', 'NAG', 'NET',
              'NH2', 'NO3', 'NTB', 'OXY', 'PEG', 'PG4', 'PGE', 'PO4', 'PXY',
              'QNC', 'SO4', 'TRS', 'IMD', 'K97', '27B']

non_standard_res = ['5GM', 'ABA', 'AIB', 'CCS', 'CSD', 'CME', 'CSO', 'DAL',
                    'DPP', 'GOA', 'IIL', 'MSE', 'NLE', 'OIL', 'OMI', 'DAL',
                    'MK8', 'SLZ', 'SMC', 'SCH', 'TIS', 'URE', 'YCM']

def is_het(residue):
    """ verifying it is hetatm and not water"""
    res = residue.id[0]
    return res != " " and res != "W"

def is_main_ligand(residue):
    """ Filter out small ligands ions and non-standard amino acids"""
    ligand_name = residue.get_resname()
    main_ligand = True
    if (ligand_name in ions) or \
            (ligand_name in small_mols) or \
            (ligand_name in non_standard_res):
        main_ligand = False
    return main_ligand


class ResidueSelect(Select):
    def __init__(self, chain, residue):
        self.chain = chain
        self.residue = residue

    def accept_chain(self, chain):
        return chain.id == self.chain.id

    def accept_residue(self, residue):
        """ Recognition of heteroatoms - Remove water molecules """
        return residue == self.residue and is_het(residue) and is_main_ligand(residue)

def extract_ligands(file_list, in_path, out_path):
    """ Extraction of the heteroatoms of .pdb files """

    for pfb_file in file_list:
        i = 1
        if pfb_file.endswith('.pdb') and not pfb_file.startswith("lig_"):
            pdb_code = pfb_file[:-4]
            # print(pdb_code)
            pdb = PDBParser(PERMISSIVE=1, QUIET=True).get_structure(pdb_code, in_path + pfb_file)
            io = PDBIO()
            io.set_structure(pdb)
            for model in pdb:
                for chain in model:
                    for residue in chain:
                        if not is_het(residue):
                            continue
                        if not is_main_ligand(residue):
                            continue
                        ligand_name = residue.get_resname()
                        # print(f"saving {chain} {residue}")
                        io.save(out_path + f"lig_{pdb_code}_{i}_{ligand_name}.pdb", ResidueSelect(chain, residue))
                        i += 1


# Main
def main():
    # Define the command line arguments
    argparser = (argparse.ArgumentParser
                 (prog="extract_ligands",
                  formatter_class=argparse.RawDescriptionHelpFormatter,
                  description="""The script extract HETATMs from pdb files and writes them as pdb files.
                     Ions, small ligands and non-standard residues are filtered out and not written """
                  ))
    argparser.add_argument('-d', '--directory', help="The directory of the files. Default: current directory")
    argparser.add_argument('-l', '--file-list',
                           help="[optional] A filename with the list of PDB files to extract the ligands from"),
    argparser.add_argument('-o', '--output',
                           help="The sub-dir path of the ligands. Default: 'main-ligands' sub-directory",
                           default=None)
    argparser.add_argument('--version', action='version', version='%(prog)s 0.1.0')
    args = argparser.parse_args()
    directory = args.directory
    file_list = args.file_list
    out_path = args.output

    # set the current directory as default
    if directory is None or directory == '.':
        directory = os.getcwd()
    if not directory.endswith("/"):
        directory = directory + "/"

    # set the output_dir default as main-ligands sub dir of directory
    if not out_path:
        out_path = directory + "/main-ligands/"

    if not out_path.endswith("/"):
        out_path = out_path + "/"

    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # prepare the list of files to work on
    if file_list is None:
        print("Info:  Processing all the pdb files in: ", directory)
        list_to_read = sorted(os.listdir(directory))
    else:
        file_list = directory + file_list
        with open(file_list, 'r') as file1:
            files = list(file1)
            list_to_read = list(map(lambda s: s.strip(), files))
        print("Info:  Processing all the pdb files in: ", file_list)

    try:
        extract_ligands(list_to_read, directory, out_path)
        print(f"Info:  Main ligands were extracted successfully ")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()


