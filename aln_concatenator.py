''' Merges alignments based on template alignments and protein family '''

import GPCR_Project_func as GPCR_Proj

def build_templates(alns):
    ''' Builds all the template files by merging the individual pairwise
    structural alignments one by one for each protein 
    
    Arguments
    - aln: Refers to [directory of alignment files|single file with all
      pairwise alignments]'''

    # Some way to get each alignment with the same protein as a list of lists,
    # such that [[A1, A2...An], [B1, B2...Bn]...]
    
    merged_aln_list = []
    for alns in ordered_alns:
        merged_aln = alns[0]
        for i, aln in enumerate(alns):
            if i != 0:
                merged_aln = GPCR_Proj.aln_merger(merged_aln, aln)

        merged_aln_list.append(merged_aln)

    return merged_aln_list

def merge_family(template, family):
    ''' Applies the alignment transform given a template to the family of
    proteins the template belongs to

    Arguments
    - template: Template sequence alignment object
    - family: list of sequence objects for the proteins in that family'''
    
    merged_family_list = []
    for aln in family:
        merged_family = GPCR_Proj.aln_merger(template, aln)
        merge_family_list.append(merged_family)

    return merged_family_list

def main():
    ''' Main function which takes in the necessary files, builds the templates
    and then merges all family files together into one alignment '''

    # Build the template

    # Get each individual family, merge with the corresponding template and add
    # to the main list

    # Save main list as a fasta file with all alignments merged
