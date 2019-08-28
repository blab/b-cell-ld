from Bio import SeqIO
from Bio import AlignIO
import numpy as np
import pandas as pd
import glob
from collections import defaultdict
import os

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--dataset', default="rabbit", help="dataset to compute on")

def get_majorminor(site):
    site = site[site.isin(['A', 'C', 'T', 'G'])]    # Remove the gaps and other invalid characters
    counts = site.value_counts()
    total = float(counts.sum())
    if total < 10 or len(counts) < 2:               # If not polymorphic, invalid.
        return None                                 # If < 10 characters left, invalid.
    counts = counts[(counts/total) >= 0.10]         # Only want sites that have minor alleles of at least 10% frequency
    if len(counts) == 2:
        return tuple(counts.index)                  # If that leaves us with two alleles (major and minor), return just the alleles
    else:
        return None                                 # Otherwise, the site is not biallelic; return None.

def r_squared(haplotypes):
    hapObs = pd.Series(haplotypes).value_counts() # Tally observed counts of haplotypes and alleles
    if len(haplotypes) < 10: # Check, yet again, that we have enough data for this site
        return np.nan
    else:
        pass

    allele_i_obs = pd.Series([h[0] for h in haplotypes]).value_counts()
    allele_j_obs = pd.Series([h[1] for h in haplotypes]).value_counts()

    # Check again that we have two alleles per site
    if len(allele_i_obs) != 2 or len(allele_j_obs) != 2:
        return None
    else:
        pass

    hapExp = {}
    N = float(hapObs.sum()) # Total N of observed haplotypes

    for i in allele_i_obs.index.values:
        for j in allele_j_obs.index.values:
            # Exp counts = row (allele0|site0) total*column (allele1|site1) total / table (allele) total
            # If this haplotype is never observed, use 0
            hapExp[i+j] = (float(allele_i_obs[i])*float(allele_j_obs[j]))/N

    chisq = 0.0
    for hap in hapExp.keys(): # Calculate chi-squared
                              # chisq = sum( (E-O)**2/E for each possible haplotype )
        if hap in hapObs:     # deal with unobserved haplotypes
            chisq += ((float(hapObs[hap]) - hapExp[hap])**2)/hapExp[hap]
        else:
            chisq += ((0.0 - hapExp[hap])**2)/hapExp[hap]

    # r^2 = chisq / (Nsamples * df) ; because all sites biallelic, df == (2-1)(2-1) == 1
    r_sq = chisq / (N)
    return r_sq

def get_rsq(file, isotype):

    df = pd.DataFrame(columns=['isotype', 'gene', 'distance','rsq'])

    gene = file.split('/')[-1].split('.')[0]
    print("gene", gene)

    align = AlignIO.read(file, 'fasta')
    align = align[1:]                                                                               # Dropping first sequence, which is germline
    align = np.array([list(rec) for rec in align], np.character, order="F")                         # Read in alignment, make a dataframe
    nsites = len(align.T)
    nseqs = len(align)
    align = pd.DataFrame(align, columns=range(nsites), index=range(nseqs))                          # Use pandas indexing to ke

    majorminor = align.apply(lambda x: get_majorminor(x))                                           # Get the (major,minor) alleles for each site if possible
    majorminor.dropna(inplace=True)                                                                 # Drop the non-biallelic sites
    validsites = majorminor.index.values                                                            # Valid sites are the biallelic ones

    for x in range(len(validsites)):                                                                # Make upper-triangle comparisons
        i = validsites[x]
        column_i = align.loc[:,i]                                                                   # Pull the alignment column
        alleles_i = majorminor.loc[i]                                                               # And the (major,minor) alleles we found earlier

        for y in range(x+1, len(validsites)):
            j = validsites[y]
            column_j = align[j]
            alleles_j = majorminor.loc[j]                                                           # Make a list of all the observed haplotypes that consist of the major or minor allele for each site

            haplotypes = [ (h[0]+h[1]).upper() for h in zip(column_i, column_j) if h[0] in alleles_i and h[1] in alleles_j]
            print(haplotypes)

            if len(haplotypes) < 10:                                                                # Do we still have enough data?
                continue
            else:
                dist = np.abs(i-j)
                r_sq = r_squared(haplotypes)
                print(r_sq)
                if r_sq is not None:
                    df.at[len(df)] = [isotype, gene, dist, r_sq]
    return df

if __name__=="__main__":
    params = parser.parse_args()

    results = pd.DataFrame(columns=['isotype', 'gene', 'distance','rsq'])

    isotypes = ['Ig']
    for isotype in isotypes:
        files = [file for file in glob.glob('fastas/' + params.dataset + '/*')]
        for file in files:
            print(file)
            df = get_rsq(file, isotype)
            print(df)
            if len(df) > 1:
                results = results.append(df)

    results.to_csv("rsq_tables/rsq_" + params.dataset + ".tsv", sep='\t', index=False)
