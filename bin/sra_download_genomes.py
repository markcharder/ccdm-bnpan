#!/usr/bin/env python3

import os
import subprocess
from Bio import Entrez
from collections import Counter, defaultdict
from io import StringIO
from typing import List, Dict, Set

import argparse
import pandas as pd
import re
import time
import xml.etree.ElementTree as ET

parser = argparse.ArgumentParser(
    description = \
        "\tSelect a diverse set of RNA sequencing samples\n\
        from NCBI based on tissue diversity.",
    add_help = False,
    formatter_class = argparse.RawTextHelpFormatter
)

parser.add_argument(
    '--help',
    '-h',
    action = 'help',
    help = '\t\t\t\tShow this help message and exit.\n\n'
)

parser.add_argument(
    '--prefetch',
    '-p',
    action = 'store_true',
    default = False,
    help = \
        "\t\t\t\tDownload the selected SRA files using prefetch\n"
)

parser.add_argument(
    '--accessions',
    '-a',
    type = str,
    default = None,
    help = \
        "\t\t\t\tA file containing a list of SRA accessions to\n\
        \t\t\tdownload, one per line. If provided, these are\n\
        \t\t\tdownloaded regardless of other selection criteria.\n"
)

parser.add_argument(
    '--ncbi-organism',
    '-no',
    type = str,
    default = "Brassica napus",
    help = \
        "\t\t\t\tThe NCBI organism to use for the Entrez search.\n\
         \t\t\tIt will be used in this way, for example:\n\
         \t\t\t\tBrassica napus[Organism]"
)

parser.add_argument(
    '--ncbi-email',
    '-ne',
    type = str,
    help = \
        "\t\t\t\tEmail address to use with NCBI Entrez API.\n\
         \t\t\tThis is recommended by NCBI, who may contact\n\
         \t\t\tusers over request issues."
)

parser.add_argument(
    '--drop-unpaired',
    '-d',
    action = 'store_true',
    default = False,
    help = \
        "\t\t\t\tDrop unpaired reads from selection.\n\
         \t\t\tThis is recommended, unless paired RNA sequencing\n\
         \t\t\tdata are very limited for the organism of interest."
)

parser.add_argument(
    '--min-bases',
    '-mb',
    type = int,
    default = 5000000000,
    help = \
        "\t\t\t\tMinimum number of bases per sample.\n"
)

parser.add_argument(
    '--batch-size',
    '-b',
    type = int,
    default = 10000,
    help = "\t\t\t\tBatch size for Entrez fetch requests."
)

parser.add_argument(
    '--rate-limit',
    '-r',
    dest="ratelim",
    type = float,
    default = 0.5,
    help = "\t\t\t\tRate limit (seconds) between Entrez fetch requests."

)

parser.add_argument(
    '--output-dir',
    '-o',
    type = str,
    default = "./sra_downloads",
    help = "\t\t\t\tOutput directory for downloaded SRA files."
)

def query_builder(
    organism: str,
    strategy: str = "wgs",
    platform: str = "Illumina",
    drop_unpaired: bool = False
    ) -> str:

    query = f"{organism}[Organism] AND {strategy}[Strategy] AND {platform}[Platform]"

    # If we are limiting to paired libraries, add this to the query.
    if drop_unpaired:
        query += " AND PAIRED[LibraryLayout]"

    return query

def search_entrez(
    organism: str,
    email: str,
    strategy: str = "wgs",
    platform: str = "Illumina",
    drop_unpaired: bool = True,
    batch_size: int = 10000,
    ratelim: float = 0.5
    ) -> List[str]:

    # Build query from user info.
    query = query_builder(organism, strategy, platform, drop_unpaired)

    # Set user email.
    Entrez.email = email

    # Create an Entrez search handle (requres Entrez from Biopython).
    search_handle = \
        Entrez.esearch(
            db="sra",
            term=query,
            retmax=1,
            usehistory="y"
        )

    # Read the results from the search handle.
    entrez_results = Entrez.read(search_handle)
    search_handle.close()

    # This gives you the number of hits found.
    count = int(entrez_results['Count'])

    # This retrieves the search session ID stored on the ncbi server.
    webenv = entrez_results["WebEnv"]

    # This retrieves the query key for the search.
    # Each search gets a unique query key.
    query_key = entrez_results["QueryKey"]

    # Print out the total number of hits to stdout.
    print(f"Total results: {count}")

    # Initialise an empty list to hold Entrez results.
    entrez_results_all = []

    # Go from 0 to the number of hits, incrementing by batch size.
    for start in range(0, count, batch_size):
        
        # The end of the batch is the min of the total count
        # or the current start + batch size.
        end = min(start + batch_size, count)

        # Print progress to stdout.
        print(f"Fetching SRA records {start}-{end}...")

        # This handle fetches runinfo for the current batch.
        fetch_handle = \
            Entrez.efetch(
                db = "sra",
                rettype = "runinfo",
                retmode = "text",
                retstart = start,
                retmax = batch_size,
                webenv = webenv,
                query_key = query_key
            )

        # Read the fetched runinfo.
        data = fetch_handle.read()
        fetch_handle.close()

        # Add this to the entrez_results_all list.
        entrez_results_all.append(data)
        
        # Sleep for a second to limit the number of retrievals per second.
        time.sleep(ratelim)
    
    return entrez_results_all

def get_runinfo(
    entrez_results_all: List[str],
    min_bases: int = 5000000000
    ) -> pd.DataFrame:

    # Start the full data list by setting it to a string
    # matching the first batch.
    full_data = entrez_results_all[0]

    # Decode to utf-8.
    if isinstance(full_data, bytes):
        full_data = full_data.decode('utf-8')

    # For each additional batch, append to the string.
    for batch in entrez_results_all[1:]:

        # Decode each batch to utf-8, innit.
        if isinstance(batch, bytes):
            batch = batch.decode('utf-8')

        lines = batch.split('\n')
        full_data += '\n'.join(lines[1:])  # Skip header

    # Read the full data string into a pandas dataframe.
    runinfo_df = pd.read_csv(StringIO(full_data))
    runinfo_df = runinfo_df[runinfo_df['bases'] >= min_bases]
    # Return the pandas dataframe.
    return runinfo_df

def download_sra(
        accessions: Set[str] = set(),
        output_dir: str = "./sra_downloads"
        ):
    downloaded = []
    failed = []

    for acc in accessions:
        print(f"Processing {acc}...")
        
        try:
            # First prefetch the SRA file
            prefetch_cmd = f"prefetch {acc} -O {output_dir}"
            subprocess.run(prefetch_cmd, shell=True, check=True, capture_output=True)
            print(f"  Prefetched {acc}")
            
        except subprocess.CalledProcessError as e:
            print(f"  {acc} not found")
            failed.append(acc)
        
        time.sleep(1)
    
    print(f"\nDownload Summary: {len(downloaded)} successful, {len(failed)} failed")
    
    return downloaded

if __name__ == "__main__":

    options = parser.parse_args()

    if options.accessions:
        os.makedirs(options.output_dir, exist_ok=True)    

        with open(options.accessions, 'r') as f:
            accessions = set(f.read().splitlines())

        download_sra(
            accessions=accessions,
            output_dir=options.output_dir
        )

    else:
        os.makedirs(options.output_dir, exist_ok=True)    

        if not options.ncbi_email:
            raise ValueError(
                "An email address must be provided for NCBI Entrez API usage."
            )
        entrez_results = \
            search_entrez(
                organism=options.ncbi_organism,
                email=options.ncbi_email,
                drop_unpaired=options.drop_unpaired,
                batch_size=options.batch_size,
                ratelim=options.ratelim
            )

        run_info_df = get_runinfo(entrez_results, min_bases=options.min_bases)

        run_info_df.to_csv(
            os.path.join(
                options.output_dir,
                "sra_run_info.tsv"
            ),
            sep="\t",
            index=False
        )

        accessions = set(run_info_df['Run'].tolist())

        download_sra(
            accessions=accessions,
            output_dir=options.output_dir
        )

        accessions = run_info_df['Run'].tolist()

