import pandas as pd
import numpy as np
import time
from math import ceil
import multiprocessing as mp
import argparse
import subprocess
import sys


parser = argparse.ArgumentParser(description='Transforms a bed file to a numpy array')
parser.add_argument('-n', '--n_cores', type=int, help='Specifies the number of cores used', required=True)
parser.add_argument('-i', '--input', type=str, help='Input .bed file (can be gzip compressed)', required=True)
parser.add_argument('-o', '--output', type=str, help='Output name (.txt or .txt.gz)', required=True)
parser.add_argument('-v', '--vcf', type=str, help='Input vcf file (to extract samples list in the same order as in the vcf)', required=False)
parser.add_argument('-s', '--samples', type=str, help='Path to sample file', required=False)

args = parser.parse_args()

if args.samples is None:

    if args.vcf is None:
        raise Exception('[PyPy] ERROR: In the absence of a sample file, a --vcf file must be specified.') 
    
    else:
        cmd = f"bcftools query -l {args.vcf}"

        try :
            call = subprocess.run(cmd, check=True, shell=True, capture_output=True)

        except subprocess.CalledProcessError as cpe:
            sys.exit('[PyPy] ERROR: The following command failed with error code %r:\n[X] => %s\n[X] (STDERR): %r' % (cpe.returncode, cpe.cmd, cpe.stderr))

        samples = np.array(call.stdout.decode('utf8').split())


# write a function to write a samples array from a file

start = time.perf_counter()

def bed_to_array(df):

    inds = []   
    for i, row in df.iterrows():
        called = np.isin(samples, row['ind'].split(','))
        inds.append(called)

    inds = np.array(inds).astype(int)
    start_ends = np.array(df[["start", "end"]], dtype=int)

    np_array = np.concatenate((start_ends, inds), axis = 1)

    return(np_array)


print("Reading in the bed file")

bed_df = pd.read_csv(args.input, sep="\t", usecols=[1, 2, 4], names=['start', 'end', 'ind'], dtype={'start': 'Int64', 'end': 'Int64', 'ind': str})

end = time.perf_counter()

print(f"Finished reading bed file in {end - start} seconds")
print("Starting multiprocessing")

num_processes = args.n_cores
size_bed = len(bed_df)
n_df = ceil(size_bed / num_processes)
list_df = [bed_df[i:i+n_df] for i in range(0,bed_df.shape[0],n_df)]

if __name__=="__main__":

    with mp.Pool(processes=num_processes) as pool:

        results = [pool.apply_async(bed_to_array, (list_df[i],)) for i in range(num_processes)]

        a = results[0].get()

        print(f'Process number 1 obtained')

        for i in range(1, num_processes):
            a = np.concatenate((a, results[i].get()))
            print(f'Process number {i + 1} concatenated')

        end = time.perf_counter()

        print(f"Finished transforming bed to array in {end - start} seconds")
        print("Exporting the numpy array")

        np.savetxt(args.output, a, fmt='%d')

end = time.perf_counter()

print(f'Time elapsed: {end - start} seconds')
