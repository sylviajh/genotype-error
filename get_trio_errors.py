#this has the threshold for variants set to 0.01 (and thus 0.99)
#and it only queries 1 out of every 100 queries
import argparse
import pyigd
import os
import multiprocessing
import time
import pandas as pd
import re
import pickle
#import cProfile

def get_dict(file):
   with open(file,'rb') as f:
      d = pickle.load(f)
   return d


def mendel(child,parent_one,parent_two):
   assert child==1 or child==2 or child ==0
   assert parent_one==1 or parent_one==2 or parent_one==0
   assert parent_one==1 or parent_one==2 or parent_one==0

   if (parent_one==2 and parent_two==2): #both parents have two alternate alleles
      if child==2: #valid assortment is child has two alternate alleles
         return True
      else:
         return False
      
   elif (parent_one==1 and parent_two==2) or (parent_one==2 and parent_two==1): #one parent heterozygous, the other has two alternate alleles
      if child==0: #invalid assortment is child has two reference alleles
         return False
      else:
         return True
      
   elif (parent_one==1 and parent_two==1): #both parents heterzygous
      return True #any assortment is valid
   
   elif (parent_one==1 and parent_two==0) or (parent_one==0 and parent_two==1): #one parent heterozygous, the other has two reference alleles
      if child==2: #invalid assortment is child has two alternate alleles
         return False
      else:
         return True
   
   elif (parent_one==2 and parent_two==0) or (parent_one==0 and parent_two==2): #one parent homozygous alternate, the other homozygous reference
      if child==1: #valid assortment is child heterozygous
         return True
      else:
         return False
      
   elif (parent_one==0 and parent_two==0): #both parents have two reference alleles
      if child==0: #valid assortment is child has two reference alleles
         return True
      else:
         return False
      
   else:
      print(f"parent one: {parent_one}, parent two: {parent_two}, child: {child}",flush=True)
      raise Exception

def get_allele(samples,hap1, hap2):
    assert hap1%2 == 0
    assert hap2%2 == 1
    assert hap1 + 1 == hap2
    if (hap1 in samples) and (hap2 in samples): #individual has two alternate alleles
      return 2
    elif (hap1 in samples) or (hap2 in samples): #individual has one alternate allele
      return 1
    else: #individual has no alternate alleles
      return 0
    
def window(IGD,index,windowsize,pickle_file):
    pos = []
    child_sample = []

    start_time = time.time()
    start_index = index * windowsize
    end_index = start_index + windowsize

    with pyigd.IGDFile(IGD) as igd_file: #open(IGD,"rb") as f:
      #igd_file = pyigd.IGDReader(f)
      total_variants = igd_file.num_variants
      if end_index > total_variants:
         end_index = total_variants

      trio_data = get_dict(pickle_file)

      for variant in range(start_index,end_index):
         samples = igd_file.get_samples(variant)[2]
         
         for child, parents in trio_data.items():
            child_allele = get_allele(samples,child * 2, child * 2 + 1)
            parent_one_allele = get_allele(samples,parents[0] * 2, parents[0] * 2 + 1)
            parent_two_allele = get_allele(samples,parents[1] * 2, parents[1] * 2 + 1)
            is_valid = mendel(child_allele,parent_one_allele,parent_two_allele)
            if not is_valid:
               pos.append(igd_file.get_position_and_flags(variant)[0])
               child_sample.append(child)
    
    errors = pd.DataFrame({'pos': pos, 'child sample': child_sample})
    errors.to_csv(f"errors_{index}.csv")
      
    end_time=time.time()
    print(f"time to parse igd and go through a window: {(end_time)-(start_time)} seconds",flush=True)
    
  
def get_input(IGD, pickle_file):
    num_windows = 72
    with pyigd.IGDFile(IGD) as igd_file: #open(IGD,"rb") as f:
      #igd_file = pyigd.IGDReader(f)
      total_variants = igd_file.num_variants
      window_size = total_variants // (num_windows - 1)
      
    index = [i for i in range(num_windows)]
    IGD_files = [IGD for _ in range(num_windows)]
    windowsizes = [window_size for _ in range(num_windows)]
    pickle_files = [pickle_file for _ in range(num_windows)]

    return zip(IGD_files,index,windowsizes,pickle_files)

def multipool(input):
    number_of_cores = 72
    with multiprocessing.Pool(number_of_cores) as pool:
        # distribute computations and collect results:
        pool.starmap(window, input)

def main():
    parser = argparse.ArgumentParser(description="getting ground truth errors using trio data")
   #  parser.add_argument('-w', '--window_size', type=int, required=True, help="window size to run (helps parallelize)")
    parser.add_argument('-i', '--IGDfile', required=True, help="igd file of data")
    parser.add_argument('-p', '--pickle_file', required =True, help = "pickle file containing dictionary where keys are children, values are parents")
    args = parser.parse_args()
    
    multipool(get_input(args.IGDfile, args.pickle_file))

if __name__ == "__main__":
    main()