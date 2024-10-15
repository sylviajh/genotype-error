#this has the threshold for variants set to 0.01 (and thus 0.99)
#and it only queries 1 out of every 100 queries
import argparse
import pyigd
import os
import multiprocessing
import time
import pandas as pd
import re
#import cProfile

def get_dict():
  import re

  trio_dict = {}
  line_count = 0

  with open("tmp.finalDoubleParent.uniqueParents.trios") as trios:
    for line in trios.readlines():
      if line_count == 0:
        line_count += 1
        continue

      split = line.split()
      for i in range(2,len(split)):
        split[i] = re.sub(r"\D", "", split[i])
      
      for j in range(2, len(split)):
        parents = [int(split[0]),int(split[1])]
        trio_dict[int(split[j])] = parents
  return trio_dict

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
      
   elif (parent_one==0 and parent_two==0): #both parents have two reference alleles
      if child==0: #valid assortment is child has two reference alleles
         return True
      else:
         return False
      
   else:
      raise Exception

def get_allele(hap_one,hap_two,id):
    if (id in hap_one) and (id in hap_two): #individual has two alternate alleles
      return 2
    elif (id in hap_one) or (id in hap_two): #individual has one alternate allele
      return 1
    else: #individual has no alternate alleles
      return 0
    
def window(IGD,index,windowsize):
    pos = []
    hap = []

    start_time = time.time()
    start_index = index * windowsize
    end_index = start_index + windowsize
    with open(IGD,"rb") as f:
      igd_file = pyigd.IGDReader(f)
      total_variants = igd_file.num_variants
      if end_index > total_variants - windowsize:
         end_index = total_variants

      trio_data = get_dict()
      #### somehow get from individual id to haplotype id ####
      got_diploid = False
      for variant in range(start_index,end_index):
         samples = igd_file.get_samples(variant)[2]
         if not got_diploid:
            hap_one = samples
         else:
            hap_two = samples
            for child, parents in trio_data.items():
               child_allele = get_allele(hap_one,hap_two,child)
               parent_one_allele = get_allele(hap_one,hap_two,parents[0])
               parent_two_allele = get_allele(hap_one,hap_two,parents[1])
               is_valid = mendel(child_allele,parent_one_allele,parent_two_allele)
               if not is_valid:
                  pos.append(igd_file.get_position_and_flags(variant)[0])
                  hap.append(child) #### this is not right since child is individual and not haplotype, must be changed
    
    errors = pd.DataFrame({'pos': pos, 'hap': hap})
    errors.to_csv(f"errors_{index}.csv")
      
    end_time=time.time()
    print(f"time to parse igd and go through a window: {(end_time)-(start_time)} seconds",flush=True)
    
  
def get_input(IGD,window_size):
    with open(IGD,"rb") as f:
      igd_file = pyigd.IGDReader(f)
      total_variants = igd_file.num_variants
      num_windows = total_variants / window_size
      if total_variants % window_size != 0:
         num_windows += 1
      
    index = [i for i in range(num_windows)]
    IGD_files = [IGD for _ in range(num_windows)]
    windowsizes = [window_size for _ in range(num_windows)]

    return zip(IGD_files,index,windowsizes)

def multipool(input):
    number_of_cores = int(os.environ['SLURM_CPUS_PER_TASK'])
    with multiprocessing.Pool(number_of_cores) as pool:
        # distribute computations and collect results:
        pool.starmap(window, input)

def main():
    parser = argparse.ArgumentParser(description="match admixed haplotypes with reference panel in windows using BK Trees")
    parser.add_argument('-w', '--window_size', type=int, required=True, help="window size to run matching with")
    parser.add_argument('-i', '--IGDfile', required=True, help="igd file")
    args = parser.parse_args()
    
    multipool(get_input(args.IGDfile,args.window_size))

if __name__ == "__main__":
    main()
