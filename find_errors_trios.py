#this has the threshold for variants set to 0.01 (and thus 0.99)
#and it only queries 1 out of every 100 queries
import pandas as pd
import BKTree
import argparse
import pyigd
import os
import multiprocessing
import time
import numpy as np
import collections
#import cProfile

def calculate_metrics(df):
    num_true_errors = 0

    num_predicted_errors1 = 0
    overlap1 = 0

    num_predicted_errors2 = 0
    overlap2 = 0

    num_predicted_errors3 = 0
    overlap3 = 0

    for i,row in df.iterrows():
        true_errors = row.error_positions
        predicted_errors1 = row.diff_markers1
        predicted_errors2 = row.diff_markers2
        predicted_errors3 = row.diff_markers3
        num_true_errors += len(true_errors)
        num_predicted_errors1 += len(predicted_errors1)
        num_predicted_errors2 += len(predicted_errors2)
        num_predicted_errors3 += len(predicted_errors3)
    
        for err in true_errors:
            if err in predicted_errors1:
                overlap1 += 1
            if err in predicted_errors2:
                overlap2 += 1
            if err in predicted_errors3:
                overlap3 += 1
    #print(num_true_errors,num_predicted_errors1,overlap1,num_predicted_errors2,overlap2,num_predicted_errors3,overlap3,flush=True)
    return num_true_errors,num_predicted_errors1,overlap1,num_predicted_errors2,overlap2,num_predicted_errors3,overlap3

def consensus_string(strings):
  consensus1 = ''
  consensus2 = ''
  consensus3 = ''
  for i in range(len(strings[0])):
      column = [s[i] for s in strings]
      counter = collections.Counter(column)
      most_common = counter.most_common(2)
      if len(most_common) == 1:
        consensus1 += most_common[0][0]
        consensus2 += most_common[0][0]
      elif most_common[1][1] == most_common[0][1]:
        consensus1 += "0"
        consensus2 += "1"
      else:
        consensus1 += most_common[0][0]
        consensus2 += most_common[0][0]
      consensus3 += most_common[0][0]
      
  return consensus1,consensus2,consensus3

def find_diffs(row):
    diff_markers1=[row['start_index']+i for i in range(len(row['query'])) if \
                    row['query'][i] != row['consensus_seq1'][i]]
    diff_markers2=[row['start_index']+i for i in range(len(row['query'])) if \
                    row['query'][i] != row['consensus_seq2'][i]]
    diff_markers3=[row['start_index']+i for i in range(len(row['query'])) if \
                    row['query'][i] != row['consensus_seq3'][i]]
   
    return pd.Series([diff_markers1,diff_markers2,diff_markers3])

def get_consensus_seq(matches):
    l_s=[''.join([str(e) for e in x]) for x in matches]
    return consensus_string(l_s)

def hamming_distance(vector1, vector2):
    return sum(el1 != el2 for el1, el2 in zip(vector1, vector2))

def collate(matches,igd_file,csv_file, start_index,end_index):
    results = pd.DataFrame.from_dict(matches)
    results[['consensus_seq1','consensus_seq2','consensus_seq3']] = results['matches'].apply(lambda x: pd.Series(get_consensus_seq(x)))
    results[['diff_markers1','diff_markers2','diff_markers3']] = results.apply(lambda x: find_diffs(x),axis=1)

    edits = pd.read_csv(csv_file)

    with (pyigd.IGDFile(igd_file) as igd):

        start_pos = igd.get_position_and_flags(start_index)[0]
        end_pos = igd.get_position_and_flags(end_index)[0]
        within = edits[(edits.pos>= start_pos) & (edits.pos < end_pos)].copy()
        
        variant_positions = []
        for i in range(start_index,end_index):
            variant_positions.append(igd.get_position_and_flags(i)[0])

        variant_positions = np.array(variant_positions)
        within['error_positions'] = within['pos'].apply(lambda x: np.searchsorted(variant_positions,x) + start_index)
        grouped_hap = within.groupby('hap')['error_positions'].apply(list)

        final_results = pd.merge(results, grouped_hap, on='hap', how='left')
        final_results['error_positions'] = final_results['error_positions'].apply(lambda x: x if isinstance(x, list) else [])

        # Calculate metrics for simulated and rephased dataframes
        return calculate_metrics(final_results)


def match(YRIgenotypes,start_index,end_index,igd_file,csv_file):
    """
    Creates matches along the window, one for each query genotype in the admixed set 
    """
    # Initialize an empty BKTreeNode
    root_node = BKTree.BKTreeNode.make_empty()

    #Insert reference vectors into the BK tree
    for i,genotype in enumerate(YRIgenotypes):
        BKTree.bk_tree_insert(root_node, ['hap{}'.format(i)], [int(x) for x in genotype], hamming_distance)

    #Look up nearest neighbors for each ADMIXgenotype
    matches=[]
    for i,query_vector in enumerate(YRIgenotypes): 
        if np.random.uniform(0, 1) <= 0.3: 
            results, dist_best, exact_matches=BKTree.bk_tree_lookup7(root_node, query_vector, hamming_distance)
            matches.append({\
                    'start_index':start_index,\
                    'end_index':end_index,\
                    'query':''.join([str(x) for x in query_vector]),\
                    'matches':[x.vector for x in results],\
                    'edit_distance':dist_best, \
                    'exact_matches':[x.vector for x in exact_matches],\
                    'hap': (i+1)\
                    })

    return collate(matches,igd_file,csv_file, start_index,end_index)

def window(YRI,start_index,end_index,csv_file):
    start_time = time.time()
    with pyigd.IGDFile(YRI) as Y:
        YRIgenotypes=[[0 for _ in range(end_index-start_index)] for _ in range(Y.num_samples)]

        for i in range(end_index-start_index): #it is important not to read the end marker itself to avoid overlaps
            for variant in Y.get_samples(i+start_index)[2]:
                YRIgenotypes[variant][i] = 1

    #profile = cProfile.Profile()
    #profile.runcall(match,YRIgenotypes,CEUgenotypes,ADMIXgenotypes,start_index,end_index,outprefix)
    #profile.dump_stats(f'profile.stats{index}')
    num_true_errors,num_predicted_errors1,overlap1,num_predicted_errors2,overlap2,num_predicted_errors3,overlap3 = match(YRIgenotypes,start_index,end_index,YRI,csv_file)
    end_time=time.time()
    print(f"time to parse igd and match a window: {(end_time)-(start_time)} seconds",flush=True)
    return num_true_errors,num_predicted_errors1,overlap1,num_predicted_errors2,overlap2,num_predicted_errors3,overlap3

def get_input(YRI,window_size,csv_file):

    with pyigd.IGDFile(YRI) as igd_file:
      df_positions = [0]
      counter = window_size+igd_file.get_position_and_flags(0)[0]
      for i in range(1,igd_file.num_variants-1):
          if (igd_file.get_position_and_flags(i+1)[0] < counter):
              continue
          else:
              df_positions.append(i+1)
              df_positions.append(i+1)
              counter += window_size

    num_windows = len(df_positions)
    positions_iter = iter(df_positions)
    Y_files = [YRI for _ in range(num_windows)]
    csv_files = [csv_file for _ in range(num_windows)]

    return zip(Y_files,positions_iter,positions_iter,csv_files)

def multipool(input):
    number_of_cores = int(os.environ['SLURM_CPUS_PER_TASK'])

    num_true_errors = 0

    num_predicted_errors1 = 0
    overlap1 = 0

    num_predicted_errors2 = 0
    overlap2 = 0

    num_predicted_errors3 = 0
    overlap3 = 0

    with multiprocessing.Pool(number_of_cores) as pool:
        # distribute computations and collect results:

        results = pool.starmap(window, input)

        for result in results:
            num_true_errors += result[0]

            num_predicted_errors1 += result[1]
            overlap1 += result[2]

            num_predicted_errors2 += result[3]
            overlap2 += result[4]

            num_predicted_errors3 += result[5]
            overlap3 += result[6]
    
    print(f"Number of true errors: {num_true_errors}",flush=True)

    print(f"Number of predicted errors with rng = 0: {num_predicted_errors1}",flush=True)
    print(f"Number of overlap with rng = 0: {overlap1}",flush=True)

    print(f"Number of predicted errors with rng = 1: {num_predicted_errors2}",flush=True)
    print(f"Number of overlap with rng = 1: {overlap2}",flush=True)

    print(f"Number of predicted errors with rng random: {num_predicted_errors3}",flush=True)
    print(f"Number of overlap with rng random: {overlap3}",flush=True)


def main():
    parser = argparse.ArgumentParser(description="match admixed haplotypes with reference panel in windows using BK Trees")
    parser.add_argument('-w', '--window_size', type=int, required=True, help="window size to run matching with")
    parser.add_argument('-i', '--IGDfile', required=True, help="igd file")
    parser.add_argument('-c','--csv_file',required=True,help="csv file of errors obtained using trio data")
    args = parser.parse_args()
    
    multipool(get_input(args.IGDfile,args.window_size,args.csv_file))

if __name__ == "__main__":
    main()
