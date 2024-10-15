# make dictionary where each key is a child in the trio data and values are the parents

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

