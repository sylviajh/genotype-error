# make dictionary where each key is a sample in the trio data and keys are the close relatives we 
# associate with it (ie the mother, father, children, and/or siblings)

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
    
    for j in range(len(split)):
      other_members = set(int(split[x]) for x in range(len(split)) if x != j)

      if trio_dict.get(int(split[j])) == None:
        trio_dict[int(split[j])] = other_members
      else:
        trio_dict[int(split[j])] = trio_dict.get(int(split[j])) | (other_members)

    # line_count += 1
    # if line_count == 5:
    #   break

