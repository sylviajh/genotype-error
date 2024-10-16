import pyigd

with open("ukb20279_c1_b0_v1.v3.igd") as f:
  igd_file = pyigd.IGDReader(f)
  print(igd_file.get_individual_ids())