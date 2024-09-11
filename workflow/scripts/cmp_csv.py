import sys
import pandas as pd

# check if two csv's are equal.  Deals with sorting and floating point precision
print(sys.argv[1], sys.argv[2])

old = pd.read_csv(sys.argv[1]).sort_values(['fileName', sys.argv[3]]).reset_index(drop=True)
new = pd.read_csv(sys.argv[2]).sort_values(['fileName', sys.argv[3]]).reset_index(drop=True)

try:
    pd.testing.assert_frame_equal(old, new)
except:
    print("failed assertion, starting pdb")
    breakpoint()
