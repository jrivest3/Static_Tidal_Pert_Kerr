#print_h5s.py 
import numpy as np
import pandas as pd
import h5py
import sys
# import os
# from stat import S_IMODE


filepath=sys.argv[1]
print(filepath)
# stat_info = os.stat(filepath)
# permissions = S_IMODE(stat_info.st_mode)
# owner_read = bool(permissions & 0o400)  # Owner read permission
# owner_write = bool(permissions & 0o200)  # Owner write permission
# owner_execute = bool(permissions & 0o100)  # Owner execute permission

# group_read = bool(permissions & 0o040)  # Group read permission
# group_write = bool(permissions & 0o020)  # Group write permission
# group_execute = bool(permissions & 0o010)  # Group execute permission

# other_read = bool(permissions & 0o004)  # Others read permission
# other_write = bool(permissions & 0o002)  # Others write permission
# other_execute = bool(permissions & 0o001)  # Others execute permission

# print(f'{owner_read=},{owner_write=},{owner_execute=},{group_read=},{group_write=},{group_execute=},{other_read=},{other_write=},{other_execute=}')


# df=pd.DataFrame(np.array(h5py.File(filepath)))

# with h5py.File(filepath, 'r') as h5_file:
#     # Assuming 'variable_1' is the dataset you want to read
#     data_array = h5_file['variable_1'][:]
#     df = pd.DataFrame(data_array)

# print(df)
Gset=h5py.File(filepath,'r')

# Gset.visit(print)
dfs=[]
def func(name):
    dfs.append(np.array(name))
# Gset.visit(func)
tau_array=Gset['tau_array']
for name in Gset.values():
    for value in name.values():
        func(value)
print(dfs[0])
# df.to_csv()
