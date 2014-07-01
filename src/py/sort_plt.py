import sys

file_name = sys.argv[1]

f = open(file_name, 'r')

num_fields = -1
num_records = -1
col_sums = []
for l in f:
    if num_fields == -1:
        num_fields = int(l)
        col_sums = [0] * num_fields
    elif num_records == -1:
        num_records = int(l)
    else:
        i = 0
        for x in l.split():
            col_sums[i] += int(x)
            i+=1

