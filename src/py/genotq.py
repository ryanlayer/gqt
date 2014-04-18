import struct

vars_per_int = 16

# Breaks a 32-bit integer into an array of 16 integers corresponding to the 16
# 2-bit ints encoded in the orignial int
def var_ubin_to_var_int(var_ubin):
    V = []
    for i in range(30,-1,-2):
        v = (var_ubin >> i) & 3
        V.append(v)
    return V

# Takes a file name
# grabs the first two ints 
# returns an array with the first int (validation int), the second int (number
# of individuals encoded in the file, the number of integers required to
# encode that number of individuals, and the open file handle
def var_ubin_validate(ubin_file):
    f = open(ubin_file,'rb')

    validate = -1
    num_inds = -1
    num_ints = -1

    try:
        #validate
        byte = f.read(4)
        curr_int = struct.unpack('I', byte)[0]
        validate = curr_int

        #get size
        byte = f.read(4)
        curr_int = struct.unpack('I', byte)[0]

        num_inds = curr_int
        num_ints = int((num_inds + vars_per_int - 1)/vars_per_int)

    except AnError as e:
        print 'Exception occurred, value:', e.value
        return [0,0,0,None]

    return [validate,num_inds,num_ints,f]

# Takes an open file halde to varub file, the number of integers in each record
# validation, and an array of ints that correspond to record numbers in the
# file.
# Returns and array of 32-bit intergers corresponding to the record numbers
# provided in the function call.
def var_ubin_get_rows(f, num_ints, R):
    V=[]

    try:
        for r in R:
            cols = []

            f.seek(8 + ((r - 1) * num_ints)*4 )

            byte = f.read(4)

            while byte != "":
                curr_int = struct.unpack('I', byte)[0]

                cols.append(curr_int)

                if len(cols) == num_ints:
                    V.append(cols)
                    break

                byte = f.read(4)
    except AnError as e:
        print 'Exception occurred, value:', e.value
        return None
    return V

# Just close the file
def var_ubin_finalize(f):
    f.close()
