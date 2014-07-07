

Consider the following genotype field with 43 values:

    2011011000000100210010101101111110111111010

The bitmaps would be:

    0100100111111011001101010010000001000000101
    0011011000000100010010101101111110111111010
    1000000000000000100000000000000000000000000
    0000000000000000000000000000000000000000000

To represent this as an int we would break it up in the 32 bit chunks and pad

    01001001111110110011010100100000 01000000101000000000000000000000
    00110110000001000100101011011111 10111111010000000000000000000000
    10000000000000001000000000000000 00000000000000000000000000000000
    00000000000000000000000000000000 00000000000000000000000000000000
    |-32---------------------------| |-32---------------------------|
                                                |--22 bit int pad---|

Giving:

    1241199904 1084227584
    906250975 3208642560   
    2147516416 0
    0 0

If we were to directly compress these values using WAH we would first need to
break up the 32 bits in the 31-bit chunks, with more padding. Given

    01001001111110110011010100100000010000001010
    |-43---------------------------------------|
    0100100111111011001101010010000001000000101000000000000000000000
    |-32---------------------------||-32---------------------------|
                                               |--22 bit int pad---|
    010010011111101100110101001000000100000010100000000000000000000000000000000000000000000000000
    |-31--------------------------||-31--------------------------||-31--------------------------|
                                               |--22 bit int pad---||-29 bit wah pad------------|

This gives us a total padding of 51, which is larger than the original value
size.  If we instead only consider the non-padded bits we would have the following:

    01001001111110110011010100100000010000001010
    |-43---------------------------------------|
    01001001111110110011010100100000010000001010000000000000000000
    |-31--------------------------||-31--------------------------|
                                                |-18 bit wah pad-|

Which goes from 1.2x overhead down to 0.4x overhead.  We can do this by
carrying the size of the input (not the size of the data structure) through the
process and stoping each step when the last field has been considered.

Starting from:

    2011011000000100210010101101111110111111010

we get the uncompress binary values:

    10000101000101000000000000010000 -> 2232680464
    10010000010001000101000101010101 -> 2420396373
    01000101010101010001000000000000 -> 1163202560 

by:
    
    char *plt = "2 0 1 1 0 1 1 0 0 0 0 0 0 1 0 0 "
                "2 1 0 0 1 0 1 0 1 1 0 1 1 1 1 1 "
                "1 0 1 1 1 1 1 1 0 1 0";

    unsigned int *ints;
    unsigned int int_len = plt_line_to_packed_ints(plt, 43, &ints);

In this case there are 43 fields, so 43 is passed.

Get the bitmaps:

    0100100111111011001101010010000001000000101
    0011011000000100010010101101111110111111010
    1000000000000000100000000000000000000000000
    0000000000000000000000000000000000000000000

by:

    unsigned int *int_bm;
    unsigned int int_bm_len = ubin_to_bitmap(ints, int_len, 86, &int_bm);

In this case there are 43 fields, but this opperation is concerned with the
number of used bits.  Uncompressed binary used 2 bits per field, for a total of
86 bits.  We still need two integers to represent each bitmap, which will leave
22 bits of padding.

    0100100111111011001101010010000001000000101000000000000000000000
    0011011000000100010010101101111110111111010000000000000000000000
    1000000000000000100000000000000000000000000000000000000000000000
    0000000000000000000000000000000000000000000000000000000000000000
    |-32---------------------------||-32---------------------------|
                                               |--22 bit int pad---|

The first step in converting these values to WAH encoded bitmaps is to break
each value up into a group of 31-bit values.  By passing the number of used bits
to this step, we can prevent the padding from growing.  

Without setting used bits (padding shown as o):

    0100100111111011001101010010000001000000101ooooooooooooooooooooo
    |-32---------------------------||-32---------------------------|
                                               |--22 bit int pad---|

    0100100111111011001101010010000
    001000000101ooooooooooooooooooo
    oo00000000000000000000000000000
    |-31--------------------------|
      |-29 bit wah pad------------|

With setting used bits (padding shown as o):

    0100100111111011001101010010000
    001000000101ooooooooooooooooooo
    |-31--------------------------|
 
Continue this for all four bitmaps we wet:

    31-bit groups:
    1:
    0100100111111011001101010010000
    0010000001010000000000000000000
    2:
    0011011000000100010010101101111
    1101111110100000000000000000000
    3:
    1000000000000000100000000000000
    0000000000000000000000000000000
    4:
    0000000000000000000000000000000
    0000000000000000000000000000000
    |-31--------------------------|

by:

    unsigned int *O;
    unsigned int O_len = map_from_32_bits_to_31_bits(&(ints[0]), 2, 43, &O);

Here, we opperate on one bitmap at a time and since there are 43 values, and
each value is represented by a single bit, we pass 43.

The compressed values are:

    1:
    00100100111111011001101010010000
    00010000001010000000000000000000
    2:
    00011011000000100010010101101111
    01101111110100000000000000000000
    3:
    01000000000000000100000000000000
    00000000000000000000000000000000
    4:
    10000000000000000000000000000010
