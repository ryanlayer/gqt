# Uncompressed plain text

Rows:
# of fields (F)
# of records (R)
record_1/field_1 record_1/field_2 ... record_1/field_F
...
record_R/field_1 record_R/field_2 ... record_R/field_F

# WAH-encode bitmaps

Everything is stored as a 32-bit int.  The bitmap index gives the offsets for
each WAH-endocde bitmap.  The values stored are index AFTER the current bitmap
ends.  So if the value is 5, then the last word in the bitmap is at position 4.

# VCF

To convert from by-variant plain text to VCF:

    ~/src/genotq/src/py/var_plt_to_vcf.py \
        -v data/10.1e4.var.txt  \
        > data/10.1e4.var.vcf

Then to convert to BCF

    bcftools view -O b data/10.1e4.var.vcf > data/10.1e4.var.bcf
