# KANBAN

## GOALS

### Urgent

### High Value

* Get a 1KG data set to validate results
    * Download 1KG VCF
    * Convert to:
        * sorted ind.txt
        * ind.txt
        * ind.txt.bgz/tabix
        * ind.wah
        * ped
        * bed
        * vcf
        * bcf

* Query time tests
    * WAH
    * RLE
    * plain text
    * grabix
    * uncompressed binary
    * vcf
    * bcf

* Funcitonal Tests
    * WAH
        * Convert to WAH
    * RLE
    * plain text
    * grabix
    * uncompressed binary
    * vcf
    * bcf

### Improvements

* generalized to 32 or 64 bit 
    * WAH
    * uncompressed binary

### Brain Storm

* Create an N-bit WAH compression scheme
* Create a compression scheme that enables efficient (parallel?) sweeps
    * This may be a very good GPU target

## QUEUE
* grabix

## IN PROGRESS

### Ryan
#### Develop
#### Document
#### Verify

## DONE
* file conversion
    * plt to ubin
    * ubin to wah
    * ubin to wahbm
* file print
    * plt
    * ubin packed int
    * ubin plt
    * wahbm plt
    * wah plt

