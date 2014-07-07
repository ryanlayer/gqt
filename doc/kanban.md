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
    * Convert to WAH
    * RLE
    * plain text
    * grabix
    * uncompressed binary
    * vcf
    * bcf

* invert major from (records -> fields/ fiels -> records)

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
* greater than query
#### Document
#### Verify

## DONE
* file conversion
    * plt to ubin
    * ubin to wah
    * ubin to wahbm
    * vcf to plt
* file print
    * plt
    * ubin packed int
    * ubin plt
    * wahbm plt
    * wah plt
