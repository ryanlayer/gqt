# Grabix-indexed bgz files

### Compress the plain text file

    bgzip -c 100.1e8.ind.txt > 100.1e8.ind.txt.bgz

Check the size difference

    ls -l 100.1e8.ind.txt 100.1e8.ind.txt.bgz

    -rw-r--r--+ 1 rl6sf  staff  117766011 Jun  3 13:56 100.1e8.ind.txt
    -rw-r--r--+ 1 rl6sf  staff   12460596 Jun 17 17:32 100.1e8.ind.txt.bgz

### Index with grabix, and use the index to grab arbitrary lines

    grabix index 100.1e8.ind.txt.bgz

    grabix grab 100.1e8.ind.txt.bgz 1

### 

