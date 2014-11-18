Messing around with software transactional memory as implemented in GCC[https://gcc.gnu.org/wiki/TransactionalMemory]

```
make
```
Then to test it out, using 50 threads,
```
$ time ./8mer_count -t50 /tmp/ben_fasta/62_palsa20140929.3kb.fa
Using 50 threads
Reading from /tmp/ben_fasta/62_palsa20140929.3kb.fa
Found AAAA AAAA 44507 times
Found AAAA AAAC 33356 times
Found AAAA AAAG 52578 times

real    27m25.245s
user    609m4.460s
sys    8m51.496s
```
Looking at `top`, it seemed to be using between 5 and 20 threads at any one time.

Using 1 thread
```
$ time ./8mer_count -t1 /tmp/ben_fasta/62_palsa20140929.3kb.fa
Using 1 threads
Reading from /tmp/ben_fasta/62_palsa20140929.3kb.fa
Found AAAA AAAA 44507 times
Found AAAA AAAC 33356 times
Found AAAA AAAG 52578 times

real    9m9.286s
user    8m2.372s
sys    1m4.976s
```

So it seems that parallelisation with STM is actually slowing things down, not speeding them up, at least in this (quite bad) implementation of kmer counting.

Two threads is even worse
```
$ time ./8mer_count -t2 /tmp/ben_fasta/62_palsa20140929.3kb.fa
Using 2 threads
Reading from /tmp/ben_fasta/62_palsa20140929.3kb.fa
Found AAAA AAAA 44507 times
Found AAAA AAAC 33356 times
Found AAAA AAAG 52578 times

real	46m15.390s
user	53m31.304s
sys	14m30.592s
```
