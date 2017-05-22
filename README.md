# poolSeq
poolSeq is an R-package that allows you to analyze and simulate Pool-Seq time-series data. Its functionality includes estimation of the effective population size and quantification of selective strength, as well as dominance. Besides simulating Pool-Seq data under a specific scenario, you can also load empirical data into R using the common sync-file format (see [PoPoolation2]).

## Installation
Before installing poolSeq you need to make sure that all the __dependencies__ are available:

* R (>= 3.3.1)
* data.table (>= 1.9.4)
* foreach (>= 1.4.2)
* stringi (>= 0.4-1)
* matrixStats (>= 0.14.2)

For now you need to install these manually. Once this is done you can proceed by __downloading__ the [latest release] of poolSeq. After the download you can __install__ poolSeq with the following R command:

```R
install.packages("/Path/To/poolSeq_0.3.0.tar.gz", repos=NULL, type="source")
```

## Usage
### Read sync-files
Synchronized (sync) files contain allele frequencies at specific genomic loci in multiple populations. Suppose you want to load a sync-file containing allele frequency trajectories of 2 populations (F0.R1, F10.R1, F20.R1, F0.R2, F10.R2, F20.R2). The following command allows you to read such file with poolSeq:

```R
mySync <- read.sync(file="/Path/to/data.sync", gen=c(0, 10, 20, 0, 10, 20), repl=c(1, 1, 1, 2, 2, 2), rising = FALSE)
```

Once the data is loaded in R you can easily get allele frequency trajectories for a specific genomic position `pos` on chromosome `chr` in replicate `repl` with `af.traj(mySync, chr, pos, repl)`. Alternatively, allele frequencies and sequence coverage of all genomic loci can be obtained for a given replicate and generation using `af(sync, repl, gen)` or `coverage(sync, repl, gen)`, respectively.

### Simulate Pool-Seq time-series data
poolSeq enables you to simulate allele frequency trajectories of unlinked loci under a Wright-Fisher model. Assuming an effective population size `Ne` of 1,000 diploid individuals, genetic drift at 10,000 unlinked loci with random starting allele frequencies `p0` over 20 generations can be simulated like this:

```R
simTraj <- wf.traj(p0=runif(10000, 0, 1), Ne=1000, t=c(0, 10, 20))
```

Pool-Seq can either be applied using all individuals of a population, or only a subset. Assuming that only 300 out of the 1,000 individuals were selected for Pool-Seq, the resulting sampling noise can be added as follows:

```R
simTraj <- matrix(sample.alleles(simTraj, size=300, mode="individuals", Ncensus=1000), nrow=nrow(simTraj), dimnames=dimnames(simTraj))
```

Then sampling noise at an average sequencing coverage of 60x can be added, drawing coverage depths from a Poisson distribution:

```R
af <- sample.alleles(simTraj, size=60, mode="coverage")
simTraj <- matrix(af$p.smpld, nrow=nrow(simTraj), dimnames=dimnames(simTraj))
simCov <- matrix(af$size, nrow=nrow(simTraj), dimnames=dimnames(simTraj))
```



[PoPoolation2]: https://sourceforge.net/projects/popoolation2/
[latest release]: https://github.com/ThomasTaus/poolSeq/releases