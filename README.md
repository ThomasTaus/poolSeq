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
`install.packages("/Path/To/poolSeq_0.3.0.tar.gz", repos=NULL, type="source")`

## Usage

[PoPoolation2]: https://sourceforge.net/projects/popoolation2/
[latest release]: https://github.com/ThomasTaus/poolSeq/releases