[![Build Status](https://travis-ci.org/torognes/swarm.svg)](https://travis-ci.org/torognes/swarm) [![codecov](https://codecov.io/gh/torognes/swarm/branch/master/graph/badge.svg)](https://codecov.io/gh/torognes/swarm)

# swarm

A robust and fast clustering method for amplicon-based studies.

The purpose of **swarm** is to provide a novel clustering algorithm
that handles massive sets of amplicons. Results of traditional
clustering algorithms are strongly input-order dependent, and rely on
an arbitrary **global** clustering threshold. **swarm** results are
resilient to input-order changes and rely on a small **local** linking
threshold *d*, representing the maximum number of differences between
two amplicons. **swarm** forms stable, high-resolution clusters, with
a high yield of biological information.

To help users, we describe
[a complete pipeline](https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline)
starting from raw fastq files, clustering with **swarm** and producing
a filtered occurrence table.

swarm 3.0 introduces:
* a much faster default algorithm,
* a reduced memory footprint,
* binaries for Windows x86-64, macOS ARM64, GNU/Linux ARM64, and
  GNU/Linux POWER8,
* an updated, hardened, and thoroughly tested code (804 carefully
  crafted black-box tests),
* strict dereplication of input sequences is now **mandatory**,
* `--seeds` option (`-w`) now outputs results sorted by decreasing
  abundance, and then by alphabetical order of sequence labels.

swarm 2.0 introduced several novelties and improvements over swarm
1.0:
* built-in breaking phase now performed automatically,
* possibility to output cluster representatives in fasta format (option
  `-w`),
* fast algorithm now used by default for *d* = 1 (linear time
  complexity),
* a new option called *fastidious* that refines *d* = 1 results and
  reduces the number of small clusters.

## Common misconceptions

**swarm** is a single-linkage clustering method, with some superficial
  similarities with other clustering methods (e.g., [Huse et al,
  2010](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2909393/)). **swarm**'s
  novelty is its iterative growth process and the use of sequence
  abundance values to delineate clusters. **swarm** properly delineates
  large clusters (high recall), and can distinguish clusters with as little as
  two differences between their centers (high precision).

**swarm** uses a local clustering threshold (*d*), not a global
  clustering threshold like other algorithms do. Users may be tempted
  to convert a 97%-global similarity threshold into a number of
  differences, and to use large *d* values. This is not a correct use
  of swarm. Clusters produced by swarm are naturally larger than *d*, and
  tests have shown that using the default *d* value (*d* = 1) gives
  good results on most datasets. Using the new fastidious option
  further improves the quality of results. For long amplicons or
  shallow sequencing, higher *d* values can be used (*d* = 2 or *d* =
  3, very rarely more).

**swarm** produces high-resolution results, especially when using *d*
  = 1. Under certain rare conditions though, a given marker may not
  evolve fast enough to distinguish molecular taxa. If it concerns
  abundant sequences, swarm may form a cluster with a large radius,
  whereas classic clustering methods will cut through randomly,
  forcing delineation where the 97%-threshold falls. So, keep in mind
  that molecular markers have limitations too.


## Quick start

**swarm** most simple usage is:

```sh
./swarm amplicons.fasta
```

That command will apply default parameters (`-d 1`) to the fasta file
`amplicons.fasta`. The fasta file must be formatted as follows:

```
>seqID1_1000
acgtacgtacgtacgt
>seqID2_25
cgtcgtcgtcgtcgt
```

where sequence identifiers are unique and end with a value indicating
the number of occurrences of the sequence (e.g., `_1000`). Alternative
format is possible with the option `-z`, please see the [user
manual](https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf). Swarm
**requires** each fasta entry to present a number of occurrences to
work properly. That crucial information can be produced during the
[dereplication](#dereplication-mandatory) step.

Use `swarm -h` to get a short help, or see the
  [user manual](https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf)
  for a complete description of input/output formats and command line
  options.

The memory footprint of **swarm** is roughly 0.6 times the size of the
input fasta file. When using the fastidious option, memory footprint
can increase significantly. See options `-c` and `-y` to control and
cap swarm's memory consumption.


## Install ##

Get the latest binaries for GNU/Linux, macOS or Windows from [the
release page](https://github.com/torognes/swarm/releases "swarm tagged
releases"). Get the source code from
[GitHub](https://github.com/torognes/swarm "swarm public repository")
using the [ZIP
button](https://github.com/torognes/swarm/archive/master.zip "swarm
zipped folder") or `git clone`, and compile `swarm` with GCC (version
4.8.5 or more recent) or with clang (version 9 or more recent):

```sh
git clone https://github.com/torognes/swarm.git
cd swarm/
make

# or, with clang
make CC="clang-9" CXX="clang++-9"
```

If you have administrator privileges, you can make **swarm**
accessible for all users. Simply copy the binary `./bin/swarm` to
`/usr/local/bin/` or to `/usr/bin/`. The man page can be installed
this way:

```sh
cd ./man/
gzip -c swarm.1 > swarm.1.gz
mv swarm.1.gz /usr/local/share/man/man1/
# or
mv swarm.1.gz /usr/share/man/man1/
```

Once installed, the man page can be accessed with the command `man
swarm`.


## Install with conda ##

(thanks to GitHub user [Gian77](https://github.com/Gian77) for
reporting this procedure)

Assuming you already have a conda set-up ([anaconda or
miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation)),
start by activate an environment with python 3:

```sh
conda activate py3
```

Make sure you have all the necessary channels for the bioconda
packages:

```sh
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

List the different versions of swarm available and install one:

```sh
conda search -c bioconda swarm
conda install -c bioconda swarm=3.0.0=hc9558a2_0
swarm --version  # check
```


## Prepare amplicon fasta files ##

To facilitate the use of **swarm**, we provide examples of shell
commands that can be use to format and check the input fasta
file. Warning, these examples may not be suitable for very large
files.

We assume your SFF or FASTQ files have been properly pair-assembled
(with [vsearch](https://github.com/torognes/vsearch) for example),
trimmed from adaptors and primers (with
[cutadapt](https://code.google.com/p/cutadapt/) for example), and
converted to fasta.


### Dereplication (mandatory) ###

In a sample, or collection of sample, a given sequence may appear
several times. That number of strictly identical occurrences
represents the *abundance* value of the sequence. Swarm requires all
fasta entries to present abundance values to be able to produce
high-resolution clusters, like this:

```
>seqID1_1000
acgtacgtacgtacgt
>seqID2_25
cgtcgtcgtcgtcgt
```

were `seqID1` has an abundance of 1,000 and `seqID2` has an abundance
of 25 (alternative formats are possible, please see the
[user manual](https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf)).

The role of the dereplication step is to identify, merge and sort
identical sequences by decreasing abundance. Here is a command using
[vsearch](https://github.com/torognes/vsearch) v1.3.3 or superior:

```sh
vsearch \
    --derep_fulllength amplicons.fasta \
    --sizeout \
    --relabel_sha1 \
    --fasta_width 0 \
    --output amplicons_linearized_dereplicated.fasta
```

The command performs the dereplication, the linearization
(`--fasta_width 0`) and the renaming with hashing values
(`--relabel_sha1`). If you can't or don't want to use vsearch, here is
an example using standard command line tools:

```sh
grep -v "^>" amplicons_linearized.fasta | \
grep -v [^ACGTacgt] | sort -d | uniq -c | \
while read abundance sequence ; do
    hash=$(printf "${sequence}" | sha1sum)
    hash=${hash:0:40}
    printf ">%s_%d_%s\n" "${hash}" "${abundance}" "${sequence}"
done | sort -t "_" -k2,2nr -k1.2,1d | \
sed -e 's/\_/\n/2' > amplicons_linearized_dereplicated.fasta
```

Amplicons containing characters other than "ACGT" are discarded. The
dereplicated amplicons receive a meaningful unique name (hash values),
and are sorted by decreasing number of occurrences and by hash values
(to guarantee a stable sorting). The use of a hashing function also
provides an easy way to compare sets of amplicons. If two amplicons
from two different sets have the same hash code, it means that the
sequences they represent are identical.

If for some reason your fasta entries don't have abundance values, and
you still want to run swarm (not recommended), you can specify a
default abundance value with **swarm**'s `--append-abundance` (`-a`)
option to be used when abundance information is missing from a
sequence.


### Launch swarm ###

Here is a typical way to use **swarm**:

```sh
./swarm -f -t 4 -w cluster_representatives.fasta amplicons.fasta > /dev/null
```

**swarm** will partition your dataset with the finest resolution
(local number of differences *d* = 1 by default, built-in elimination
of potential chained clusters, fastidious processing) using 4
CPU-cores. cluster representatives will be written to a new fasta file,
other results will be discarded (`/dev/null`).

See the
[user manual](https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf)
for details on swarm's options and parameters.


## Frequently asked questions ##

To facilitate the use of **swarm**, we provide examples of options or
shell commands that can be use to parse **swarm**'s output. We assume
that the amplicon fasta file was prepared as describe above
(linearization and dereplication).


### Refine swarm clusters ###

The chain-breaking, which used to be performed in a second step in
**swarm** 1.0, is now built-in and performed by default. It is
possible to deactivate it with the `--no-otu-breaking` option, but it
is not recommended. The fastidious option is recommended when using
*d* = 1, as it will reduce the number of small clusters while maintaining
a high clustering resolution. The principle of the fastidious option
is described in the figure below:


![](https://github.com/frederic-mahe/swarm/blob/master/figures/swarm_2.0_fastidious_reduced.png)


### Count the number of amplicons per cluster ###

You might want to check the size distribution of clusters (number of
amplicons in each cluster), and count the number of singletons (clusters
containing only one amplicon). It can be easily done with the
`--statistics-file filename` option. Each line in the output file
represents a cluster and provides different metrics. See the manual for a
complete description.


### Get the seed sequence for each cluster ###

It is frequent for subsequent analyses to keep only one representative
amplicon per cluster (usually the seed) to reduce the computational
burden. That operation is easily done with **swarm** by using the `-w
filename` option.


### Get fasta sequences for all amplicons in a cluster ###

For each cluster, get the fasta sequences for all amplicons. Warning, this
loop can generate a very large number of files. To limit the number of
files, a test can be added to exclude swarms with less than *n*
elements. See
[this wiki page](https://github.com/torognes/swarm/wiki/Get-fasta-sequences-for-all-amplicons-in-a-cluster)
for more examples.

```sh
INPUT_SWARM="amplicons.swarms"
INPUT_FASTA="amplicons.fasta"
OUTPUT_FOLDER="swarms_fasta"
AMPLICONS=$(mktemp)
mkdir "${OUTPUT_FOLDER}"
while read CLUSTER ; do
    tr " " "\n" <<< "${CLUSTER}" | sed -e 's/^/>/' > "${AMPLICONS}"
    seed=$(head -n 1 "${AMPLICONS}")
    grep -A 1 -F -f "${AMPLICONS}" "${INPUT_FASTA}" | \
        sed -e '/^--$/d' > "./${OUTPUT_FOLDER}/${seed/>/}.fasta"
done < "${INPUT_SWARM}"
rm "${AMPLICONS}"
```


## Citation ##

To cite **swarm**, please refer to:

**Swarm: robust and fast clustering method for amplicon-based studies.**<br />
Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014)<br />
PeerJ 2:e593 doi: [10.7717/peerj.593](http://dx.doi.org/10.7717/peerj.593)

**Swarm v2: highly-scalable and high-resolution amplicon clustering.**<br />
Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2015)<br />
PeerJ 3:e1420 doi: [10.7717/peerj.1420](http://dx.doi.org/10.7717/peerj.1420)

**Swarm v3: towards tera-scale amplicon clustering.**<br />
Mahé F, Czech L, Stamatakis A, Quince C, de Vargas C, Dunthorn M, Rognes T. (2021)<br />
Bioinformatics doi: [10.1093/bioinformatics/btab493](https://doi.org/10.1093/bioinformatics/btab493)


## Acknowledgments ##

Many thanks to the following people for their valuable contributions:

* Lucas Czech
* Etienne Platini


## Contact ##

You are welcome to:

* submit suggestions and bug-reports at: https://github.com/torognes/swarm/issues
* send a pull request on: https://github.com/torognes/swarm/
* compose a friendly e-mail to: Frédéric Mahé <frederic.mahe@cirad.fr> and Torbjørn Rognes <torognes@ifi.uio.no>


## Third-party pipelines ##

**swarm** is available in third-party pipelines:

* [FROGS](https://github.com/geraldinepascal/FROGS): a
  [Galaxy](https://galaxyproject.org/)/CLI workflow designed to
  produce a cluster count matrix from high depth sequencing amplicon
  data.
* [LotuS (v1.30)](http://psbweb05.psb.ugent.be/lotus/): extremely fast
  cluster building, annotation, phylogeny and abundance matrix pipeline,
  based on raw sequencer output.
* [QIIME (v1.9)](http://qiime.org/): a multi-purpose pipeline for
  performing microbiome analysis from raw DNA sequencing data.


## Alternatives ##

If you want to try alternative free and open-source clustering
methods, here are some links:

* [vsearch](https://github.com/torognes/vsearch)
* [Oligotyping](http://merenlab.org/projects/oligotyping/)
* [DNAclust](http://dnaclust.sourceforge.net/)
* [Sumaclust](http://metabarcoding.org/sumatra)
* [Crunchclust](https://code.google.com/p/crunchclust/)


## Roadmap ##

swarm adheres to [semantic versioning 2.0.0](https://semver.org/):

> Given a version number MAJOR.MINOR.PATCH, increment the:
>
> MAJOR version when you make incompatible API changes,
> MINOR version when you add functionality in a backwards compatible manner, and
> PATCH version when you make backwards compatible bug fixes.

swarm 3.1.x:
- use more C++11 and STL features,
- eliminate most of clang-tidy's warnings,
- measure the effect of code modernization on run-time performances

swarm 3.1.y:
- refactor to reduce cyclomatic complexity (simpler and shorter functions),
- reduce/eliminate linuxisms to improve portability

swarm 3.2.0:
- swarm can be compiled natively on a BSD or a Windows system

swarm 4.0.0:
- drop compatibility with GCC 4 and GCC 5 (late-2024 GCC 8 will become
  the new *de facto* standard for HPC centers),
- start using C++14 and C++17 features,
- rename option `-n` to `--no-cluster-breaking` (API change)


## Version history ##

### version 3.1.x ###

**swarm** 3.1.x fixes a minor bug, improves code and documentation,
and eliminates compilation warnings and static analysis warnings:
- add: more compilation checks (`shadow`, `useless-cast`,
  `conversion`, `sign-conversion`),
- fix: minor bug introduced in version 3.1.1 (alloc-dealloc-mismatch;
  bug had no impact on clustering results),
- fix: 50 warnings triggered by newly added compilation checks,
- fix: x clang-tidy warnings (from 2,629 warnings, down to x),
- improve: documentation for output option `--network_file`,
- improve: build target platform detection,
- improve: code modernization for long-term maintenance,
- improve: general performance? (case d = 0)

### version 3.1.4 ###

**swarm** 3.1.4 fixes a minor bug, eliminates compilation warnings and
static analysis warnings, and improves code:
- fix: add checks to prevent silent overflow of short unsigned integers,
- fix: compilation warnings with GCC 13 and clang 18,
- fix: 1,040 clang-tidy warnings (from 3,669 warnings, down to 2,629),
- improve: code modernization for long-term maintenance,
- improve: double the maximal number of threads (from 256 threads to 512),
- improve: make `-DNDEBUG` the default compilation behavior,
- performance: stable for all modes, except a 6 to 10% increase in
  memory footprint when d > 2

### version 3.1.3 ###

**swarm** 3.1.3 fixes a few minor bugs, removes warnings, and improves code
and documentation:
- fix: bug introduced in version 3.1.1, that caused swarm to allocate way too
much memory when d > 1 (bug had no impact on clustering results),
- fix: off-by-one error when allocating memory for a Bloom filter (bug had no
impact on clustering results),
- fix: compilation warning with GCC 12 (and more recent) when using link-time
optimization,
- fix: compilation warning with clang 13 (and more recent): unused set
variable,
- fix: five clang-tidy warnings (readability-braces-around-statements),
- fix: minor code refactoring,
- improve: more uniform vocabulary throughout swarm's documentation (code,
help, manpage, README, companion scripts and wiki),
- improve: code coverage of our test suite (swarm-tests).

### version 3.1.2 ###

**swarm** 3.1.2 fixes a bug with fastidious mode introduced in version
3.1.1, that could cause Swarm to crash. Probably due to allocating too
much memory.

### version 3.1.1 ###

**swarm** 3.1.1 eliminates a risk of segmentation fault with extremely
long sequence headers. Documentation and error messages have been
improved, and code cleaning continued.

### version 3.1 ###

**swarm** 3.1 includes a fix for a bug in the 16-bit SIMD alignment
code that was exposed with a combination of d>1, long sequences, and
very high gap penalties. The code has also been been cleaned up,
tested and improved substantially, and it is now fully C++11
compliant. Support for macOS on Apple Silicon (ARM64) has been added.

### version 3.0 ###

**swarm** 3.0 is much faster when _d_ = 1, and consumes less memory.
Strict dereplication is now mandatory.

### version 2.2.2 ###

**swarm** 2.2.2 fixes a bug causing Swarm to wait forever in very rare
cases when multiple threads were used.

### version 2.2.1 ###

**swarm** 2.2.1 fixes a memory allocation bug for d=1.

### version 2.2.0 ###

**swarm** 2.2.0 fixes several problems and improves usability.
Corrected output to structure and uclust files when using fastidious
mode. Corrected abundance output in some cases. Added check for
duplicated sequences and fixed check for duplicated sequence
IDs. Checks for empty sequences. Sorts sequences by additional fields
to improve stability. Improves compatibility with compilers and
operating systems.  Outputs sequences in upper case. Allows 64-bit
abundances. Shows message when waiting for input from stdin. Improves
error messages and warnings. Improves checking of command line
options. Fixes remaining errors reported by test suite. Updates
documentation.

### version 2.1.13 ###

**swarm** 2.1.13 removes a bug in the progress bar when writing seeds.

### version 2.1.12 ###

**swarm** 2.1.12 removes a debugging message.

### version 2.1.11 ###

**swarm** 2.1.11 fixes two bugs related to the SIMD implementation
of alignment that might result in incorrect alignments and scores.
The bug only applies when d>1.

### version 2.1.10 ###

**swarm** 2.1.10 fixes two bugs related to gap penalties of
alignments.  The first bug may lead to wrong aligments and similarity
percentages reported in UCLUST (.uc) files. The second bug makes Swarm
use a slightly higher gap extension penalty than specified. The
default gap extension penalty used have actually been 4.5 instead of
4.

### version 2.1.9 ###

**swarm** 2.1.9 fixes a problem when compiling with GCC version 6.

### version 2.1.8 ###

**swarm** 2.1.8 fixes a rare bug triggered when clustering extremely
short undereplicated sequences. Also, alignment parameters are not
shown when d=1.

### version 2.1.7 ###

**swarm** 2.1.7 fixes more problems with seed output. Ignore CR
  characters in FASTA files. Improved help and error messsages.

### version 2.1.6 ###

**swarm** 2.1.6 fixes problems with older compilers that do not have
the x86intrin.h header file. It also fixes a bug in the output of seeds
with the `-w` option when d>1.

### version 2.1.5 ###

**swarm** 2.1.5 fixes minor bugs.

### version 2.1.4 ###

**swarm** 2.1.4 fixes minor bugs in the swarm algorithm used for d=1.

### version 2.1.3 ###

**swarm** 2.1.3 adds checks of numeric option arguments.

### version 2.1.2 ###

**swarm** 2.1.2 adds the -a (--append-abundance) option to set a
default abundance value to be used when abundance information is
missing from the input file. If this option is not specified, missing
abundance information will result in a fatal error. The error message
in that case is improved.

### version 2.1.1 ###

**swarm** 2.1.1 fixes a bug with the fastidious option that caused it
to ignore some connections between heavy and light swarms.

### version 2.1.0 ###

**swarm** 2.1.0 marks the first official release of swarm 2.

### version 2.0.7 ###

**swarm** 2.0.7 writes abundance information in usearch style when using
options `-w` (`--seeds`) in combination with `-z` (`--usearch-abundance`).

### version 2.0.6 ###

**swarm** 2.0.6 fixes a minor bug.

### version 2.0.5 ###

**swarm** 2.0.5 improves the implementation of the fastidious option
and adds options to control memory usage of the Bloom filter (`-y` and
`-c`). In addition, an option (`-w`) allows to output cluster
representatives sequences with updated abundances (sum of all
abundances inside each cluster). This version also enables dereplication
when `d = 0`.

### version 2.0.4 ###

**swarm** 2.0.4 includes a fully parallelized fastidious option.

### version 2.0.3 ###

**swarm** 2.0.3 includes a working fastidious option.

### version 2.0.2 ###

**swarm** 2.0.2 fixes SSSE3 problems.

### version 2.0.1 ###

**swarm** 2.0.1 is a development release that partially implements the
fastidious option.

### version 2.0.0 ###

**swarm** 2.0.0 simplifies the usage of swarm by using the fast
algorithm and the built-in cluster breaking by default. Some options are
changed and some new output options are introduced.

### version 1.2.21 ###

**swarm** 1.2.21 is supposed to fix some problems related to the use of the
SSSE3 CPU instructions which are not always available.

### version 1.2.20 ###

**swarm** 1.2.20 presents a production-ready version of the
alternative algorithm (option `-a`), with optional built-in cluster
breaking (option `-n`). That alternative algorithmic approach (usable
only with *d* = 1) is considerably faster than currently used
clustering algorithms, and can deal with datasets of 100 million
unique amplicons or more in a few hours. Of course, results are
rigourously identical to the results previously produced with
swarm. That release also introduces new options to control swarm
output (options `-i` and `-l`).

### version 1.2.19 ###

**swarm** 1.2.19 fixes a problem related to abundance information when
  the sequence identifier includes multiple underscore characters.

### version 1.2.18 ###

**swarm** 1.2.18 reenables the possibility of reading sequences from
  `stdin` if no file name is specified on the command line. It also
  fixes a bug related to CPU features detection.

### version 1.2.17 ###

**swarm** 1.2.17 fixes a memory allocation bug introduced in
  version 1.2.15.

### version 1.2.16 ###

**swarm** 1.2.16 fixes a bug in the abundance sort introduced in
  version 1.2.15.

### version 1.2.15 ###

**swarm** 1.2.15 sorts the input sequences in order of decreasing
  abundance unless they are detected to be sorted already. When using
  the alternative algorithm for *d* = 1 it also sorts all subseeds in
  order of decreasing abundance.

### version 1.2.14 ###

**swarm** 1.2.14 fixes a bug in the output with the swarm breaker
  option (`-b`) when using the alternative algorithm (`-a`).

### version 1.2.13 ###

**swarm** 1.2.13 updates the citation.

### version 1.2.12 ###

**swarm** 1.2.12 improves speed of new search strategy for *d* = 1.

### version 1.2.11 ###

**swarm** 1.2.11 corrects the number of differences reported in the
  break swarms output.

### version 1.2.10 ###

**swarm** 1.2.10 allows amplicon abundances to be specified using the
  usearch style in the sequence header (e.g. `>id;size=1`) when the
  `-z` option is chosen. Also fixes the bad URL shown in the previous
  version of swarm.

### version 1.2.9 ###

**swarm** 1.2.9 includes a parallelized variant of the new search
  strategy for *d* = 1. It seems to be fairly scalable up to about 16
  threads for longer reads (~400bp), while up to about 8 threads for
  shorter reads (~150bp). Using about 50% more threads than available
  physical cores is recommended. This version also includes the *d*
  parameter in the beginning of the mothur-style output (e.g.,
  `swarm\_1`). Also, in the break_swarms output the real number of
  differences between the seed and the amplicon is indicated in the
  last column.

### version 1.2.8 ###

**swarm** 1.2.8 fixes an error with the gap extension
  penalty. Previous versions effectively used a gap penalty twice as
  large as intended. This version also introduces an experimental new
  search strategy in the case where *d* = 1 that appears to be almost
  linear and faster at least for datasets of about half a million
  sequences or more. The new strategy can be turned on with the `-a`
  option.

### version 1.2.7 ###

**swarm** 1.2.7 incorporates a few small changes and improvements to
  make it ready for integration into QIIME.

### version 1.2.6 ###

**swarm** 1.2.6 add an option (`-r` or `--mothur`) to format the
  output file as a mothur-compatible list file instead of the native
  swarm format.  When **swarm** encounters an illegal character in the
  input sequences it will now report the illegal character and the
  line number.

### version 1.2.5 ###

**swarm** 1.2.5 can be run on CPUs without the POPCNT feature. It
  automatically checks whether the CPU feature is available and uses
  the appropriate code.  The code that avoids POPCNT is just slightly
  slower. Only basic SSE2 is now required.

### version 1.2.4 ###

**swarm** 1.2.4 changes the name of the new option from
  `--break_swarms` to `--break-swarms` for consistency with other
  options, and also adds a companion script `swarm_breaker.py` to
  refine swarm results (`scripts` folder).

### version 1.2.3 ###

**swarm** 1.2.3 adds an option (`-b` or `--break_swarms`) to output
  all pairs of amplicons to `stderr`. The data can be used for
  post-processing of the results to refine the swarms. The syntax of
  the inline assembly code is also changed for compatibility with more
  compilers.

### version 1.2.2 ###

**swarm** 1.2.2 fixes an issue with incorrect values in the statistics
  file (maximum generation and radius of swarms). This version is also
  a bit faster.

### version 1.2.1 ###

**swarm** 1.2.1 removes the need for a SSE4.1 capable CPU and should
  now be able to run on most servers, desktops and laptops.


### version 1.2.0 ###

**swarm** 1.2.0 introduces a pre-filtering of similar amplicons based
  on *k*-mers. This eliminates most of the time-consuming pairwise
  alignments and greatly improves speed. The speedup can be more than
  100-fold compared to previous swarm versions when using a single
  thread with a large set of amplicons. Using multiple threads induces
  a computational overhead, but becomes more and more efficient as the
  size of the amplicon set increases.


### version 1.1.1 ###

**swarm** now works on Apple computers. This version also corrects an
  issue in the pairwise global alignment step that could lead to
  sub-optimal alignments. Slightly different alignments may result
  relative to previous version, giving slightly different swarms.


### version 1.1.0 ###

**swarm** 1.1.0 introduces new optimizations and is 20% faster than
  the previous version on our test dataset. It also introduces two new
  output options: `statistics` and `uclust-like` format.

By specifying the `-s` option to **swarm** it will now output detailed
  statistics about each swarm to a specified file. It will print the
  number of unique amplicons, the number of occurrences, the name of
  the seed and its abundance, the number of singletons (amplicons with
  an abundance of 1), the number of iterations and the maximum radius
  of the swarm (i.e. number of differences between the seed and the
  furthermost amplicon). When using input data sorted by decreasing
  abundance, the seed is the most abundant amplicon in the swarm.

Some pipelines use the
  [uclust output format](http://www.drive5.com/uclust/uclust_userguide_1_1_579.html#_Toc257997686
  "page describing the uclust output format") as input for subsequent
  analyses. **swarm** can now output results in this format to a
  specified file with the `-u` option.
