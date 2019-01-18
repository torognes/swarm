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
a filtered OTU table.

swarm 2.0 introduces several novelties and improvements over swarm
1.0:
* built-in breaking phase now performed automatically,
* possibility to output OTU representatives in fasta format (option
  `-w`),
* fast algorithm now used by default for *d* = 1 (linear time
  complexity),
* a new option called *fastidious* that refines *d* = 1 results and
  reduces the number of small OTUs,

## Common misconceptions

**swarm** is a single-linkage clustering method, with some superficial
  similarities with other clustering methods (e.g.,
  [Huse et al, 2010](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2909393/)). **swarm**'s
  novelty is its iterative growth process and the use of sequence
  abundance values to delineate OTUs. **swarm** properly delineates
  large OTUs (high recall), and can distinguish OTUs with as little as
  two differences between their centers (high precision).

**swarm** uses a local clustering threshold (*d*), not a global
  clustering threshold like other algorithms do. Users may be tempted
  to convert a 97%-global similarity threshold into a number of
  differences, and to use large *d* values. This is not a correct use
  of swarm. OTUs produced by swarm are naturally larger than *d*, and
  tests have shown that using the default *d* value (*d* = 1) gives
  good results on most datasets. Using the new fastidious option
  further improves the quality of results. For long amplicons or
  shallow sequencing, higher *d* values can be used (*d* = 2 or *d* =
  3, very rarely more).

**swarm** produces high-resolution results, especially when using *d*
  = 1. Under certain rare conditions though, a given marker may not
  evolve fast enough to distinguish molecular taxa. If it concerns
  abundant sequences, swarm may form an OTU with a large radius,
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
format is possible with the option `-z`, please see the
[user manual](https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf). Swarm
**requires** each fasta entry to present a number of occurrences to
work properly. That crucial information can be produced during the
[dereplication](#dereplication-mandatory) step.

Use `swarm -h` to get a short help, or see the
  [user manual](https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf)
  for a complete description of input/output formats and command line
  options.

The memory footprint of **swarm** is roughly 1.6 times the size of the
input fasta file. When using the fastidious option, memory footprint
can increase significantly. See options `-c` and `-y` to control and
cap swarm's memory consumption.


## Install ##

Get the latest binaries for GNU/Linux or macOS from
[the release page](https://github.com/torognes/swarm/releases "swarm
tagged releases"). Get the source code from
[GitHub](https://github.com/torognes/swarm "swarm public repository")
using the
[ZIP button](https://github.com/torognes/swarm/archive/master.zip
"swarm zipped folder") or git, and compile swarm:

```sh
git clone https://github.com/torognes/swarm.git
cd swarm/src/
make
cd ../bin/
```

If you have administrator privileges, you can make **swarm**
accessible for all users. Simply copy the binary to `/usr/bin/`. The
man page can be installed this way:

```sh
cd ./man/
gzip -c swarm.1 > swarm.1.gz
mv swarm.1.gz /usr/share/man/man1/
```

Once installed, the man page can be accessed with the command `man
swarm`.


## Prepare amplicon fasta files ##

To facilitate the use of **swarm**, we provide examples of shell
commands that can be use to format and check the input fasta
file. Warning, these examples may not be suitable for very large
files.

We assume your SFF or FASTQ files have been properly pair-assembled
(with [pear](https://github.com/xflouris/PEAR) for example), trimmed
from adaptors and primers (with
[cutadapt](https://code.google.com/p/cutadapt/) for example), and
converted to fasta.


### Linearization ###

**swarm** accepts wrapped fasta files as well as linear fasta
files. However, linear fasta files where amplicons are written on two
lines (one line for the fasta header, one line for the sequence) are
easier to manipulate. For instance, many post-clustering queries can
be easily done with `grep` when fasta files are linear. You can use
the following code to linearize your fasta files. Code tested with GNU
Awk 4.0.1.

```sh
awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' amplicons.fasta > amplicons_linearized.fasta
```


### Dereplication (mandatory) ###

In a sample, or collection of sample, a given sequence is likely to
appear several times. That number of strictly identical occurrences
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
you still want to run swarm, you can easily add fake abundance values:

```sh
sed '/^>/ s/$/_1/' amplicons.fasta > amplicons_with_abundances.fasta
```

Alternatively, you may specify a default abundance value with
**swarm**'s `--append-abundance` (`-a`) option to be used when
abundance information is missing from a sequence.


### Launch swarm ###

Here is a typical way to use **swarm**:

```sh
./swarm -f -t 4 -w OTU_representatives.fasta amplicons.fasta > /dev/null
```

**swarm** will partition your dataset with the finest resolution
(local number of differences *d* = 1 by default, built-in elimination
of potential chained OTUs, fastidious processing) using 4
CPU-cores. OTU representatives will be written to a new fasta file,
other results will be discarded (`/dev/null`).

See the
[user manual](https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf)
for details on swarm's options and parameters.


## Frequently asked questions ##

To facilitate the use of **swarm**, we provide examples of options or
shell commands that can be use to parse **swarm**'s output. We assume
that the amplicon fasta file was prepared as describe above
(linearization and dereplication).


### Refine swarm OTUs ###

The chain-breaking, which used to be performed in a second step in
**swarm** 1.0, is now built-in and performed by default. It is
possible to deactivate it with the `--no-otu-breaking` option, but it
is not recommended. The fastidious option is recommended when using
*d* = 1, as it will reduce the number of small OTUs while maintaining
a high clustering resolution. The principle of the fastidious option
is described in the figure below:


![](https://github.com/frederic-mahe/swarm/blob/master/figures/swarm_2.0_fastidious_reduced.png)


### Count the number of amplicons per OTU ###

You might want to check the size distribution of OTUs (number of
amplicons in each OTU), and count the number of singletons (OTUs
containing only one amplicon). It can be easily done with the
`--statistics-file filename` option. Each line in the output file
represents an OTU and provides different metrics. See the manual for a
complete description.


### Get the seed sequence for each OTU ###

It is frequent for subsequent analyses to keep only one representative
amplicon per OTU (usually the seed) to reduce the computational
burden. That operation is easily done with **swarm** by using the `-w
filename` option.


### Get fasta sequences for all amplicons in a OTU ###

For each OTU, get the fasta sequences for all amplicons. Warning, this
loop can generate a very large number of files. To limit the number of
files, a test can be added to exclude swarms with less than *n*
elements. See
[this wiki page](https://github.com/torognes/swarm/wiki/Get-fasta-sequences-for-all-amplicons-in-a-OTU)
for more examples.

```sh
INPUT_SWARM="amplicons.swarms"
INPUT_FASTA="amplicons.fasta"
OUTPUT_FOLDER="swarms_fasta"
AMPLICONS=$(mktemp)
mkdir "${OUTPUT_FOLDER}"
while read OTU ; do
    tr " " "\n" <<< "${OTU}" | sed -e 's/^/>/' > "${AMPLICONS}"
    seed=$(head -n 1 "${AMPLICONS}")
    grep -A 1 -F -f "${AMPLICONS}" "${INPUT_FASTA}" | \
        sed -e '/^--$/d' > "./${OUTPUT_FOLDER}/${seed/>/}.fasta"
done < "${INPUT_SWARM}"
rm "${AMPLICONS}"
```


## Troubleshooting ##

If **swarm** exits with an error message saying `This program
requires a processor with SSE2`, your computer is too old to run
**swarm** (or based on a non x86-64 architecture). **swarm** only runs
on CPUs with the SSE2 instructions, i.e. most Intel and AMD CPUs
released since 2004.


## Citation ##

To cite **swarm**, please refer to:

Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014)
Swarm: robust and fast clustering method for amplicon-based studies.
PeerJ 2:e593 doi: [10.7717/peerj.593](http://dx.doi.org/10.7717/peerj.593)

Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2015)
Swarm v2: highly-scalable and high-resolution amplicon clustering.
PeerJ 3:e1420 doi: [10.7717/peerj.1420](http://dx.doi.org/10.7717/peerj.1420)


## Contact ##

You are welcome to:

* submit suggestions and bug-reports at: https://github.com/torognes/swarm/issues
* send a pull request on: https://github.com/torognes/swarm/
* compose a friendly e-mail to: Frédéric Mahé <mahe@rhrk.uni-kl.de> and Torbjørn Rognes <torognes@ifi.uio.no>


## Third-party pipelines ##

**swarm** is available in third-party pipelines:

* [FROGS](https://github.com/geraldinepascal/FROGS): a
  [Galaxy](https://galaxyproject.org/)/CLI workflow designed to
  produce an OTU count matrix from high depth sequencing amplicon
  data.
* [LotuS (v1.30)](http://psbweb05.psb.ugent.be/lotus/): extremely fast
  OTU building, annotation, phylogeny and abundance matrix pipeline,
  based on raw sequencer output.
* [QIIME (v1.9)](http://qiime.org/): a multi-purpose pipeline for
  performing microbiome analysis from raw DNA sequencing data.


## Alternatives ##

If you want to try alternative free and open-source clustering
methods, here are some links:

* [VSEARCH](https://github.com/torognes/vsearch)
* [Oligotyping](http://merenlab.org/projects/oligotyping/)
* [DNAclust](http://dnaclust.sourceforge.net/)
* [Sumaclust](http://metabarcoding.org/sumatra)
* [Crunchclust](https://code.google.com/p/crunchclust/)


## Version history ##

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
`-c`). In addition, an option (`-w`) allows to output OTU
representatives sequences with updated abundances (sum of all
abundances inside each OTU). This version also enables dereplication
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
algorithm and the built-in OTU breaking by default. Some options are
changed and some new output options are introduced.

### version 1.2.21 ###

**swarm** 1.2.21 is supposed to fix some problems related to the use of the
SSSE3 CPU instructions which are not always available.

### version 1.2.20 ###

**swarm** 1.2.20 presents a production-ready version of the
alternative algorithm (option `-a`), with optional built-in OTU
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
