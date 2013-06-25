# swarm #

A robust and fast clustering method for amplicon-based studies.

The purpose of **swarm** is to provide a novel clustering algorithm to handle large sets of amplicons. Traditional clustering algorithms results are strongly input-order dependent, and rely on an arbitrary **global** clustering threshold. **swarm** results are resilient to input-order changes and rely on a small **local** linking threshold *d*, the maximum number of differences between two amplicons. **swarm** forms stable high-resolution clusters, with a high yield of biological information.

Table of Content
================

* [Quick start](#quick_start)
* [Install](#install)
* [Prepare amplicon fasta files](#prepare_amplicon)
   * [Linearization](#linearization)
   * [Dereplication](#dereplication)
   * [Launch swarm](#launch)
* [Parse swarm results](#parse)
   * [Count the number of amplicons per OTU](#OTU_sizes)
   * [Get the seed sequence for each swarm](#extract_seeds)
   * [Get fasta sequences for all amplicons in a swarm](#extract_all)
* [Troubleshooting](#troubleshooting)
* [New features](#features)
   * [version 1.2.2](#version122)
   * [version 1.2.1](#version121)
   * [version 1.2.0](#version120)
   * [version 1.1.1](#version111)
   * [version 1.1.0](#version110)
       * [Statistics](#stats)
       * [Uclust-like output format](#uclust)

<a name="quick_start"/>
## Quick start ##

**swarm** most simple usage is (with default parameters, use `-h` to get help or see the user manual for details):

```
./swarm amplicons.fasta
```

<a name="install"/>
## Install ##

Get the source code and a **swarm** binary from [GitHub](https://github.com/torognes/swarm "swarm public repository") using the [ZIP button](https://github.com/torognes/swarm/archive/master.zip "swarm zipped folder") or git:

```
git clone https://github.com/torognes/swarm.git
cd swarm/
```

Use the command `make` to compile **swarm** from scratch.

If you have administrator privileges, you can make **swarm** accessible for all users. Simply copy the binary to `/usr/bin/`. The man page can be installed this way:

```
gzip -c swarm.1 > swarm.1.gz
mv swarm.1.gz /usr/share/man/man1/
```

Once installed, the man page can be accessed with the command `man swarm`.

<a name="prepare_amplicon"/>
## Prepare amplicon fasta files ##

To facilitate the use of **swarm**, we provide examples of shell commands that can be use to format and check the input fasta file (warning, this may not be suitable for very large files). The amplicon clipping step (adaptor and primer removal) and the filtering step are not discussed here.

<a name="linearization"/>
### Linearization ###

Amplicons written on two lines are easier to manipulate: one line for the fasta header, one line for the sequence (tested with GNU Awk 4.0.1).

```
awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {print}' amplicons.fasta > amplicons_linearized.fasta
```

<a name="dereplication"/>
### Dereplication ###

To speed up the clustering process, strictly identical amplicons should be merged. This step is not mandatory, but it is an important time saver, especially for highly redundant high-throughput sequencing results.

```
grep -v "^>" amplicons_linearized.fasta | sort -d | uniq -c |
while read abundance sequence ; do
    hash=$(sha1sum <<< "${sequence}")
    hash=${hash:0:40}
    printf ">%s_%d_%s\n" "${hash}" "${abundance}" "${sequence}"
done | sort -t "_" -k2,2nr -k1.2,1d |
sed -e 's/\_/\n/2' > amplicons_linearized_dereplicated.fasta
```

The dereplicated amplicons receive a meaningful unique name (hash codes), and are sorted by decreasing number of copies and by hash codes (to guarantee a stable sorting). The use of a hashing function also provides an easy way to compare sets of amplicons. If two amplicons from two different sets have the same hash code, it means that the sequences they represent are identical.

<a name="launch"/>
### Launch swarm ###

If you want **swarm** to partition your dataset with the finest resolution (a local number of differences *d* = 1) on a quadricore CPU:

```
./swarm -d 1 -t 4 amplicons.fasta > amplicons.swarms
```

See the user manual (man page and PDF) for details on **swarm**'s options and parameters.

<a name="parse"/>
## Parse swarm results ##

To facilitate the use of **swarm**, we provide examples of shell commands that can be use to parse **swarm**'s output. We assume that the amplicon fasta file was prepared as describe above (linearization and dereplication).

<a name="OTU_sizes"/>
### Count the number of amplicons per OTU ###

You might want to check the size distribution of OTU (number of amplicons in each OTU), and count the number of singletons (OTUs containing only one amplicon).

```
awk '{print NF}' amplicons.swarms | sort -n | uniq -c
awk 'NF == 1 {sum+=1} END {print sum}' amplicons.swarms
```

The number of amplicons in each OTU and several other metrics are available in the statistics file produced by swarm when using the -s option.

<a name="extract_seeds"/>
### Get the seed sequence for each OTU ###

It is frequent for subsequent analyses to keep only one representative amplicon per OTU (usually the seed) to reduce the computational burden. That operation is easily done with **swarm** results.

```
SEEDS=$(mktemp)
cut -d " " -f 1 amplicons.swarms | sed -e 's/^/>/' > "${SEEDS}"
grep -A 1 -F -f "${SEEDS}" amplicons.fasta | sed -e '/^--$/d' > amplicons_seeds.fasta
rm "${SEEDS}"
```

<a name="extract_all"/>
### Get fasta sequences for all amplicons in a OTU ###

For each OTU, get the fasta sequences for all amplicons. Warning, this loop can generate a very large number of files. To limit the number of files, a test can be added to exclude swarms with less than *n* elements.

```
INPUT_SWARM="amplicons.swarms"
INPUT_FASTA="amplicons.fasta"
OUTPUT_FOLDER="swarms_fasta"
AMPLICONS=$(mktemp)
mkdir "${OUTPUT_FOLDER}"
while read swarm ; do
    tr " " "\n" <<< "${swarm}" | sed -e 's/^/>/' > "${AMPLICONS}"
    seed=$(head -n 1 "${AMPLICONS}")
    grep -A 1 -F -f "${AMPLICONS}" "${INPUT_FASTA}" | sed -e '/^--$/d' > "./${OUTPUT_FOLDER}/${seed/>/}.fasta"
done < "${INPUT_SWARM}"
rm "${AMPLICONS}"
```

<a name="troubleshooting"/>
## Troubleshooting ##

If **swarm** exits with an error message saying `Error: Illegal character in sequence`, at least one of your amplicon sequence contains a character other than ACGT (or acgt). This command will help you to find the concerned sequences and their line numbers:

```
awk '!/^>/ && /[^ACGTacgt]/ {print NR, $0}' amplicons.fasta
```

If **swarm** exists with an error message saying `This program requires a processor with SSE2`, your computer is too old to run **swarm** (or based on a non x86-64 architecture). **swarm** only runs on CPUs with the SSE2 and POPCNT instructions. The POPCNT instruction is the most recently introduced: November 2008 for Intel CPUs and October 2011 for AMD CPUs. **swarm** should be able to run on most Intel or AMD CPU released since.

<a name="features"/>
## New features##

<a name="version122"/>
### version 1.2.2 ###

**swarm** 1.2.2 fixes an issue with incorrect values in the statistics file (maximum generation and radius of swarms). This version is also a bit faster.

<a name="version121"/>
### version 1.2.1 ###

**swarm** 1.2.1 removes the need for a SSE4.1 capable CPU and should now be able to run on most servers, desktops and laptops.

<a name="version120"/>
### version 1.2.0 ###

**swarm** 1.2.0 introduces a pre-filtering of similar amplicons based on *k*-mers. This eliminates most of the time-consuming pairwise alignments and greatly improves speed. The speedup can be more than 100-fold compared to previous swarm versions when using a single thread with a large set of amplicons. Using multiple threads induces a computational overhead, but becomes more and more efficient as the size of the amplicon set increases.

<a name="version111"/>
### version 1.1.1 ###

**swarm** now works on Apple computers. This version also corrects an issue in the pairwise global alignment step that could lead to sub-optimal alignments. Slightly different alignments may result relative to previous version, giving slightly different swarms.

<a name="version110"/>
### version 1.1.0 ###

**swarm** 1.1.0 introduces new optimizations and is 20% faster than the previous version on our test dataset. It also introduces two new output options: statistics and uclust-like format.

<a name="stats"/>
#### Statistics ####

By specifying the `-s` option to **swarm** it will now output detailed statistics about each swarm to a specified file. It will print the number of unique amplicons, the number of copies, the name of the seed and its abundance, the number of singletons (amplicons with an abundance of 1), the number of iterations and the maximum radius of the swarm (i.e. number of differences between the seed and the furthermost amplicon). When using input data sorted by decreasing abundance, the seed is the most abundant amplicon in the swarm.

<a name="uclust"/>
#### Uclust-like output format ####

Some pipelines use the [uclust output format](http://www.drive5.com/uclust/uclust_userguide_1_1_579.html#_Toc257997686 "page describing the uclust output format") as input for subsequent analyses. **swarm** can now output results in this format to a specified file with the `-u` option.
