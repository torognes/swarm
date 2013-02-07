# swarm #

A robust and fast clustering method for amplicon-based studies.

The purpose of **swarm** is to provide a novel clustering algorithm to handle large sets of amplicons. Traditional clustering algorithms results are strongly input-order dependent, and rely on an arbitrary **global** clustering threshold. **swarm** results are resilient to input-order changes and rely on a small **local** linking threshold *d*, the maximum number of differences between two amplicons. **swarm** forms stable high-resolution clusters, with a high yield of biological information.

## Warning ##

**swarm** only runs on CPUs with SSE4.1 instructions. These instructions were introduced by Intel in November 2007 for servers and January 2008 for desktop and portable CPUs. It has been supported by AMD CPUs since October 2011. **swarm** should be able to run on any Intel or AMD CPU released since.

## Quick start ##

**swarm** most simple usage is (with default parameters, see user manual for details):

```
./swarm amplicons.fasta
```

## Install ##

Get the source code and a binary from [GitHub](https://github.com/torognes/swarm "swarm public repository"):

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

## Prepare amplicon fasta files ##

To facilitate the use of **swarm**, we provide examples of shell commands that can be use to format and check the input fasta file (warning, this may not be suitable for very large files). The amplicon clipping step (adaptor and primer removal) and the filtering step are not discussed here.

### linearization ###

Amplicons written on two lines are easier to manipulate (one line for the fasta header, one line for the sequence):

```
awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {print}' amplicons.fasta > amplicons_linearized.fasta
```

### dereplication ###

To speed up the clustering process, strictly identical amplicons should be merged. This step is not mandatory, but it is an important time saver, especially for highly redundant high-throughput sequencing results.

```
grep -v "^>" amplicons_linearized.fasta | sort -d | uniq -c |
while read abundance sequence ; do
    hash=$(echo ${sequence} | sha1sum)
    hash=${hash:0:40}
    printf ">%s_%d_%s\n" "${hash}" "${abundance}" "${sequence}"
done | sort -t "_" -k2,2nr | sed -e 's/\_/\n/2' > amplicons_linearized_dereplicated.fasta
```

The dereplicated amplicons receive a meaningful unique name, and are sorted by decreasing number of copies. The use of a hashing function also provides an easy way to compare sets of amplicons. If two amplicons from two different sets have the same hash, it means that they have identical sequences.

### Search for duplicated amplicon names ###

**swarm** does not check if your amplicons have unique names. To avoid ambiguous results, test your fasta file with this bash command (it should print nothing, unless you have duplicated amplicon names):

```
grep -o -E "^>\S+" amplicons.fasta | tr -d ">" | sort -d | uniq -d
```

See the user manual (man page and PDF) for details on **swarm**'s options and parameters.
