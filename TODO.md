# swarm's future #

As sequencing get cheaper, **swarm** must gain speed to be able to
tackle with ever increasingly larger data sets. In that document, we
list some algorithmic changes that could increase swarm's speed.

Table of Content
================

* [More k-mer based filtering](#kmers)
* [Nucleotidic profile based filtering](#nucleotidic_profile)
* [Use prefix and suffix trees](#prefix-suffix)

<a name="kmers"/>
## More k-mer based filtering ##

In older versions of swarm, most of the computation time was spent on
pairwise global alignments.

In version 1.2 and superior, swarm implements a pre-filtering of
similar amplicons based on k-mers. That pre-filtering relies on the
idea the similarity of two amplicons can be quickly assessed by
comparing their k-mer profiles (i.e. the presence-absence of
nucleotidic strings). If two divergent amplicons can accidentally have
similar k-mer profiles, two similar amplicons are bound to have
similar k-mer profiles. By selecting amplicons with at most *k*
differences in their k-mer profiles, we have a fast and exact way to
avoid many unnecessary and costly pairwise alignments.

The efficiency of that pre-filtering is proportional to the length of
the k-mers. The longer the k-mers, the more efficient the
filtering. However, the cost of storing and comparing the k-mer
profiles increases 4 times more rapidly than *k*.

Swarm uses *k* = 5, and k-mer profile comparisons now represent most
of the computation time. Can we in turn limit the number of k-mer
profile comparisons?

One possibility is to apply a 4-mer pre-pre-filtering to avoid
unnecessary 5-mer vector comparisons. A 4-mer vector comparison
contains less information, but it is 4 times shorter and 4 times
faster to compute. So a speed gain if possible if the 4-mer filtering
step reduces the search-space for the 5-mer filtering by at least 25%.

If that approach turns out to be efficient, why not add a 3-mer
filtering before the 4-mer filtering? The memory consumption would
increase by 30% and the k-mer vector would look like that: 1 x 64 bits
+ 4 x 64 bits + 16 x 16 bits.

The maximal theoretical speed-up is 16x (if 3-mer filtering has the
same efficiency than 5-mer filtering). An observed 2x or 4x speed-up
would be great and would make swarm's speed comparable to usearch's.

<a name="nucleotidic_profile"/>
## Nucleotidic profile based filtering ##

The idea is to compare the nucleotidic composition of two amplicons to
avoid a more costly comparison. The nucleotidic composition of each
amplicon can be stored in four 16-bit long integers, one integer for
each nucleotide. We are searching for amplicons with at most *d*
differences, how does it translates when comparing nucleotidic
compositions?

```
threshold = (d * 2) - length_difference
```

If a comparison yields a value greater than the above threshold, then
there is no need to go further for that specific pair. For a pair of
amplicons *X* and *Y* he value is:

```
abs(Xa -Ya) + abs(Xc -Yc) + abs(Xg -Yg) + abs(Xt -Yt)
```

Adding the amplicon length difference in the mix should improve a bit
the efficiency of the filtering. Tests are needed to verify if
nucleotidic composition comparisons gives a significant speed up.

<a name="prefix-suffix"/>
## Use prefix and suffix trees ##

That approach targets specifically the situation where we are
searching for amplicons with only one difference (`d = 1`).

A prefix tree is a data structure that allows for a fast retrieval of
all amplicons starting with the same sequence (i.e. the same
prefix). A suffix tree does the same for sequence endings.

The algorithm is radically different from the one we are using now,
and would probably consume (much) more memory.

Let's assume that we can construct the prefix-suffix tree in linear
time. The sequences we are storing use a small alphabet (ACGT), and we
can expect a high level of redundancy (and more compact trees).

For a given amplicon ACGTACGT, let's search for amplicons with one
difference. It is an iterative process that will scan each position of
the sequence. Here is the description of what happens for the position 4:

- list all amplicons with the prefix ACGT (listA),
- list all amplicons with the suffix ACGT (listB),
  - if length difference with seed amplicon is 1 and amplicon is in
    listA then we have found all amplicons with an insertion at that point,
- list all amplicons with the suffix CGT (listC),
  - if length difference with seed amplicon is (-1 or 0) and amplicon
    is in listA then we have found all amplicons with a deletion or a
    mutation at that point,

Repeating that process for all position in the sequence (except for
the terminal positions which are special cases) guarantees to find all
micro-variants of a given amplicons. The process creates many lists,
but these lists can be reduced by removing already clustered amplicons
(and the seed). We are searching for amplicons present in both listA
and listB or listA and listC. Near the extremities of our seed
amplicon, the lists can be huge (many amplicons starting with the same
small prefix). The question is: can we manipulate large lists, and can
we intersect lists in a fast and efficient way?
