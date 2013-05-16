/*
    SWARM

    Copyright (C) 2012-2013 Torbjorn Rognes and Frederic Mahe

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

#include "swarm.h"

//#define DEBUG

#define SEPCHAR ' '
#define BITS 8

unsigned long count_comparisons_8;
unsigned long count_comparisons_16;
unsigned long count_comparisons_again;

unsigned long targetcount;
unsigned long * targetindices;
unsigned long * targetampliconids;
unsigned long * scores;
unsigned long * diffs;
unsigned long * alignlengths;
unsigned long * qgramamps;
unsigned long * qgramdiffs;
unsigned long * qgramindices;

struct ampliconinfo_s
{
  unsigned ampliconid;
  unsigned diffestimate; /* a lower bound estimate on distance from initial seed */
  unsigned swarmid;
  unsigned generation;
  unsigned radius; /* actual diff from initial seed */
} * amps;

unsigned swarmed;
unsigned seeded;

void fprint_id(FILE * stream, unsigned long x)
{
  char * h = db_getheader(x);
  char c = *h++;
  while (c && (c != ' '))
    {
      fputc(c, stream);
      c = *h++;
    }
}

void fprint_id_noabundance(FILE * stream, unsigned long x)
{
  char * h = db_getheader(x);
  char c = *h++;
  while (c && (c != ' ') && (c != '_'))
    {
      fputc(c, stream);
      c = *h++;
    }
}

void algo_run()
{
  count_comparisons_8 = 0;
  count_comparisons_16 = 0;
  count_comparisons_again = 0;

  unsigned long searches = 0;
  unsigned long estimates = 0;

  unsigned long largestswarm = 0;
  unsigned long swarmsize = 0;

  unsigned long amplicons = db_getsequencecount();
  unsigned long longestamplicon = db_getlongestsequence();
  
  amps = (ampliconinfo_s *) xmalloc(amplicons * sizeof(struct ampliconinfo_s));

  targetampliconids = (unsigned long *) xmalloc(amplicons * sizeof(unsigned long));
  targetindices = (unsigned long *) xmalloc(amplicons * sizeof(unsigned long));
  scores = (unsigned long *) xmalloc(amplicons * sizeof(unsigned long));
  diffs = (unsigned long *) xmalloc(amplicons * sizeof(unsigned long));
  alignlengths = (unsigned long *) xmalloc(amplicons * sizeof(unsigned long));

  qgramamps = (unsigned long *) xmalloc(amplicons * sizeof(unsigned long));
  qgramdiffs = (unsigned long *) xmalloc(amplicons * sizeof(unsigned long));
  qgramindices = (unsigned long *) xmalloc(amplicons * sizeof(unsigned long));

  unsigned long * hits = (unsigned long *) xmalloc(amplicons * sizeof(unsigned long));

  unsigned long diff_saturation = MIN(255 / penalty_mismatch,
				      255 / (penalty_gapopen + penalty_gapextend));

  unsigned char * dir = 0;
  unsigned long * hearray = 0;

  if (uclustfile)
  {
    dir = (unsigned char *) xmalloc(longestamplicon*longestamplicon);
    hearray = (unsigned long *) xmalloc(2 * longestamplicon * sizeof(unsigned long));
  }

  /* set ampliconid for all */
  for(unsigned long i=0; i<amplicons; i++)
  {
    amps[i].ampliconid = i;
  }


  /* always search in 8 bit mode unless resolution is very high */
  
  unsigned long bits;

  if ((unsigned long)resolution <= diff_saturation)
    bits = 8;
  else
    bits = 16;
 
  seeded = 0;
  swarmed = 0;

  unsigned long swarmid = 0;
  unsigned long maxgen = 0;

  while (seeded < amplicons)
    {

      /* process each initial seed */
      
      swarmid++;

      unsigned long amplicons_copies = 0;
      unsigned long singletons = 0;
      unsigned long hitcount = 0;
      unsigned long diffsum = 0;
      unsigned long maxradius = 0;
      unsigned long generation = 1;
      unsigned long seedindex;

      seedindex = seeded;
      seeded++;

      amps[seedindex].swarmid = swarmid;
      amps[seedindex].generation = generation;
      amps[seedindex].radius = 0;
     
      unsigned long seedampliconid = amps[seedindex].ampliconid;
      hits[hitcount++] = seedampliconid;
      
#ifdef DEBUG
      fprintf(stderr, "Seed amplicon id: %lu\n", seedampliconid);
#endif

      unsigned long abundance = db_getabundance(seedampliconid);
      amplicons_copies += abundance;
      if (abundance == 1)
	singletons++;

      swarmsize = 1;
      swarmed++;


      /* find diff estimates between seed and each amplicon in pool */

      targetcount = 0;

      unsigned long listlen = amplicons - swarmed;

      for(unsigned long i=0; i<listlen; i++)
	qgramamps[i] = amps[swarmed+i].ampliconid;
      qgram_diff_parallel(seedampliconid, listlen, qgramamps, qgramdiffs);
      estimates += listlen;
      
      for(unsigned long i=0; i < listlen; i++)
	{
	  unsigned poolampliconid = qgramamps[i];
	  long diff = qgramdiffs[i];
	  amps[swarmed+i].diffestimate = diff;
	  if (diff <= resolution)
	    {
	      targetindices[targetcount] = swarmed+i;
	      targetampliconids[targetcount] = poolampliconid;
	      targetcount++;
	    }
	}
      
      if (targetcount > 0)
      {
	search_do(seedampliconid, targetcount, targetampliconids, 
		  scores, diffs, alignlengths, bits);
	searches++;

	if (bits == 8)
	  count_comparisons_8 += targetcount;
	else
	  count_comparisons_16 += targetcount;

	for(unsigned long t=0; t<targetcount; t++)
	{
	  unsigned diff = diffs[t];

	  if (diff <= (unsigned long) resolution)
	  {
	    unsigned i = targetindices[t];

	    /* move the target (i) to the position (swarmed) of the first */
	    /* unswarmed amplicon in the pool */
	    
#ifdef DEBUG
	    unsigned long targetampliconid = amps[i].ampliconid;

	    fprintf(stderr,
		    "Moving from pos %u to %u. Seeded: %u Swarmed: %u ta: %lu pa: %u\n",
		    i, swarmed, seeded, swarmed, targetampliconid, amps[swarmed].ampliconid);
#endif

	    if (swarmed < i)
	    {
	      struct ampliconinfo_s temp = amps[i];
	      for(unsigned j=i; j>swarmed; j--)
	      {
		amps[j] = amps[j-1];
	      }
	      amps[swarmed] = temp;
	    }

	    amps[swarmed].swarmid = swarmid;
	    amps[swarmed].generation = generation;
	    amps[swarmed].radius = diff;
	    if (diff > maxradius)
	      maxradius = diff;

	    unsigned poolampliconid = amps[swarmed].ampliconid;
	    hits[hitcount++] = poolampliconid;

	    diffsum += diff;
	    
	    abundance = db_getabundance(poolampliconid);
	    amplicons_copies += abundance;
	    if (abundance == 1)
	      singletons++;

	    swarmsize++;

	    swarmed++;
	  }
	}  


	while (seeded < swarmed)
	{

	  /* process each subseed */

	  unsigned subseedampliconid;
	  unsigned subseedradius;
	  
	  unsigned long subseedindex;
	  
	  subseedindex = seeded;
	  subseedampliconid = amps[subseedindex].ampliconid;
	  subseedradius = amps[subseedindex].radius;
	  generation = amps[subseedindex].generation;
	  if (generation > maxgen)
	    maxgen = generation;

	  seeded++;
	  
#ifdef DEBUG
	  fprintf(stderr, "Subseed amplicon id: %u Generation: %lu\n", subseedampliconid, generation);
#endif

	  targetcount = 0;
	  
	  unsigned long listlen=0;
	  for(unsigned long i=swarmed; i<amplicons; i++)
	    {
	      unsigned long targetampliconid = amps[i].ampliconid;
	      if (amps[i].diffestimate <= 2*generation*resolution)
		{
		  qgramamps[listlen] = targetampliconid;
		  qgramindices[listlen] = i;
		  listlen++;
		}
	    }
      	  qgram_diff_parallel(subseedampliconid, listlen, qgramamps, qgramdiffs);
	  estimates += listlen;
      
	  for(unsigned long i=0; i < listlen; i++)
	    if ((long)qgramdiffs[i] <= resolution)
	      {
		targetindices[targetcount] = qgramindices[i];
		targetampliconids[targetcount] = qgramamps[i];
		targetcount++;
	      }
	  
	  if (targetcount > 0)
	  {
	    search_do(subseedampliconid, targetcount, targetampliconids, 
		      scores, diffs, alignlengths, bits);
	    searches++;

	    if (bits == 8)
	      count_comparisons_8 += targetcount;
	    else
	      count_comparisons_16 += targetcount;
	    
	    for(unsigned long t=0; t<targetcount; t++)
	    {
	      unsigned diff = diffs[t];
	      
	      if (diff <= (unsigned long) resolution)
	      {
		unsigned i = targetindices[t];
		
		/* move the target (i) to the position (swarmed) of the first */
		/* unswarmed amplicon in the pool */
		
		/* then move the target further into the swarmed but unseeded part */
		/* of the list, so that the swarmed amplicons are ordered by id */

		/* find correct position in list */

		unsigned long targetampliconid = amps[i].ampliconid;
		unsigned pos = swarmed;

		while ((pos > seeded) && (amps[pos-1].ampliconid > targetampliconid) &&
		       (amps[pos-1].generation > generation))
		  pos--;

#ifdef DEBUG
		fprintf(stderr,
			"Moving from pos %u to %u. Seeded: %u Swarmed: %u ta: %lu pa: %u\n",
			i, pos, seeded, swarmed, targetampliconid, amps[pos].ampliconid);
#endif

		if (pos < i)
		{
		  struct ampliconinfo_s temp = amps[i];
		  for(unsigned j=i; j>pos; j--)
		    {
		      amps[j] = amps[j-1];
		    }
		  amps[pos] = temp;
		}

		amps[pos].swarmid = swarmid;
		amps[pos].generation = generation + 1;
		amps[pos].radius = subseedradius + diff;
		if (diff > maxradius)
		  maxradius = diff;

		unsigned poolampliconid = amps[pos].ampliconid;
		hits[hitcount++] = poolampliconid;
		diffsum += diff;

		abundance = db_getabundance(poolampliconid);
		amplicons_copies += abundance;
		if (abundance == 1)
		  singletons++;

		swarmsize++;

		swarmed++;
	      }
	    }  
	  }
	}
      }
      
      if (swarmsize > largestswarm)
	largestswarm = swarmsize;

      for(unsigned long i=0; i<hitcount; i++)
	{
	  if (i>0)
	    fputc(SEPCHAR, outfile);
	  fprint_id(outfile, hits[i]);
	}
      fputs("\n", outfile);

      if (uclustfile)
	{
	  fprintf(uclustfile, "C\t%lu\t%lu\t*\t*\t*\t*\t*\t",
		  swarmid-1, swarmsize);
	  fprint_id(uclustfile, seedampliconid);
	  fprintf(uclustfile, "\t*\n");
	  
	  fprintf(uclustfile, "S\t%lu\t%lu\t*\t*\t*\t*\t*\t",
		  swarmid-1, db_getsequencelen(seedampliconid));
	  fprint_id(uclustfile, seedampliconid);
	  fprintf(uclustfile, "\t*\n");
	  fflush(uclustfile);

	  for(unsigned long i=1; i<hitcount; i++)
	    {
	      unsigned long hit = hits[i];
	      
	      char * dseq = db_getsequence(hit);
	      char * dend = dseq + db_getsequencelen(hit);
	      char * qseq = db_getsequence(seedampliconid);
	      char * qend = qseq + db_getsequencelen(seedampliconid);

	      unsigned long nwscore = 0;
	      unsigned long nwdiff = 0;
	      char * nwalignment = NULL;
	      unsigned long nwalignmentlength = 0;

	      nw(dseq, dend, qseq, qend, 
		 score_matrix_63, gapopen, gapextend, 
		 & nwscore, & nwdiff, & nwalignmentlength, & nwalignment,
		 dir, hearray, 0, 0);
	      
	      double percentid = 100.0 * (nwalignmentlength - nwdiff) / nwalignmentlength;
	      
	      fprintf(uclustfile, "H\t%lu\t%lu\t%.1f\t+\t0\t0\t%s\t",
		      swarmid-1, db_getsequencelen(hit), percentid, 
		      nwdiff > 0 ? nwalignment : "=");
	      
	      fprint_id(uclustfile, hit);
	      fprintf(uclustfile, "\t");
	      fprint_id(uclustfile, seedampliconid);
	      fprintf(uclustfile, "\n");
	      fflush(uclustfile);

	      if (nwalignment)
		free(nwalignment);
	    }

	}
      

      if (statsfile)
	{
	  unsigned long maxgenerations = generation > 0 ? generation - 1 : 0;
	  abundance = db_getabundance(seedampliconid);

	  fprintf(statsfile, "%lu\t%lu\t", swarmsize, amplicons_copies);
	  fprint_id_noabundance(statsfile, seedampliconid);
	  fprintf(statsfile, "\t%lu\t%lu\t%lu\t%lu\n", 
		  abundance, singletons, maxgenerations, maxradius);
	}
    }
  

  if (uclustfile)
  {
    free(dir);
    free(hearray);
  }

  free(qgramdiffs);
  free(qgramamps);
  free(qgramindices);
  free(hits);
  free(alignlengths);
  free(diffs);
  free(scores);
  free(targetindices);
  free(targetampliconids);
  free(amps);


  fprintf(stderr, "\n");

  fprintf(stderr, "Number of swarms:  %lu\n", swarmid);

  fprintf(stderr, "Largest swarm:     %lu\n", largestswarm);

  fprintf(stderr, "Max generations:   %lu\n", maxgen);

  fprintf(stderr, "\n");

  fprintf(stderr, "Estimates:         %lu\n", estimates);

  fprintf(stderr, "Searches:          %lu\n", searches);

  fprintf(stderr, "\n");

  fprintf(stderr, "Comparisons (8b):  %lu (%.2lf%%)\n", count_comparisons_8, (200.0 * count_comparisons_8 / amplicons / (amplicons+1)));

  fprintf(stderr, "Comparisons (16b): %lu (%.2lf%%)\n", count_comparisons_16, (200.0 * count_comparisons_16 / amplicons / (amplicons+1)));

  fprintf(stderr, "Comparisons (cnv): %lu (%.2lf%%)\n", count_comparisons_again, (200.0 * count_comparisons_again / amplicons / (amplicons+1)));
  
  fprintf(stderr, "Comparisons (tot): %lu (%.2lf%%)\n", count_comparisons_8 + count_comparisons_16 + count_comparisons_again, (200.0 * (count_comparisons_8 + count_comparisons_16 + count_comparisons_again) / amplicons / (amplicons+1)));

}

