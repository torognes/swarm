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

#define SEPCHAR ' '
#define BITS 8

unsigned long count_comparisons_8;
unsigned long count_comparisons_16;
unsigned long count_comparisons_again;

unsigned long targetcount;
unsigned long * targets;
unsigned long * scores;
unsigned long * diffs;
unsigned long * alignlengths;

unsigned long targetcount2;
unsigned long * targets2;
unsigned long * scores2;
unsigned long * diffs2;
unsigned long * alignlengths2;

unsigned long * swarmids;
unsigned long * seeded;
unsigned long * genids;


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

void search_again(unsigned long seed)
{
  /* find targets with saturated score */
  targetcount2 = 0;
  for(unsigned long i=0; i<targetcount; i++)
    if (scores[i] == 255)
      targets2[targetcount2++] = targets[i];

  /* redo search */
  search_do(seed, targetcount2, targets2, scores2, diffs2, alignlengths2, 16);
  count_comparisons_again += targetcount2;

  /* replace scores and diffs in original lists */
  unsigned long i = 0;
  for(unsigned long j=0; j<targetcount2; j++)
    {
      /* find corresponding target in original list */
      while (targets[i] != targets2[j])
	i++;

      /* replace scores and diffs */
      scores[i] = scores2[j];
      diffs[i] = diffs2[j];
      alignlengths[i] = alignlengths[2];
    }
}

long mindiff2(unsigned long a, unsigned long b)
{
  long a_len = db_getsequencelen(a);
  long b_len = db_getsequencelen(b);

  long * a_composition = db_getseqinfo(a)->composition;
  long * b_composition = db_getseqinfo(b)->composition;

  long mindiff2 = labs(a_len - b_len);
  for(int i=0; i<4; i++)
    mindiff2 += labs(a_composition[i] - b_composition[i]);

  return mindiff2;
}


void algo_run()
{
  count_comparisons_8 = 0;
  count_comparisons_16 = 0;
  count_comparisons_again = 0;
  
  unsigned long largestswarm = 0;
  unsigned long swarmsize = 0;

  unsigned long listlength = db_getsequencecount();
  
  targets = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));
  scores = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));
  diffs = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));
  alignlengths = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));

  targets2 = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));
  scores2 = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));
  diffs2 = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));
  alignlengths2 = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));

  swarmids = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));
  seeded = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));

  unsigned long * hits = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));

  unsigned long diff_saturation = MIN(255 / penalty_mismatch,
				      255 / (penalty_gapopen + penalty_gapextend));

  /* reset swarmid to zero for all */
  for(unsigned long i=0; i<listlength; i++)
    swarmids[i] = 0;
  
  /* reset seeded to zero for all */
  for(unsigned long i=0; i<listlength; i++)
    seeded[i] = 0;
  
  /* generation id's = distance, in generations, from inital seed */
  genids = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));

  /* radius = actual number of differences from initial seed */
  unsigned long * radius = 
    (unsigned long *) xmalloc(listlength * sizeof(unsigned long));

  unsigned long swarmid = 1;
  unsigned long maxgen = 0;

  unsigned long bits;
  unsigned long bits2;

  
  while (1)
    {

      unsigned long amplicons_copies = 0;
      unsigned long singletons = 0;
      unsigned long hitcount = 0;
      unsigned long diffsum = 0;
      unsigned long alignlensum = 0;
      unsigned long maxradius = 0;

      /* always start in 8 bit mode unless resolution is very high*/

      if ((unsigned long)resolution <= diff_saturation)
	bits = 8;
      else
	bits = 16;

      /* find first non-swarmed sequence */
      
      long seed = -1;
      for(unsigned long i=0; i<listlength; i++)
	if (swarmids[i] == 0)
	  {
	    seed = i;
	    break;
	  }
      
      if (seed < 0)
	break;

      unsigned long generation = 0;

      hits[hitcount++] = seed;

      swarmids[seed] = swarmid;
      seeded[seed] = 1;
      genids[seed] = 0;
      radius[seed] = 0;
      swarmsize = 1;
      
      unsigned long abundance = db_getabundance(seed);
      amplicons_copies += abundance;
      if (abundance == 1)
	singletons++;

      unsigned long lastgen = 0;

      /* find all non-swarmed sequences , prepare for search */
      /* start at seed + 1 */

      targetcount = 0;
      for(unsigned long t=seed+1; t<listlength; t++)
	if (swarmids[t] == 0)
	  targets[targetcount++] = t;

      if (targetcount > 0)
      {
	search_do(seed, targetcount, targets, scores, diffs, alignlengths, bits);

	if (bits == 8)
	  count_comparisons_8 += targetcount;
	else
	  count_comparisons_16 += targetcount;

	for(unsigned long i=0; i<targetcount; i++)
	{
	  if (diffs[i] <= (unsigned long) resolution)
	  {
	    unsigned long t = targets[i];
	    swarmids[t] = swarmid;
	    genids[t] = generation + 1;
	    radius[t] = radius[seed] + diffs[i];
	    if (radius[t] > maxradius)
	      maxradius = radius[t];

	    hits[hitcount++] = t;
	    diffsum += diffs[i];
	    alignlensum += alignlengths[i];
	    
	    abundance = db_getabundance(t);
	    amplicons_copies += abundance;
	    if (abundance == 1)
	      singletons++;

	    lastgen = generation + 1;
	    swarmsize++;
	  }
	}  

	generation = 1;
	
	while (generation <= lastgen)
	{
	
	  if ((unsigned long)resolution <= diff_saturation)
	    bits2 = 8;
	  else
	    bits2 = 16;
	  

	  if (generation > maxgen)
	    maxgen = generation;

	  /* find first potential seed with genid = generation, swarmids=swarmid,
	     seeded = 0 */

	  long subseed = -1;

	  for(unsigned long i=0; i<targetcount; i++)
	  {
	    unsigned long t = targets[i];
	    if ((swarmids[t] == swarmid) && 
		(genids[t] == generation) && 
		(seeded[t] == 0))
	    {
	      subseed = t;
	      break;
	    }
	  }
	

	  if (subseed >= 0)
	  {
	    
	    seeded[subseed] = 1;
	    
	    /* switch to 16 bits if necessary */

	    if ((bits == 8) && 
		(2*generation*resolution > diff_saturation))
	      {
		bits = 16;
		search_again(seed);
	      }

	    /* prepare for new subsearch */
	    
	    targetcount2 = 0;
	    for(unsigned long i=0; i<targetcount; i++)
	    {
	      unsigned long t = targets[i];
	      if ((swarmids[t] == 0)
		  && (diffs[i] <= 2*generation*resolution) 
		  && (mindiff2(t, subseed) <= 2*resolution)
		  )
		targets2[targetcount2++] = t;
	    }
	    
	    if (targetcount2 > 0)
	    {
	      search_do(subseed, targetcount2, targets2, scores2, diffs2, alignlengths2, bits2);
	      if (bits2 == 8)
		count_comparisons_8 += targetcount2;
	      else
		count_comparisons_16 += targetcount2;
	      
	      for(unsigned long i=0; i<targetcount2; i++)
	      {
		if (diffs2[i] <= (unsigned long) resolution)
		{
		  unsigned long t = targets2[i];
		  swarmids[t] = swarmid;
		  genids[t] = generation + 1;
		  radius[t] = radius[subseed] + diffs2[i];
		  if (radius[t] > maxradius)
		    maxradius = radius[t];

		  hits[hitcount++] = t;
		  diffsum += diffs2[i];
		  alignlensum += alignlengths2[i];

		  abundance = db_getabundance(t);
		  amplicons_copies += abundance;
		  if (abundance == 1)
		    singletons++;

		  lastgen = generation + 1;
		  swarmsize++;
		}
	      }  
	    }
	  }
	  else
	  {
	    generation++;
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
	  fprint_id(uclustfile, seed);
	  fprintf(uclustfile, "\t*\n");
	  
	  fprintf(uclustfile, "S\t%lu\t%lu\t*\t*\t*\t*\t*\t",
		  swarmid-1, db_getsequencelen(seed));
	  fprint_id(uclustfile, seed);
	  fprintf(uclustfile, "\t*\n");
	  fflush(uclustfile);

	  for(unsigned long i=1; i<hitcount; i++)
	    {
	      unsigned long hit = hits[i];
	      
	      char * dseq = db_getsequence(hit);
	      char * dend = dseq + db_getsequencelen(hit);
	      char * qseq = db_getsequence(seed);
	      char * qend = qseq + db_getsequencelen(seed);

	      unsigned long nwscore = 0;
	      unsigned long nwdiff = 0;
	      char * nwalignment = NULL;
	      unsigned long nwalignmentlength = 0;

	      nw(dseq, dend, qseq, qend, 
		 score_matrix_63, gapopen, gapextend, 
		 & nwscore, & nwdiff, & nwalignmentlength, & nwalignment);
	      
	      double percentid = 100.0 * (nwalignmentlength - nwdiff) / nwalignmentlength;
	      
	      fprintf(uclustfile, "H\t%lu\t%lu\t%.1f\t+\t0\t0\t%s\t",
		      swarmid-1, db_getsequencelen(hit), percentid, nwalignment);

	      fprint_id(uclustfile, hit);
	      fprintf(uclustfile, "\t");
	      fprint_id(uclustfile, seed);
	      fprintf(uclustfile, "\n");
	      fflush(uclustfile);

	      if (nwalignment)
		free(nwalignment);
	    }

	}
      

      if (statsfile)
	{
	  unsigned long maxgenerations = generation > 0 ? generation - 1 : 0;
	  abundance = db_getabundance(seed);

	  fprintf(statsfile, "%lu\t%lu\t", swarmsize, amplicons_copies);
	  fprint_id_noabundance(statsfile, seed);
	  fprintf(statsfile, "\t%lu\t%lu\t%lu\t%lu\n", 
		  abundance, singletons, maxgenerations, maxradius);
	}


      swarmid++;
    }
  
  free(hits);

  free(genids);
  free(seeded);

  free(swarmids);

  free(targets2);
  free(scores2);
  free(diffs2);
  free(alignlengths2);

  free(targets);
  free(scores);
  free(diffs);
  free(alignlengths);

  fprintf(stderr, "\n");

  fprintf(stderr, "Number of swarms:  %lu\n", swarmid-1);

  fprintf(stderr, "Largest swarm:     %lu\n", largestswarm);

  fprintf(stderr, "Largest radius:    %lu\n", maxgen);

  fprintf(stderr, "\n");

  fprintf(stderr, "Comparisons (8b):  %lu (%.2lf%%)\n", count_comparisons_8, (200.0 * count_comparisons_8 / listlength / listlength));

  fprintf(stderr, "Comparisons (16b): %lu (%.2lf%%)\n", count_comparisons_16, (200.0 * count_comparisons_16 / listlength / listlength));

  fprintf(stderr, "Comparisons (cnv): %lu (%.2lf%%)\n", count_comparisons_again, (200.0 * count_comparisons_again / listlength / listlength));
  
  fprintf(stderr, "Comparisons (tot): %lu (%.2lf%%)\n", count_comparisons_8 + count_comparisons_16 + count_comparisons_again, (200.0 * (count_comparisons_8 + count_comparisons_16 + count_comparisons_again) / listlength / listlength));
  
}

