/*
    SWARM

    Copyright (C) 2012 Torbjorn Rognes and Frederic Mahe

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
//#define SEPCHAR '\t'
#define BITS 8

unsigned long count_comparisons_8;
unsigned long count_comparisons_16;
unsigned long count_comparisons_again;

unsigned long targetcount;
unsigned long * targets;
unsigned long * scores;
unsigned long * diffs;

unsigned long targetcount2;
unsigned long * targets2;
unsigned long * scores2;
unsigned long * diffs2;

unsigned long * swarmids;
unsigned long * seeded;
unsigned long * genids;

void search_again(unsigned long seed)
{
  /* find targets with saturated score */
  targetcount2 = 0;
  for(unsigned long i=0; i<targetcount; i++)
    if (scores[i] == 255)
      targets2[targetcount2++] = targets[i];

  /* redo search */
  search_do(seed, targetcount2, targets2, scores2, diffs2, 16);
  count_comparisons_again += targetcount2;

  /* replace scores and diffs in original lists */
  unsigned long i = 0;
  for(unsigned long j=0; j<targetcount2; j++)
    {
      /* find corresponding target in original list */
      while (targets[i] != targets2[j])
	i++;
      /* replace scores and diffs */

#if 0
      fprintf(stderr, "Target: %lu  Score: %lu -> %lu  Diff: %lu -> %lu\n",
	      targets[i], scores[i], scores2[j], diffs[i], diffs2[j]);
#endif

      scores[i] = scores2[j];
      diffs[i] = diffs2[j];
    }
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

  targets2 = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));
  scores2 = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));
  diffs2 = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));

  swarmids = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));
  seeded = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));

  unsigned long diff_saturation = MIN(255 / penalty_mismatch,
				      255 / (penalty_gapopen + penalty_gapextend));

  /* reset swarmid to zero for all */
  for(unsigned long i=0; i<listlength; i++)
    swarmids[i] = 0;
  
  /* reset seeded to zero for all */
  for(unsigned long i=0; i<listlength; i++)
    seeded[i] = 0;
  
  /* generation id's = distance from inital seed */
  genids = (unsigned long *) xmalloc(listlength * sizeof(unsigned long));

  unsigned long swarmid = 1;

  unsigned long maxgen = 0;

  unsigned long bits;
  unsigned long bits2;

  while (1)
    {

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

      //      printf("Swarm %4lu  Generation %4lu  Hit %4lu  %s (Seed)\n", swarmid, generation, seed, db_getheader(seed));
      
      fputs(db_getheader(seed), out);

      swarmids[seed] = swarmid;
      seeded[seed] = 1;
      genids[seed] = 0;
      swarmsize = 1;
      
      unsigned long lastgen = 0;

      /* find all non-swarmed sequences , prepare for search */
      /* start at seed + 1 */

      targetcount = 0;
      for(unsigned long t=seed+1; t<listlength; t++)
	if (swarmids[t] == 0)
	  targets[targetcount++] = t;

      if (targetcount > 0)
      {
	search_do(seed, targetcount, targets, scores, diffs, bits);
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
	    //	      printf("Swarm %4lu  Generation %4lu  Hit %4lu  %s\n", swarmid, generation, t, db_getheader(t));
	    fputc(SEPCHAR, out);
	    fputs(db_getheader(t), out);
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
	      if ((swarmids[t] == 0) && (diffs[i] <= 2*generation*resolution))
		targets2[targetcount2++] = t;
	    }
	    
	    if (targetcount2 > 0)
	    {
	      search_do(subseed, targetcount2, targets2, scores2, diffs2, bits2);
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
		  // printf("Swarm %4lu  Generation %4lu  Hit %4lu  %s\n", swarmid, generation, t, db_getheader(t));
		  fputc(SEPCHAR, out);
		  fputs(db_getheader(t), out);
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
      
      fputs("\n", out);
      
      if (swarmsize > largestswarm)
	largestswarm = swarmsize;

      swarmid++;
    }
  
  free(genids);
  free(seeded);

  free(swarmids);

  free(targets2);
  free(scores2);
  free(diffs2);

  free(targets);
  free(scores);
  free(diffs);

  fprintf(stderr, "\n");

  fprintf(stderr, "Number of swarms:  %lu\n", swarmid-1);

  fprintf(stderr, "Largest swarm:     %lu\n", largestswarm);

  fprintf(stderr, "Generation depth:  %lu\n", maxgen);

  fprintf(stderr, "Comparisons (8b):  %lu (%.2lf%%)\n", count_comparisons_8, (200.0 * count_comparisons_8 / listlength / listlength));

  fprintf(stderr, "Comparisons (16b): %lu (%.2lf%%)\n", count_comparisons_16, (200.0 * count_comparisons_16 / listlength / listlength));

  fprintf(stderr, "Comparisons (cnv): %lu (%.2lf%%)\n", count_comparisons_again, (200.0 * count_comparisons_again / listlength / listlength));
  
  fprintf(stderr, "Comparisons (tot): %lu (%.2lf%%)\n", count_comparisons_8 + count_comparisons_16 + count_comparisons_again, (200.0 * (count_comparisons_8 + count_comparisons_16 + count_comparisons_again) / listlength / listlength));
  
}

