/*
    SWARM

    Copyright (C) 2012-2017 Torbjorn Rognes and Frederic Mahe

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

#define BITS 8

static unsigned long count_comparisons_8;
static unsigned long count_comparisons_16;

static unsigned long targetcount;
static unsigned long * targetindices;
static unsigned long * targetampliconids;
static unsigned long * scores;
static unsigned long * diffs;
static unsigned long * alignlengths;
static unsigned long * qgramamps;
static unsigned long * qgramdiffs;
static unsigned long * qgramindices;

static struct ampliconinfo_s
{
  unsigned ampliconid;
  unsigned diffestimate; /* lower bound estimate of dist from initial seed */
  unsigned swarmid;
  unsigned generation;
  unsigned radius; /* actual diff from initial seed */
} * amps;

static unsigned long swarmed;
static unsigned long seeded;

void algo_run()
{
  count_comparisons_8 = 0;
  count_comparisons_16 = 0;

#ifdef VERBOSE
  unsigned long searches = 0;
  unsigned long estimates = 0;
#endif

  unsigned long largestswarm = 0;

  unsigned long maxgenerations = 0;

  unsigned long amplicons = db_getsequencecount();
  unsigned long longestamplicon = db_getlongestsequence();

  db_qgrams_init();

  qgram_diff_init();

  amps = (struct ampliconinfo_s *) xmalloc(amplicons * sizeof(struct ampliconinfo_s));

  targetampliconids = (unsigned long *) xmalloc(amplicons * 
                                                sizeof(unsigned long));
  targetindices = (unsigned long *) xmalloc(amplicons * sizeof(unsigned long));
  scores = (unsigned long *) xmalloc(amplicons * sizeof(unsigned long));
  diffs = (unsigned long *) xmalloc(amplicons * sizeof(unsigned long));
  alignlengths = (unsigned long *) xmalloc(amplicons * sizeof(unsigned long));

  qgramamps = (unsigned long *) xmalloc(amplicons * sizeof(unsigned long));
  qgramdiffs = (unsigned long *) xmalloc(amplicons * sizeof(unsigned long));
  qgramindices = (unsigned long *) xmalloc(amplicons * sizeof(unsigned long));

  unsigned long * hits = (unsigned long *) xmalloc(amplicons *
                                                   sizeof(unsigned long));

  unsigned long diff_saturation = MIN(255 / penalty_mismatch,
                                      255 / (penalty_gapopen + 
                                             penalty_gapextend));

  unsigned char * dir = 0;
  unsigned long * hearray = 0;

  if (uclustfile)
    {
      dir = (unsigned char *) xmalloc(longestamplicon*longestamplicon);
      hearray = (unsigned long *) xmalloc(2 * longestamplicon *
                                          sizeof(unsigned long));
    }

  /* set ampliconid for all */
  for(unsigned long i=0; i<amplicons; i++)
    {
      amps[i].ampliconid = i;
    }

  /* always search in 8 bit mode unless resolution is very high */
  
  unsigned long bits;

  if ((unsigned long)opt_differences <= diff_saturation)
    bits = 8;
  else
    bits = 16;
 
  seeded = 0;
  swarmed = 0;

  unsigned long swarmid = 0;
  
  progress_init("Clustering:       ", amplicons);
  while (seeded < amplicons)
    {

      /* process each initial seed */
      
      swarmid++;

      unsigned long swarmsize = 0;
      unsigned long amplicons_copies = 0;
      unsigned long singletons = 0;
      unsigned long hitcount = 0;
      unsigned long maxradius = 0;
      unsigned long maxgen = 0;
      unsigned long seedindex;

      seedindex = seeded;
      seeded++;

      amps[seedindex].swarmid = swarmid;
      amps[seedindex].generation = 0;
      amps[seedindex].radius = 0;
     
      unsigned long seedampliconid = amps[seedindex].ampliconid;
      hits[hitcount++] = seedampliconid;
      
      unsigned long abundance = db_getabundance(seedampliconid);
      amplicons_copies += abundance;
      if (abundance == 1)
        singletons++;

      swarmsize = 1;
      swarmed++;


      /* find diff estimates between seed and each amplicon in pool */

      targetcount = 0;

      unsigned long listlen = 0;

      for(unsigned long i=0; i < amplicons-swarmed; i++)
        {
          unsigned ampid = amps[swarmed+i].ampliconid;
          if ((opt_no_otu_breaking) || (db_getabundance(ampid) <= abundance))
            {
              qgramamps[i] = ampid;
              listlen++;
            }
        }

      qgram_diff_fast(seedampliconid, listlen, qgramamps, qgramdiffs);

#ifdef VERBOSE
      estimates += listlen;
#endif
      
      for(unsigned long i=0; i < listlen; i++)
        {
          unsigned poolampliconid = qgramamps[i];
          long diff = qgramdiffs[i];
          amps[swarmed+i].diffestimate = diff;
          if (diff <= opt_differences)
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
#ifdef VERBOSE
          searches++;
#endif

          if (bits == 8)
            count_comparisons_8 += targetcount;
          else
            count_comparisons_16 += targetcount;

          for(unsigned long t=0; t<targetcount; t++)
            {
#if 0
              printf("seed: %lu target: %lu score: %lu "
                     "diffs: %lu alignlen: %lu bits: %lu\n",
                     seedampliconid,
                     targetampliconids[t],
                     scores[t],
                     diffs[t],
                     alignlengths[t],
                     bits);
#endif

              unsigned diff = diffs[t];

              if (diff <= (unsigned long) opt_differences)
                {
                  unsigned i = targetindices[t];

                  /* move the target (i) to the position (swarmed)
                     of the first unswarmed amplicon in the pool */
                  
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
                  amps[swarmed].generation = 1;
                  if (maxgen < 1)
                    maxgen = 1;
                  amps[swarmed].radius = diff;
                  if (diff > maxradius)
                    maxradius = diff;

                  unsigned poolampliconid = amps[swarmed].ampliconid;
                  hits[hitcount++] = poolampliconid;

                  if (opt_internal_structure)
                    {
                      fprint_id_noabundance(internal_structure_file, seedampliconid);
                      fprintf(internal_structure_file, "\t");
                      fprint_id_noabundance(internal_structure_file, poolampliconid);
                      fprintf(internal_structure_file, "\t%u", diff);
                      fprintf(internal_structure_file, "\t%lu\t1", swarmid);
                      fprintf(internal_structure_file, "\n");
                    }

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
              unsigned long subseedgeneration;
              unsigned long subseedabundance;
          
              subseedindex = seeded;
              subseedampliconid = amps[subseedindex].ampliconid;
              subseedradius = amps[subseedindex].radius;
              subseedgeneration = amps[subseedindex].generation;
              subseedabundance = db_getabundance(subseedampliconid);

              seeded++;
          
              targetcount = 0;
          
              unsigned long listlen=0;
              for(unsigned long i=swarmed; i<amplicons; i++)
                {
                  unsigned long targetampliconid = amps[i].ampliconid;
                  if ((amps[i].diffestimate <= subseedradius + opt_differences) &&
                      ((opt_no_otu_breaking) || 
                       (db_getabundance(targetampliconid)
                        <= subseedabundance)))
                    {
                      qgramamps[listlen] = targetampliconid;
                      qgramindices[listlen] = i;
                      listlen++;
                    }
                }

              qgram_diff_fast(subseedampliconid, listlen, qgramamps, 
                              qgramdiffs);

#ifdef VERBOSE
              estimates += listlen;
#endif

              for(unsigned long i=0; i < listlen; i++)
                if ((long)qgramdiffs[i] <= opt_differences)
                  {
                    targetindices[targetcount] = qgramindices[i];
                    targetampliconids[targetcount] = qgramamps[i];
                    targetcount++;
                  }
          
              if (targetcount > 0)
                {
                  search_do(subseedampliconid, targetcount, targetampliconids, 
                            scores, diffs, alignlengths, bits);
#ifdef VERBOSE
                  searches++;
#endif

                  if (bits == 8)
                    count_comparisons_8 += targetcount;
                  else
                    count_comparisons_16 += targetcount;
            
                  for(unsigned long t=0; t<targetcount; t++)
                    {
                      unsigned diff = diffs[t];
              
                      if (diff <= (unsigned long) opt_differences)
                        {
                          unsigned i = targetindices[t];
                
                          /* find correct position in list */

                          /* move the target (i) to the position (swarmed)
                             of the first unswarmed amplicon in the pool
                             then move the target further into the swarmed
                             but unseeded part of the list, so that the
                             swarmed amplicons are ordered by id */

                          unsigned long targetampliconid = amps[i].ampliconid;
                          unsigned pos = swarmed;

                          while ((pos > seeded) &&
                                 (amps[pos-1].ampliconid > targetampliconid) &&
                                 (amps[pos-1].generation > subseedgeneration))
                            pos--;

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
                          amps[pos].generation = subseedgeneration + 1;
                          if (maxgen < amps[pos].generation)
                            maxgen = amps[pos].generation;
                          amps[pos].radius = subseedradius + diff;
                          if (amps[pos].radius > maxradius)
                            maxradius = amps[pos].radius;

                          unsigned poolampliconid = amps[pos].ampliconid;
                          hits[hitcount++] = poolampliconid;

                          if (opt_internal_structure)
                            {
                              fprint_id_noabundance(internal_structure_file, subseedampliconid);
                              fprintf(internal_structure_file, "\t");
                              fprint_id_noabundance(internal_structure_file, poolampliconid);
                              fprintf(internal_structure_file, "\t%u", diff);
                              fprintf(internal_structure_file, "\t%lu\t%lu", swarmid, subseedgeneration + 1);
                              fprintf(internal_structure_file, "\n");
                            }

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

      if (maxgen > maxgenerations)
        maxgenerations = maxgen;


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
                 score_matrix_63, penalty_gapopen, penalty_gapextend,
                 & nwscore, & nwdiff, & nwalignmentlength, & nwalignment,
                 dir, hearray, 0, 0);
              
              double percentid = 100.0 * (nwalignmentlength - 
                                          nwdiff) / nwalignmentlength;
              
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
          abundance = db_getabundance(seedampliconid);

          fprintf(statsfile, "%lu\t%lu\t", swarmsize, amplicons_copies);
          fprint_id_noabundance(statsfile, seedampliconid);
          fprintf(statsfile, "\t%lu\t%lu\t%lu\t%lu\n", 
                  abundance, singletons, maxgen, maxradius);
        }
      progress_update(seeded);
    }
  progress_done();
  
  if (uclustfile)
    {
      free(dir);
      free(hearray);
    }


  /* output results */
  
  if (amplicons > 0)
    {
      char sep_amplicons;
      char sep_swarms;

      if (opt_mothur)
        {
          /* mothur list file output */
          sep_amplicons = ',';
          sep_swarms = '\t';
          fprintf(outfile, "swarm_%ld\t%lu\t", opt_differences, swarmid);
        }
      else
        {
          /* native swarm output */
          sep_amplicons = SEPCHAR;  /* usually a space */
          sep_swarms = '\n';
        }

      fprint_id(outfile, amps[0].ampliconid);
      long previd = amps[0].swarmid;

      for (unsigned long i=1; i<amplicons; i++)
        {
          long id = amps[i].swarmid;
          if (id == previd)
            fputc(sep_amplicons, outfile);
          else
            fputc(sep_swarms, outfile);
          fprint_id(outfile, amps[i].ampliconid);
          previd = id;
        }

      fputc('\n', outfile);
    }


  /* dump seeds in fasta format with sum of abundances */

  if ((opt_seeds) && (amplicons > 0))
    {
      progress_init("Writing seeds:    ", amplicons);

      unsigned long mass = 0;
      unsigned previd = amps[0].swarmid;
      unsigned seed = amps[0].ampliconid;
      mass += db_getabundance(seed);

      for (unsigned long i=1; i<amplicons; i++)
        {
          unsigned id = amps[i].swarmid;

          if (id != previd)
            {
              fprintf(fp_seeds, ">");
              fprint_id_with_new_abundance(fp_seeds, seed, mass);
              fprintf(fp_seeds, "\n");
              db_fprintseq(fp_seeds, seed, 0);

              mass = 0;
              seed = amps[i].ampliconid;
            }

          mass += db_getabundance(amps[i].ampliconid);
          previd = id;
          progress_update(i);
        }

      fprintf(fp_seeds, ">");
      fprint_id_with_new_abundance(fp_seeds, seed, mass);
      fprintf(fp_seeds, "\n");
      db_fprintseq(fp_seeds, seed, 0);

      progress_done();
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

  db_qgrams_done();

  qgram_diff_done();

  fprintf(logfile, "\n");

  fprintf(logfile, "Number of swarms:  %lu\n", swarmid);

  fprintf(logfile, "Largest swarm:     %lu\n", largestswarm);

  fprintf(logfile, "Max generations:   %lu\n", maxgenerations);

#ifdef VERBOSE
  fprintf(logfile, "\n");

  fprintf(logfile, "Estimates:         %lu\n", estimates);

  fprintf(logfile, "Searches:          %lu\n", searches);

  fprintf(logfile, "\n");

  fprintf(logfile, "Comparisons (8b):  %lu (%.2lf%%)\n",
          count_comparisons_8, (200.0 * count_comparisons_8 / 
                                amplicons / (amplicons+1)));

  fprintf(logfile, "Comparisons (16b): %lu (%.2lf%%)\n",
          count_comparisons_16, (200.0 * count_comparisons_16 /
                                 amplicons / (amplicons+1)));

  fprintf(logfile, "Comparisons (tot): %lu (%.2lf%%)\n",
          count_comparisons_8 + count_comparisons_16,
          (200.0 * (count_comparisons_8 + count_comparisons_16) /
           amplicons / (amplicons+1)));
#endif

}

