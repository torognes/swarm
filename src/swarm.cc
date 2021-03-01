/*
    SWARM

    Copyright (C) 2012-2021 Torbjorn Rognes and Frederic Mahe

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

/* OPTIONS */

std::string progname;  // unused variable?

struct Parameters p;
std::string opt_log;
std::string opt_output_file;
std::string opt_seeds;
std::string opt_statistics_file;
std::string opt_uclust_file;

int64_t opt_append_abundance;
int64_t opt_boundary;
int64_t opt_ceiling;
int64_t opt_gap_extension_penalty;
int64_t opt_gap_opening_penalty;
int64_t opt_match_reward;
int64_t opt_mismatch_penalty;
bool opt_mothur {false};
bool opt_no_otu_breaking {false};
int64_t opt_threads;
bool opt_usearch_abundance {false};

int64_t penalty_factor;
int64_t penalty_gapextend;
int64_t penalty_gapopen;
int64_t penalty_mismatch;

/* Other variables */

int64_t mmx_present {0};
int64_t sse_present {0};
int64_t sse2_present {0};
int64_t sse3_present {0};
int64_t ssse3_present {0};
int64_t sse41_present {0};
int64_t sse42_present {0};
int64_t popcnt_present {0};
int64_t avx_present {0};
int64_t avx2_present {0};

static uint64_t dbsequencecount {0};

uint64_t duplicates_found {0};

FILE * outfile {nullptr};
FILE * statsfile {nullptr};
FILE * uclustfile {nullptr};
FILE * logfile {nullptr};
FILE * internal_structure_file {nullptr};
FILE * fp_seeds {nullptr};
FILE * network_file {nullptr};

const std::array<char, 32> sym_nt =
  {'-', 'A', 'C', 'G', 'T', ' ', ' ', ' ',
   ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
   ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
   ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};

const std::string DASH_FILENAME {"-"};

const std::vector<std::string> header_message
  {swarm_version,
   "\n",
   "Copyright (C) 2012-2021 Torbjorn Rognes and Frederic Mahe\n",
   "https://github.com/torognes/swarm\n",
   "\n",
   "Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M (2014)\n",
   "Swarm: robust and fast clustering method for amplicon-based studies\n",
   "PeerJ 2:e593 https://doi.org/10.7717/peerj.593\n",
   "\n",
   "Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M (2015)\n",
   "Swarm v2: highly-scalable and high-resolution amplicon clustering\n",
   "PeerJ 3:e1420 https://doi.org/10.7717/peerj.1420\n",
   "\n"};

const std::vector<std::string> args_usage_message
  /*0         1         2         3         4         5         6         7          */
  /*01234567890123456789012345678901234567890123456789012345678901234567890123456789 */
  {"Usage: swarm [OPTIONS] [FASTAFILE]\n",
   "\n",
   "General options:\n",
   " -h, --help                          display this help and exit\n",
   " -t, --threads INTEGER               number of threads to use (1)\n",
   " -v, --version                       display version information and exit\n",
   " -x, --disable-sse3                  disable SSE3 and later x86 instructions\n",
   "\n",
   "Clustering options:\n",
   " -d, --differences INTEGER           resolution (1)\n",
   " -n, --no-otu-breaking               never break OTUs (not recommended!)\n",
   "\n",
   "Fastidious options (only when d = 1):\n",
   " -b, --boundary INTEGER              min mass of large OTUs (3)\n",
   " -c, --ceiling INTEGER               max memory in MB for Bloom filter (unlim.)\n",
   " -f, --fastidious                    link nearby low-abundance swarms\n",
   " -y, --bloom-bits INTEGER            bits used per Bloom filter entry (16)\n",
   "\n",
   "Input/output options:\n",
   " -a, --append-abundance INTEGER      value to use when abundance is missing\n",
   " -i, --internal-structure FILENAME   write internal OTU structure to file\n",
   " -j, --network_file FILENAME         dump sequence network to file\n",
   " -l, --log FILENAME                  log to file, not to stderr\n",
   " -o, --output-file FILENAME          output result to file (stdout)\n",
   " -r, --mothur                        output using mothur-like format\n",
   " -s, --statistics-file FILENAME      dump OTU statistics to file\n",
   " -u, --uclust-file FILENAME          output using UCLUST-like format to file\n",
   " -w, --seeds FILENAME                write OTU representatives to FASTA file\n",
   " -z, --usearch-abundance             abundance annotation in usearch style\n",
   "\n",
   "Pairwise alignment advanced options (only when d > 1):\n",
   " -m, --match-reward INTEGER          reward for nucleotide match (5)\n",
   " -p, --mismatch-penalty INTEGER      penalty for nucleotide mismatch (4)\n",
   " -g, --gap-opening-penalty INTEGER   gap open penalty (12)\n",
   " -e, --gap-extension-penalty INTEGER gap extension penalty (4)\n",
#ifndef __WIN32
   "\n",
   "See 'man swarm' for more details.\n",
#endif
   "\n"
  };


#ifdef __x86_64__

void cpuid(unsigned int f1,
           unsigned int f2,
           unsigned int & a,
           unsigned int & b,
           unsigned int & c,
           unsigned int & d);

void cpu_features_detect();

void cpu_features_show();

void cpuid(unsigned int f1,
           unsigned int f2,
           unsigned int & a,
           unsigned int & b,
           unsigned int & c,
           unsigned int & d)
{
  __asm__ __volatile__ ("cpuid"
                        : "=a" (a), "=b" (b), "=c" (c), "=d" (d)
                        : "a" (f1), "c" (f2));
}

void cpu_features_detect()
{
  constexpr unsigned int post_pentium {7};  // new cpus: a & 0xff > 6
  constexpr unsigned int bit_mmx {23};
  constexpr unsigned int bit_sse {25};
  constexpr unsigned int bit_sse2 {26};
  constexpr unsigned int bit_sse3 {0};
  constexpr unsigned int bit_ssse3 {9};
  constexpr unsigned int bit_sse41 {19};
  constexpr unsigned int bit_sse42 {20};
  constexpr unsigned int bit_popcnt {23};
  constexpr unsigned int bit_avx {28};

  unsigned int a {0};
  unsigned int b {0};
  unsigned int c {0};
  unsigned int d {0};

  cpuid(0, 0, a, b, c, d);
  unsigned int maxlevel = a & UINT8_MAX;

  if (maxlevel >= 1)
  {
    cpuid(1, 0, a, b, c, d);
    mmx_present    = (d >> bit_mmx) & 1;
    sse_present    = (d >> bit_sse) & 1;
    sse2_present   = (d >> bit_sse2) & 1;
    sse3_present   = (c >> bit_sse3) & 1;
    ssse3_present  = (c >> bit_ssse3) & 1;
    sse41_present  = (c >> bit_sse41) & 1;
    sse42_present  = (c >> bit_sse42) & 1;
    popcnt_present = (c >> bit_popcnt) & 1;
    avx_present    = (c >> bit_avx) & 1;

    if (maxlevel >= post_pentium)
    {
      cpuid(post_pentium, 0, a, b, c, d);
      constexpr unsigned int bit_avx2 {5};
      avx2_present   = (b >> bit_avx2) & 1;
    }
  }
}

void cpu_features_show()
{
  fprintf(logfile, "CPU features:     ");
  if (mmx_present != 0){
    fprintf(logfile, " mmx");
  }
  if (sse_present != 0) {
    fprintf(logfile, " sse");
  }
  if (sse2_present != 0) {
    fprintf(logfile, " sse2");
  }
  if (sse3_present != 0) {
    fprintf(logfile, " sse3");
  }
  if (ssse3_present != 0) {
    fprintf(logfile, " ssse3"); // Supplemental SSSE3, introduced in 2006
  }
  if (sse41_present != 0) {
    fprintf(logfile, " sse4.1");
  }
  if (sse42_present != 0) {
    fprintf(logfile, " sse4.2");
  }
  if (popcnt_present != 0) {
    fprintf(logfile, " popcnt");
  }
  if (avx_present != 0) {
    fprintf(logfile, " avx");
  }
  if (avx2_present != 0) {
    fprintf(logfile, " avx2");
  }
  fprintf(logfile, "\n");
}

#endif

auto args_long(char * str, const char * option) -> int64_t;
void args_show();
void show(const std::vector<std::string> & message);
void args_init(int argc, char **argv);
void open_files();
void close_files();

auto args_long(char * str, const char * option) -> int64_t
{
  constexpr int base_value {10};
  char * endptr {nullptr};
  int64_t temp = strtol(str, & endptr, base_value);
  if (*endptr != 0)
    {
      fprintf(stderr, "\nInvalid numeric argument for option %s\n", option);
      exit(1);
    }
  return temp;
}

void args_show()
{
#ifdef __x86_64__
  cpu_features_show();
#endif

  fprintf(logfile, "Database file:     %s\n", p.input_filename.c_str());
  fprintf(logfile, "Output file:       %s\n", opt_output_file.c_str());
  if (! opt_statistics_file.empty()) {
    fprintf(logfile, "Statistics file:   %s\n", opt_statistics_file.c_str());
  }
  if (! opt_uclust_file.empty()) {
    fprintf(logfile, "Uclust file:       %s\n", opt_uclust_file.c_str());
  }
  if (! p.opt_internal_structure.empty()) {
    fprintf(logfile, "Int. struct. file  %s\n", p.opt_internal_structure.c_str());
  }
  if (! p.opt_network_file.empty()) {
    fprintf(logfile, "Network file       %s\n", p.opt_network_file.c_str());
  }
  fprintf(logfile, "Resolution (d):    %" PRId64 "\n", p.opt_differences);
  fprintf(logfile, "Threads:           %" PRId64 "\n", opt_threads);

  if (p.opt_differences > 1)
    {
      fprintf(logfile,
              "Scores:            match: %" PRId64 ", mismatch: %" PRId64 "\n",
              opt_match_reward, opt_mismatch_penalty);
      fprintf(logfile,
              "Gap penalties:     opening: %" PRId64 ", extension: %" PRId64 "\n",
              opt_gap_opening_penalty, opt_gap_extension_penalty);
      fprintf(logfile,
              "Converted costs:   mismatch: %" PRId64 ", gap opening: %" PRId64 ", "
              "gap extension: %" PRId64 "\n",
              penalty_mismatch, penalty_gapopen, penalty_gapextend);
    }
  fprintf(logfile, "Break OTUs:        %s\n",
          opt_no_otu_breaking ? "No" : "Yes");
  if (p.opt_fastidious) {
    fprintf(logfile, "Fastidious:        Yes, with boundary %" PRId64 "\n",
            opt_boundary);
  }
  else {
    fprintf(logfile, "Fastidious:        No\n");
  }
}

void show(const std::vector<std::string> & message)
{
  for (const auto & m : message) {
    fprintf(logfile, "%s", m.c_str());
  }
}

void args_init(int argc, char **argv)
{
  /* Set defaults */
  constexpr unsigned int min_ceiling {8};
  constexpr unsigned int max_ceiling {1 << 30};  // 1,073,741,824 (MiB of RAM)
  constexpr unsigned int append_abundance_default {0};
  constexpr unsigned int boundary_default {3};
  constexpr unsigned int ceiling_default {0};
  constexpr unsigned int gap_extension_penalty_default {4};
  constexpr unsigned int gap_opening_penalty_default {12};
  constexpr unsigned int match_reward_default {5};
  constexpr unsigned int mismatch_penalty_default {4};
  constexpr unsigned int threads_default {1};
  constexpr unsigned int max_threads {256};

  progname = argv[0];

  opt_append_abundance = append_abundance_default;
  opt_boundary = boundary_default;
  opt_ceiling = ceiling_default;
  opt_gap_extension_penalty = gap_extension_penalty_default;
  opt_gap_opening_penalty = gap_opening_penalty_default;
  opt_match_reward = match_reward_default;
  opt_mismatch_penalty = mismatch_penalty_default;
  opt_output_file = DASH_FILENAME;
  opt_threads = threads_default;
  opterr = 1;  // unused variable?

  char short_options[] = "a:b:c:d:e:fg:hi:j:l:m:no:p:rs:t:u:vw:xy:z";

  /* unused short option letters: kq */

  static struct option long_options[] =
  {
    {"append-abundance",      required_argument, nullptr, 'a' },
    {"boundary",              required_argument, nullptr, 'b' },
    {"ceiling",               required_argument, nullptr, 'c' },
    {"differences",           required_argument, nullptr, 'd' },
    {"gap-extension-penalty", required_argument, nullptr, 'e' },
    {"fastidious",            no_argument,       nullptr, 'f' },
    {"gap-opening-penalty",   required_argument, nullptr, 'g' },
    {"help",                  no_argument,       nullptr, 'h' },
    {"internal-structure",    required_argument, nullptr, 'i' },
    {"log",                   required_argument, nullptr, 'l' },
    {"network-file",          required_argument, nullptr, 'j' },
    {"match-reward",          required_argument, nullptr, 'm' },
    {"no-otu-breaking",       no_argument,       nullptr, 'n' },
    {"output-file",           required_argument, nullptr, 'o' },
    {"mismatch-penalty",      required_argument, nullptr, 'p' },
    {"mothur",                no_argument,       nullptr, 'r' },
    {"statistics-file",       required_argument, nullptr, 's' },
    {"threads",               required_argument, nullptr, 't' },
    {"uclust-file",           required_argument, nullptr, 'u' },
    {"version",               no_argument,       nullptr, 'v' },
    {"seeds",                 required_argument, nullptr, 'w' },
    {"disable-sse3",          no_argument,       nullptr, 'x' },
    {"bloom-bits",            required_argument, nullptr, 'y' },
    {"usearch-abundance",     no_argument,       nullptr, 'z' },
    {nullptr,                 0,                 nullptr, 0 }
  };

  constexpr int n_options {26};
  std::array<int, n_options> used_options {{}};
  used_options.fill(0);  // set int values to zero by default

  int option_index {0};
  int c {0};

  while ((c = getopt_long(argc, argv, short_options, long_options, &option_index)) != -1)
  {

    /* check if any option is specified more than once */

    if ((c >= 'a') && (c <= 'z'))
      {
        auto optindex = static_cast<unsigned int>(c - 'a');  // c - 'a' cannot be negative
        if (used_options[optindex] == 1)
          {
            int longoptindex {0};
            while (long_options[longoptindex].name != nullptr)
              {
                if (long_options[longoptindex].val == c) {
                  break;
                }
                longoptindex++;
              }

            fprintf(stderr,
                    "Error: Option -%c or --%s specified more than once.\n",
                    c,
                    long_options[longoptindex].name);
            exit(1);
          }
        used_options[optindex] = 1;
      }

    switch(c)
      {
      case 'a':
        /* append-abundance */
        opt_append_abundance = args_long(optarg, "-a or --append-abundance");
        break;

      case 'b':
        /* boundary */
        opt_boundary = args_long(optarg, "-b or --boundary");
        break;

      case 'c':
        /* ceiling */
        opt_ceiling = args_long(optarg, "-c or --ceiling");
        break;

      case 'd':
        /* differences (resolution) */
        p.opt_differences = args_long(optarg, "-d or --differences");
        break;

      case 'e':
        /* gap extension penalty */
        opt_gap_extension_penalty = args_long(optarg, "-e or --gap-extension-penalty");
        break;

      case 'f':
        /* fastidious */
        p.opt_fastidious = true;
        break;

      case 'g':
        /* gap-opening-penalty */
        opt_gap_opening_penalty = args_long(optarg, "-g or --gap-opening-penalty");
        break;

      case 'h':
        /* help */
        p.opt_help = true;
        break;

      case 'i':
        /* internal-structure */
        p.opt_internal_structure = optarg;
        break;

      case 'j':
        /* network-file */
        p.opt_network_file = optarg;
        break;

      case 'l':
        /* log */
        opt_log = optarg;
        break;

      case 'm':
        /* match-reward */
        opt_match_reward = args_long(optarg, "-m or --match-reward");
        break;

      case 'n':
        /* no-otu-breaking */
        opt_no_otu_breaking = true;
        break;

      case 'o':
        /* output-file */
        opt_output_file = optarg;
        break;

      case 'p':
        /* mismatch-penalty */
        opt_mismatch_penalty = args_long(optarg, "-p or --mismatch-penalty");
        break;

      case 'r':
        /* mothur */
        opt_mothur = true;
        break;

      case 's':
        /* statistics-file */
        opt_statistics_file = optarg;
        break;

      case 't':
        /* threads */
        opt_threads = args_long(optarg, "-t or --threads");
        break;

      case 'u':
        /* uclust-file */
        opt_uclust_file = optarg;
        break;

      case 'v':
        /* version */
        p.opt_version = true;
        break;

      case 'w':
        /* seeds */
        opt_seeds = optarg;
        break;

      case 'x':
        /* disable-sse3 */
        p.opt_disable_sse3 = true;
        break;

      case 'y':
        /* bloom-bits */
        p.opt_bloom_bits = args_long(optarg, "-y or --bloom-bits");
        break;

      case 'z':
        /* usearch-abundance */
        opt_usearch_abundance = true;
        break;

      default:
        show(header_message);
        show(args_usage_message);
        exit(1);
    }
  }

  if (optind < argc) {
    p.input_filename = argv[optind];
  }

  if ((opt_threads < 1) || (opt_threads > max_threads))
    {
      fprintf(stderr, "\nError: Illegal number of threads specified with "
              "-t or --threads, must be in the range 1 to %u.\n", max_threads);
      exit(1);
    }

  if ((p.opt_differences < 0) || (p.opt_differences > UINT8_MAX)) {
    fatal("Illegal number of differences specified with -d or --differences, "
          "must be in the range 0 to 255.");
  }

  if (p.opt_fastidious && (p.opt_differences != 1)) {
    fatal("Fastidious mode (specified with -f or --fastidious) only works "
          "when the resolution (specified with -d or --differences) is 1.");
  }

  if (p.opt_disable_sse3 && (p.opt_differences < 2)) {
    fatal("Option --disable-sse3 or -x has no effect when d < 2. "
          "(SSE3 instructions are only used when d > 1).");
  }

  if (! p.opt_fastidious)
    {
      constexpr unsigned int boundary_index {1};
      constexpr unsigned int ceiling_index {2};
      constexpr unsigned int bloom_bits_index {24};

      if (used_options[boundary_index] != 0) {
        fatal("Option -b or --boundary specified without -f or --fastidious.\n");
      }
      if (used_options[ceiling_index] != 0) {
        fatal("Option -c or --ceiling specified without -f or --fastidious.\n");
      }
      if (used_options[bloom_bits_index] != 0) {
        fatal("Option -y or --bloom-bits specified without -f or --fastidious.\n");
      }
    }

  if (p.opt_differences < 2)
    {
      constexpr unsigned int match_reward_index {12};
      constexpr unsigned int mismatch_penalty_index {15};
      constexpr unsigned int gap_opening_penalty_index {6};
      constexpr unsigned int gap_extension_penalty_index {4};

      if (used_options[match_reward_index] != 0) {
        fatal("Option -m or --match-reward specified when d < 2.");
      }
      if (used_options[mismatch_penalty_index] != 0) {
        fatal("Option -p or --mismatch-penalty specified when d < 2.");
      }
      if (used_options[gap_opening_penalty_index] != 0) {
        fatal("Option -g or --gap-opening-penalty specified when d < 2.");
      }
      if (used_options[gap_extension_penalty_index] != 0) {
        fatal("Option -e or --gap-extension-penalty specified when d < 2.");
      }
    }

  if (opt_gap_opening_penalty < 0) {
    fatal("Illegal gap opening penalty specified with -g or "
          "--gap-opening-penalty, must not be negative.");
  }

  if (opt_gap_extension_penalty < 0) {
    fatal("Illegal gap extension penalty specified with -e or "
          "--gap-extension-penalty, must not be negative.");
  }

  if ((opt_gap_opening_penalty + opt_gap_extension_penalty) < 1) {
    fatal("Illegal gap penalties specified, the sum of the gap open and "
          "the gap extension penalty must be at least 1.");
  }

  if (opt_match_reward < 1) {
    fatal("Illegal match reward specified with -m or --match-reward, "
          "must be at least 1.");
  }

  if (opt_mismatch_penalty < 1) {
    fatal("Illegal mismatch penalty specified with -p or --mismatch-penalty, "
          "must be at least 1.");
  }

  if (opt_boundary < 2) {
    fatal("Illegal boundary specified with -b or --boundary, "
          "must be at least 2.");
  }

  if ((used_options[2] != 0) && ((opt_ceiling < min_ceiling) ||
                                 (opt_ceiling > max_ceiling))) {
    fatal("Illegal memory ceiling specified with -c or --ceiling, "
          "must be in the range 8 to 1,073,741,824 MB.");
  }

  constexpr unsigned int min_bits_per_entry {2};
  constexpr unsigned int max_bits_per_entry {64};
  if ((p.opt_bloom_bits < min_bits_per_entry) ||
      (p.opt_bloom_bits > max_bits_per_entry)) {
    fatal("Illegal number of Bloom filter bits specified with -y or "
          "--bloom-bits, must be in the range 2 to 64.");
  }

  if ((used_options[0] != 0) && (opt_append_abundance < 1)) {
    fatal("Illegal abundance value specified with -a or --append-abundance, "
          "must be at least 1.");
  }

  if ((! p.opt_network_file.empty()) && (p.opt_differences != 1)) {
    fatal("A network file can only written when d=1.");
  }
}

void open_files()
{
  /* open files */

  if (! opt_log.empty())
    {
      logfile = fopen_output(opt_log.c_str());
      if (logfile == nullptr) {
        fatal("Unable to open log file for writing.");
      }
    }

  outfile = fopen_output(opt_output_file.c_str());
  if (outfile == nullptr) {
    fatal("Unable to open output file for writing.");
  }

  if (! opt_seeds.empty())
    {
      fp_seeds = fopen_output(opt_seeds.c_str());
      if (fp_seeds == nullptr) {
        fatal("Unable to open seeds file for writing.");
      }
    }

  if (! opt_statistics_file.empty())
    {
      statsfile = fopen_output(opt_statistics_file.c_str());
      if (statsfile == nullptr) {
        fatal("Unable to open statistics file for writing.");
      }
    }

  if (! opt_uclust_file.empty())
    {
      uclustfile = fopen_output(opt_uclust_file.c_str());
      if (uclustfile == nullptr) {
        fatal("Unable to open uclust file for writing.");
      }
    }

  if (! p.opt_internal_structure.empty())
    {
      internal_structure_file = fopen_output(p.opt_internal_structure.c_str());
      if (internal_structure_file == nullptr) {
        fatal("Unable to open internal structure file for writing.");
      }
    }

  if (! p.opt_network_file.empty())
    {
      network_file = fopen_output(p.opt_network_file.c_str());
      if (network_file == nullptr) {
        fatal("Unable to open network file for writing.");
      }
    }
}

void close_files() {
  const std::vector<FILE *> file_handles
    {network_file, internal_structure_file,
     uclustfile, statsfile, fp_seeds, outfile,
     logfile};
  for (auto * const file_handle : file_handles) {
    if (file_handle != nullptr) {
      fclose(file_handle);
    }
  }
}

auto main(int argc, char** argv) -> int
{
  logfile = stderr;

  arch_srandom(1);

  args_init(argc, argv);

#ifdef __x86_64__
  cpu_features_detect();

  if (sse2_present == 0) {
    fatal("This program requires a processor with SSE2 instructions.\n");
  }

  if (p.opt_disable_sse3)
    {
      sse3_present = 0;
      ssse3_present = 0;
      sse41_present = 0;
      sse42_present = 0;
      popcnt_present = 0;
      avx_present = 0;
      avx2_present = 0;
    }
#endif

  open_files();

  if (p.opt_version || p.opt_help)
    {
      show(header_message);
      if (p.opt_help) {
        show(args_usage_message);
      }
      close_files();
      exit(0);
    }

  penalty_mismatch = 2 * opt_match_reward + 2 * opt_mismatch_penalty;
  penalty_gapopen = 2 * opt_gap_opening_penalty;
  penalty_gapextend = opt_match_reward + 2 * opt_gap_extension_penalty;

  penalty_factor = gcd(gcd(penalty_mismatch, penalty_gapopen), penalty_gapextend);

  penalty_mismatch /= penalty_factor;
  penalty_gapopen /= penalty_factor;
  penalty_gapextend /= penalty_factor;

  int64_t diff_saturation_16 = (std::min((UINT16_MAX / penalty_mismatch),
                                         (UINT16_MAX - penalty_gapopen)
                                         / penalty_gapextend));

  if (p.opt_differences > diff_saturation_16)
    fatal("Resolution (d) too high for the given scoring system");

  show(header_message);

  args_show();

  fprintf(logfile, "\n");

  db_read(p.input_filename.c_str(), p);

  fprintf(logfile, "Database info:     %" PRIu64 " nt", db_getnucleotidecount());
  fprintf(logfile, " in %u sequences,", db_getsequencecount());
  fprintf(logfile, " longest %u nt\n", db_getlongestsequence());

  dbsequencecount = db_getsequencecount();

  score_matrix_init();

  switch (p.opt_differences)
    {
    case 0:
      dereplicate(p);
      break;

    case 1:
      algo_d1_run(p);
      break;

    default:
      algo_run(p);
      break;
    }

  score_matrix_free();

  db_free();

  close_files();
}
