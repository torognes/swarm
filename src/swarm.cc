/*
    SWARM

    Copyright (C) 2012-2024 Torbjorn Rognes and Frederic Mahe

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
#include "algo.h"
#include "algod1.h"
#include "db.h"
#include "derep.h"
#include "utils/alignment_parameters.h"
#include "utils/fatal.h"
#include "utils/gcd.h"
#include "utils/input_output.h"
#include "utils/opt_boundary.h"
#include "utils/opt_no_cluster_breaking.h"
#include "utils/seqinfo.h"
#include "x86_cpu_features.h"
#include "zobrist.h"
#include <algorithm>  // std::min()
#include <array>
#include <cassert>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t
#include <cstdio>  // FILE, fclose, stderr  // refactoring: replace with <fstream>
#include <cstdlib>  // std::exit
#include <getopt.h>  // getopt_long, optarg, optind, opterr, struct
                     // option (no_argument, required_argument)
#include <limits>
#include <string>
// #include <unistd.h>  // getopt_long... FAIL
#include <vector>

#ifndef PRId64
#ifdef _WIN32
#define PRId64 "I64d"
#else
constexpr char PRId64[] = "ld";
#endif
#endif


constexpr auto int_max = std::numeric_limits<int>::max();
constexpr auto int16_max = std::numeric_limits<int16_t>::max();
static_assert(int_max > int16_max, "Your compiler uses very short integers.");
const std::string swarm_version {"3.1.4"};


/* OPTIONS */

struct Parameters parameters;
std::string opt_log;

int64_t opt_boundary;
bool opt_no_cluster_breaking {false};
int64_t opt_threads;

int64_t penalty_mismatch;
int64_t penalty_gapextend;
int64_t penalty_gapopen;


/* fine names and command line options */

std::FILE * outfile {nullptr};
std::FILE * statsfile {nullptr};
std::FILE * uclustfile {nullptr};
std::FILE * logfile {stderr};  // cstdio stderr macro is expanded to type std::FILE*
std::FILE * internal_structure_file {nullptr};
std::FILE * fp_seeds {nullptr};
std::FILE * network_file {nullptr};

constexpr int n_options {26};

std::string short_options = "a:b:c:d:e:fg:hi:j:l:m:no:p:rs:t:u:vw:xy:z"; /* unused: kq */

// refactoring: add option -q (no-cluster-breaking)
const std::array<struct option, 25> long_options = {
  { // struct option { name, has_arg, flag, val }
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
  }
};


const std::vector<std::string> header_message
  {"Swarm ", swarm_version,
   "\n",
   "Copyright (C) 2012-2024 Torbjorn Rognes and Frederic Mahe\n",
   "https://github.com/torognes/swarm\n",
   "\n",
   "Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M (2014)\n",
   "Swarm: robust and fast clustering method for amplicon-based studies\n",
   "PeerJ 2:e593 https://doi.org/10.7717/peerj.593\n",
   "\n",
   "Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M (2015)\n",
   "Swarm v2: highly-scalable and high-resolution amplicon clustering\n",
   "PeerJ 3:e1420 https://doi.org/10.7717/peerj.1420\n",
   "\n",
   "Mahe F, Czech L, Stamatakis A, Quince C, de Vargas C, Dunthorn M, Rognes T (2022)\n",
   "Swarm v3: towards tera-scale amplicon clustering\n",
   "Bioinformatics 38:1, 267-269 https://doi.org/10.1093/bioinformatics/btab493\n",
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
   "\n",
   "Clustering options:\n",
   " -d, --differences INTEGER           resolution (1)\n",
   " -n, --no-otu-breaking               never break clusters (not recommended!)\n",
   "\n",
   "Fastidious options (only when d = 1):\n",
   " -b, --boundary INTEGER              min mass of large clusters (3)\n",
   " -c, --ceiling INTEGER               max memory in MB for Bloom filter (unlim.)\n",
   " -f, --fastidious                    link nearby low-abundance swarms\n",
   " -y, --bloom-bits INTEGER            bits used per Bloom filter entry (16)\n",
   "\n",
   "Input/output options:\n",
   " -a, --append-abundance INTEGER      value to use when abundance is missing\n",
   " -i, --internal-structure FILENAME   write internal cluster structure to file\n",
   " -j, --network_file FILENAME         dump sequence network to file\n",
   " -l, --log FILENAME                  log to file, not to stderr\n",
   " -o, --output-file FILENAME          output result to file (stdout)\n",
   " -r, --mothur                        output using mothur-like format\n",
   " -s, --statistics-file FILENAME      dump cluster statistics to file\n",
   " -u, --uclust-file FILENAME          output using UCLUST-like format to file\n",
   " -w, --seeds FILENAME                write cluster representatives to FASTA file\n",
   " -z, --usearch-abundance             abundance annotation in usearch style\n",
   "\n",
   "Pairwise alignment advanced options (only when d > 1):\n",
   " -m, --match-reward INTEGER          reward for nucleotide match (5)\n",
   " -p, --mismatch-penalty INTEGER      penalty for nucleotide mismatch (4)\n",
   " -g, --gap-opening-penalty INTEGER   gap open penalty (12)\n",
   " -e, --gap-extension-penalty INTEGER gap extension penalty (4)\n",
   " -x, --disable-sse3                  disable SSE3 and later x86 instructions\n",
#ifndef __WIN32
   "\n",
   "See 'man swarm' for more details.\n",
#endif
   "\n"
  };


auto args_long(char * str, const char * option) -> int64_t
{
  static constexpr int base_value {10};
  char * endptr {nullptr};
  const int64_t temp = strtol(str, & endptr, base_value);
  if (*endptr != 0)
    {
      fatal(error_prefix, "Invalid numeric argument for option ", option, ".\n\n",
            "Frequent causes are:\n",
            " - a missing space between an argument and the next option,\n",
            " - a long option name not starting with a double dash\n",
            "   (swarm accepts '--help' or '-h', but not '-help')\n\n",
            "Please see 'swarm --help' for more details.");
    }
  return temp;
}


auto args_show() -> void  // refactoring: pass a ref to parameters as argument
{
#ifdef __x86_64__
  cpu_features_detect(parameters);
  cpu_features_test(parameters);
  cpu_features_show(parameters);
#endif

  std::fprintf(logfile, "Database file:     %s\n", parameters.input_filename.c_str());
  std::fprintf(logfile, "Output file:       %s\n", parameters.opt_output_file.c_str());
  if (not parameters.opt_statistics_file.empty()) {
    std::fprintf(logfile, "Statistics file:   %s\n", parameters.opt_statistics_file.c_str());
  }
  if (not parameters.opt_uclust_file.empty()) {
    std::fprintf(logfile, "Uclust file:       %s\n", parameters.opt_uclust_file.c_str());
  }
  if (not parameters.opt_internal_structure.empty()) {
    std::fprintf(logfile, "Int. struct. file  %s\n", parameters.opt_internal_structure.c_str());
  }
  if (not parameters.opt_network_file.empty()) {
    std::fprintf(logfile, "Network file       %s\n", parameters.opt_network_file.c_str());
  }
  std::fprintf(logfile, "Resolution (d):    %" PRId64 "\n", parameters.opt_differences);
  std::fprintf(logfile, "Threads:           %" PRId64 "\n", opt_threads);

  if (parameters.opt_differences > 1)
    {
      std::fprintf(logfile,
              "Scores:            match: %" PRId64 ", mismatch: %" PRId64 "\n",
              parameters.opt_match_reward, parameters.opt_mismatch_penalty);
      std::fprintf(logfile,
              "Gap penalties:     opening: %" PRId64 ", extension: %" PRId64 "\n",
              parameters.opt_gap_opening_penalty, parameters.opt_gap_extension_penalty);
      std::fprintf(logfile,
              "Converted costs:   mismatch: %" PRId64 ", gap opening: %" PRId64 ", "
              "gap extension: %" PRId64 "\n",
              parameters.penalty_mismatch, parameters.penalty_gapopen, parameters.penalty_gapextend);
    }
  std::fprintf(logfile, "Break clusters:        %s\n",
          parameters.opt_no_cluster_breaking ? "No" : "Yes");
  if (parameters.opt_fastidious) {
    std::fprintf(logfile, "Fastidious:        Yes, with boundary %" PRId64 "\n",
            parameters.opt_boundary);
  }
  else {
    std::fprintf(logfile, "Fastidious:        No\n");
  }
  std::fprintf(logfile, "\n");
}


auto show(const std::vector<std::string> & message) -> void
{
  for (const auto & message_element : message) {
    std::fprintf(logfile, "%s", message_element.c_str());
  }
}


auto args_init(int argc, char **argv) -> std::array<bool, n_options>
{
  /* Set defaults */
  static constexpr unsigned int threads_default {1};
  std::array<bool, n_options> used_options {{}};  // value initialization sets values to 'false'

  opt_boundary = parameters.opt_boundary;
  opt_threads = threads_default;
  opterr = 1;  // unused variable? get_opt option?

  int option_character {0};

  while (true)
  {
    int option_index {0};
    option_character = getopt_long(argc, argv, short_options.c_str(), long_options.data(), &option_index);

    if (option_character == -1) {
      break;
    }

    /* check if any option is specified more than once */

    if ((option_character >= 'a') and (option_character <= 'z'))
      {
        assert(option_character - 'a' >= 0);
        auto optindex = static_cast<unsigned int>(option_character - 'a');
        if (used_options[optindex])
          {
            auto longoptindex = 0UL;
            for (const auto& long_option: long_options) {
              assert(long_option.name != nullptr);
              if (long_option.name == nullptr) {
                break; // refactoring: unreachable?
              }

              if (long_option.val == option_character) {
                break;
              }
              ++longoptindex;
            }

            fatal(error_prefix, "Option -", static_cast<char>(option_character),
                  " or --", long_options[longoptindex].name,
                  " specified more than once.");
          }
        used_options[optindex] = true;
      }

    switch(option_character)
      {
      case 'a':  // refactoring: replace with enum class? static_cast<Short_option>(option_character), but what about unexpected characters such as 'k' or 'q'?
        /* append-abundance */
        parameters.opt_append_abundance = args_long(optarg, "-a or --append-abundance");
        break;

      case 'b':
        /* boundary */
        parameters.opt_boundary = args_long(optarg, "-b or --boundary");
        opt_boundary = parameters.opt_boundary;
        break;

      case 'c':
        /* ceiling */
        parameters.opt_ceiling = args_long(optarg, "-c or --ceiling");
        break;

      case 'd':
        /* differences (resolution) */
        parameters.opt_differences = args_long(optarg, "-d or --differences");
        break;

      case 'e':
        /* gap extension penalty */
        parameters.opt_gap_extension_penalty = args_long(optarg, "-e or --gap-extension-penalty");
        break;

      case 'f':
        /* fastidious */
        parameters.opt_fastidious = true;
        break;

      case 'g':
        /* gap-opening-penalty */
        parameters.opt_gap_opening_penalty = args_long(optarg, "-g or --gap-opening-penalty");
        break;

      case 'h':
        /* help */
        parameters.opt_help = true;
        break;

      case 'i':
        /* internal-structure */
        parameters.opt_internal_structure = optarg;
        break;

      case 'j':
        /* network-file */
        parameters.opt_network_file = optarg;
        break;

      case 'l':
        /* log */
        opt_log = optarg;
        break;

      case 'm':
        /* match-reward */
        parameters.opt_match_reward = args_long(optarg, "-m or --match-reward");
        break;

      case 'n':
        /* no-cluster-breaking */
        parameters.opt_no_cluster_breaking = true;
        opt_no_cluster_breaking = parameters.opt_no_cluster_breaking;
        break;

      case 'o':
        /* output-file */
        parameters.opt_output_file = optarg;
        break;

      case 'p':
        /* mismatch-penalty */
        parameters.opt_mismatch_penalty = args_long(optarg, "-p or --mismatch-penalty");
        break;

      case 'r':
        /* mothur */
        parameters.opt_mothur = true;
        break;

      case 's':
        /* statistics-file */
        parameters.opt_statistics_file = optarg;
        break;

      case 't':
        /* threads */
        opt_threads = args_long(optarg, "-t or --threads");
        break;

      case 'u':
        /* uclust-file */
        parameters.opt_uclust_file = optarg;
        break;

      case 'v':
        /* version */
        parameters.opt_version = true;
        break;

      case 'w':
        /* seeds */
        parameters.opt_seeds = optarg;
        break;

      case 'x':
        /* disable-sse3 */
        parameters.opt_disable_sse3 = true;
        break;

      case 'y':
        /* bloom-bits */
        parameters.opt_bloom_bits = args_long(optarg, "-y or --bloom-bits");
        break;

      case 'z':
        /* usearch-abundance */
        parameters.opt_usearch_abundance = true;
        break;

      default:
        show(header_message);
        show(args_usage_message);
        fatal();
    }
  }

  if (optind < argc) {  // external variable defined in unistd.h for
                        // use with the getopt function
    parameters.input_filename = argv[optind];
  }

  // scoring system
  parameters.penalty_mismatch = 2 * parameters.opt_match_reward + 2 * parameters.opt_mismatch_penalty;
  penalty_mismatch = parameters.penalty_mismatch;
  parameters.penalty_gapopen = 2 * parameters.opt_gap_opening_penalty;
  penalty_gapopen = parameters.penalty_gapopen;
  parameters.penalty_gapextend = parameters.opt_match_reward + 2 * parameters.opt_gap_extension_penalty;
  penalty_gapextend = parameters.penalty_gapextend;

  const int64_t penalty_factor {gcd(gcd(parameters.penalty_mismatch, parameters.penalty_gapopen), parameters.penalty_gapextend)};

  // clang: risk of DivideZero, but that would require gcd(0, 0) which is not possible
  parameters.penalty_mismatch /= penalty_factor;
  penalty_mismatch = parameters.penalty_mismatch;
  parameters.penalty_gapopen /= penalty_factor;
  penalty_gapopen = parameters.penalty_gapopen;
  parameters.penalty_gapextend /= penalty_factor;
  penalty_gapextend = parameters.penalty_gapextend;

  return used_options;
}


auto args_check(const std::array<bool, n_options> & used_options) -> void {
  static constexpr auto uint8_max = std::numeric_limits<uint8_t>::max();
  static constexpr auto uint16_max = std::numeric_limits<uint16_t>::max();
  static constexpr unsigned int min_bits_per_entry {2};
  static constexpr unsigned int max_bits_per_entry {64};
  static constexpr unsigned int min_ceiling {8};
  static constexpr unsigned int max_ceiling {1U << 30U};  // 1,073,741,824 (MiB of RAM)
  static constexpr unsigned int max_threads {512};
  // meaning of the used_options values
  static constexpr unsigned int append_abundance_index {0};
  static constexpr unsigned int boundary_index {1};
  static constexpr unsigned int ceiling_index {2};
  static constexpr unsigned int gap_extension_penalty_index {4};
  static constexpr unsigned int gap_opening_penalty_index {6};
  static constexpr unsigned int match_reward_index {12};
  static constexpr unsigned int mismatch_penalty_index {15};
  static constexpr unsigned int bloom_bits_index {24};

  if ((opt_threads < 1) or (opt_threads > max_threads))
    {
      fatal(error_prefix, "Illegal number of threads specified with "
            "-t or --threads, must be in the range 1 to ", max_threads, ".");
    }

  if ((parameters.opt_differences < 0) or (parameters.opt_differences > uint8_max)) {
    fatal(error_prefix, "Illegal number of differences specified with -d or --differences, "
          "must be in the range 0 to ", uint8_max, ".");
  }

  if (parameters.opt_fastidious and (parameters.opt_differences != 1)) {
    fatal(error_prefix, "Fastidious mode (specified with -f or --fastidious) only works "
          "when the resolution (specified with -d or --differences) is 1.");
  }

  if (parameters.opt_disable_sse3 and (parameters.opt_differences < 2)) {
    fatal(error_prefix, "Option --disable-sse3 or -x has no effect when d < 2 "
          "(SSE3 instructions are only used when d > 1).");
  }

  if (not parameters.opt_fastidious)
    {
      if (used_options[boundary_index]) {
        fatal(error_prefix, "Option -b or --boundary specified without -f or --fastidious.");
      }
      if (used_options[ceiling_index]) {
        fatal(error_prefix, "Option -c or --ceiling specified without -f or --fastidious.");
      }
      if (used_options[bloom_bits_index]) {
        fatal(error_prefix, "Option -y or --bloom-bits specified without -f or --fastidious.");
      }
    }

  if (parameters.opt_differences < 2)
    {
      if (used_options[match_reward_index]) {
        fatal(error_prefix, "Option -m or --match-reward specified when d < 2.");
      }
      if (used_options[mismatch_penalty_index]) {
        fatal(error_prefix, "Option -p or --mismatch-penalty specified when d < 2.");
      }
      if (used_options[gap_opening_penalty_index]) {
        fatal(error_prefix, "Option -g or --gap-opening-penalty specified when d < 2.");
      }
      if (used_options[gap_extension_penalty_index]) {
        fatal(error_prefix, "Option -e or --gap-extension-penalty specified when d < 2.");
      }
    }

  if (parameters.opt_gap_opening_penalty < 0) {
    fatal(error_prefix, "Illegal gap opening penalty specified with -g or "
          "--gap-opening-penalty, must not be negative.");
  }

  if (parameters.opt_gap_extension_penalty < 0) {
    fatal(error_prefix, "Illegal gap extension penalty specified with -e or "
          "--gap-extension-penalty, must not be negative.");
  }

  if ((parameters.opt_gap_opening_penalty + parameters.opt_gap_extension_penalty) < 1) {
    fatal(error_prefix, "Illegal gap penalties specified, the sum of the gap open and "
          "the gap extension penalty must be at least 1.");
  }

  if (parameters.opt_match_reward < 1) {
    fatal(error_prefix, "Illegal match reward specified with -m or --match-reward, "
          "must be at least 1.");
  }

  if (parameters.opt_mismatch_penalty < 1) {
    fatal(error_prefix, "Illegal mismatch penalty specified with -p or --mismatch-penalty, "
          "must be at least 1.");
  }

  if (parameters.opt_boundary < 2) {
    fatal(error_prefix, "Illegal boundary specified with -b or --boundary, "
          "must be at least 2.");
  }

  if ((used_options[ceiling_index]) and ((parameters.opt_ceiling < min_ceiling) or
                                         (parameters.opt_ceiling > max_ceiling))) {
    fatal(error_prefix, "Illegal memory ceiling specified with -c or --ceiling, "
          "must be in the range 8 to 1,073,741,824 MB.");
  }

  if ((parameters.opt_bloom_bits < min_bits_per_entry) or
      (parameters.opt_bloom_bits > max_bits_per_entry)) {
    fatal(error_prefix, "Illegal number of Bloom filter bits specified with -y or "
          "--bloom-bits, must be in the range 2 to 64.");
  }

  if ((used_options[append_abundance_index]) and (parameters.opt_append_abundance < 1)) {
    fatal(error_prefix, "Illegal abundance value specified with -a or --append-abundance, "
          "must be at least 1.");
  }

  if ((not parameters.opt_network_file.empty()) and (parameters.opt_differences != 1)) {
    fatal(error_prefix, "A network file can only written when d = 1.");
  }

  if (parameters.opt_version) {
    show(header_message);
    std::exit(EXIT_SUCCESS);
  }

  if (parameters.opt_help) {
    show(header_message);
    show(args_usage_message);
    std::exit(EXIT_SUCCESS);
  }

  // scoring system check
  const int64_t diff_saturation_16 = (std::min((uint16_max / parameters.penalty_mismatch),
                                               (uint16_max - parameters.penalty_gapopen)
                                               / parameters.penalty_gapextend));

  if (parameters.opt_differences > diff_saturation_16) {
    fatal(error_prefix, "Resolution (d) too high for the given scoring system.");
  }

  if (parameters.penalty_mismatch > uint8_max) {
    fatal(error_prefix, "Alignment scoring system yielded a mismatch penalty greater than 255, "
          "please use different parameter values.");
  }
}


auto open_files() -> void
{
  // special case (always '-')??
  outfile = fopen_output(parameters.opt_output_file.c_str());
  if (outfile == nullptr) {
    fatal(error_prefix, "Unable to open output file for writing.");
  }

  /* open files */

  if (not opt_log.empty())
    {
      logfile = fopen_output(opt_log.c_str());
      if (logfile == nullptr) {
        fatal(error_prefix, "Unable to open log file for writing.");
      }
    }

  if (not parameters.opt_seeds.empty())
    {
      fp_seeds = fopen_output(parameters.opt_seeds.c_str());
      if (fp_seeds == nullptr) {
        fatal(error_prefix, "Unable to open seeds file for writing.");
      }
    }

  if (not parameters.opt_statistics_file.empty())
    {
      statsfile = fopen_output(parameters.opt_statistics_file.c_str());
      if (statsfile == nullptr) {
        fatal(error_prefix, "Unable to open statistics file for writing.");
      }
    }

  if (not parameters.opt_uclust_file.empty())
    {
      uclustfile = fopen_output(parameters.opt_uclust_file.c_str());
      if (uclustfile == nullptr) {
        fatal(error_prefix, "Unable to open uclust file for writing.");
      }
    }

  if (not parameters.opt_internal_structure.empty())
    {
      internal_structure_file = fopen_output(parameters.opt_internal_structure.c_str());
      if (internal_structure_file == nullptr) {
        fatal(error_prefix, "Unable to open internal structure file for writing.");
      }
    }

  if (not parameters.opt_network_file.empty())
    {
      network_file = fopen_output(parameters.opt_network_file.c_str());
      if (network_file == nullptr) {
        fatal(error_prefix, "Unable to open network file for writing.");
      }
    }
}


auto close_files() -> void {
  const std::vector<std::FILE *> file_handles
    {network_file, internal_structure_file,
     uclustfile, statsfile, fp_seeds, outfile,
     logfile};
  for (auto * const file_handle : file_handles) {
    if (file_handle != nullptr) {
      std::fclose(file_handle);
    }
  }
}


auto main(int argc, char** argv) -> int
{
  // initialization and checks
  const auto used_options = args_init(argc, argv);
  args_check(used_options);
  open_files();
  show(header_message);
  args_show();

  // parse fasta input
  std::vector<char> data_v;  // refactoring: std::string fails? .data() -> const char *
  std::vector<struct seqinfo_s> seqindex_v;
  std::vector<uint64_t> zobrist_tab_base_v;
  std::vector<uint64_t> zobrist_tab_byte_base_v;
  db_read(parameters.input_filename.c_str(),
          parameters,
          data_v,
          seqindex_v,
          zobrist_tab_base_v,
          zobrist_tab_byte_base_v);

  // clustering
  switch (parameters.opt_differences)
    {
    case 0:
      dereplicate(parameters);
      break;

    case 1:
      algo_d1_run(parameters);
      break;

    default:
      algo_run(parameters, seqindex_v);
      break;
    }

  // clean up
  zobrist_exit();
  db_free();
  close_files();
}
