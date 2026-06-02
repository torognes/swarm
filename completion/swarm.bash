# bash completion for swarm        -*- shell-script -*-
#
# Provides tab-completion of swarm's options and file arguments.
#
# Installation (handled automatically by `make install`):
#   copy this file to one of the bash-completion lookup directories, e.g.
#     /usr/share/bash-completion/completions/swarm
#   or source it from your ~/.bashrc:
#     source /path/to/swarm.bash

_swarm()
{
    local cur prev
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # Fallback for systems without the bash-completion package: a minimal
    # _filedir that completes files and directories.
    if ! declare -F _filedir >/dev/null 2>&1; then
        _filedir()
        {
            COMPREPLY=( $(compgen -f -- "$cur") )
        }
    fi

    # All recognised options, short and long.
    local all_opts="\
-a -b -c -d -e -f -g -h -i -j -l -m -n -o -p -r -s -t -u -v -w -x -y -z \
--append-abundance --boundary --ceiling --differences --gap-extension-penalty \
--fastidious --gap-opening-penalty --help --internal-structure --log \
--network-file --match-reward --no-otu-breaking --output-file --mismatch-penalty \
--mothur --statistics-file --threads --uclust-file --version --seeds \
--disable-sse3 --bloom-bits --usearch-abundance"

    # Options whose argument is a file name -> complete with files.
    case "$prev" in
        -i|--internal-structure|\
        -j|--network-file|\
        -l|--log|\
        -o|--output-file|\
        -s|--statistics-file|\
        -u|--uclust-file|\
        -w|--seeds)
            _filedir
            return 0
            ;;
        # Options whose argument is an integer -> no value suggestion.
        -a|--append-abundance|\
        -b|--boundary|\
        -c|--ceiling|\
        -d|--differences|\
        -e|--gap-extension-penalty|\
        -g|--gap-opening-penalty|\
        -m|--match-reward|\
        -p|--mismatch-penalty|\
        -t|--threads|\
        -y|--bloom-bits)
            return 0
            ;;
    esac

    # Complete option names when the current word starts with a dash.
    if [[ "$cur" == -* ]]; then
        COMPREPLY=( $(compgen -W "$all_opts" -- "$cur") )
        return 0
    fi

    # Otherwise complete the positional FASTA file argument with file names.
    _filedir
    return 0
}
complete -F _swarm swarm
