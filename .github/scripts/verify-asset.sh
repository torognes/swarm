#!/bin/sh
# Check one packaged swarm release asset before it is published.
#
#   usage: verify-asset.sh <assets-dir> <version>
#   ASSET       asset base name, e.g. linux-aarch64        (required)
#   MAX_GLIBC   highest glibc symbol version the binary may
#               need                                       (default: 2.34)
#
# readelf rather than ldd or objdump: it parses an ELF of any architecture
# on any host, which is what makes the cross-built aarch64 and ppc64le
# assets checkable on an x86_64 runner. The 2.34 default is the floor the
# hand-built 3.1.6 to 3.2.0 assets had; a job that silently fell out of the
# pinned container would raise it, and nobody would notice until a user on
# an older distribution reported a failure to start.
set -eu

assets_dir="${1:?usage: verify-asset.sh <assets-dir> <version>}"
version="${2:?usage: verify-asset.sh <assets-dir> <version>}"
: "${ASSET:?ASSET is required}"
MAX_GLIBC="${MAX_GLIBC:-2.34}"

dirname="swarm-${version}-${ASSET}"
workdir=$(mktemp -d)
trap 'rm -rf "${workdir}"' EXIT

if [ -f "${assets_dir}/${dirname}.zip" ] ; then
  unzip -q "${assets_dir}/${dirname}.zip" -d "${workdir}"
else
  tar xzf "${assets_dir}/${dirname}.tar.gz" -C "${workdir}"
fi
root="${workdir}/${dirname}"

status=0
fail() { echo "::error::${ASSET}: $*" >&2 ; status=1 ; }

# --- contents -------------------------------------------------------------
for file in man/swarm.1 completion/swarm.bash completion/_swarm \
            README.md CITATION.cff LICENSE ; do
  test -f "${root}/${file}" || fail "${file} is missing"
done

# The hand-packed assets carried these, because they were built on macOS;
# a Linux or macOS job here must not reintroduce them.
test "$(find "${root}" -name '._*' | wc -l)" -eq 0 || fail "AppleDouble files in the archive"

binary="${root}/bin/swarm"
test -f "${binary}" || binary="${root}/bin/swarm.exe"
test -f "${binary}" || { fail "no binary in bin/" ; exit 1 ; }

echo "${ASSET}: $(wc -c < "${binary}") byte binary"

# The version in CITATION.cff travels with the archive, so it is the one
# thing inside that can still be checked against the asset name.
grep -q "^version: *${version}\$" "${root}/CITATION.cff" ||
  fail "CITATION.cff does not say version ${version}"

# --- the binary itself ----------------------------------------------------
# Identify the format by its magic number rather than by the asset name, so
# that a matrix entry wired to the wrong compiler is caught here.
magic=$(od -An -tx1 -N4 "${binary}" | tr -d ' \n')
case "${magic}" in
  7f454c46) format=elf ;;
  cffaedfe|cefaedfe|cafebabe|bebafeca) format=macho ;;
  4d5a*) format=pe ;;
  *) format="unknown (${magic})" ;;
esac

case "${ASSET}" in
  linux-*)  expected_format=elf ;;
  macos-*)  expected_format=macho ;;
  win-*)    expected_format=pe ;;
  *)        expected_format="${format}" ;;
esac
test "${format}" = "${expected_format}" ||
  fail "expected a ${expected_format} binary and found ${format}"

case "${format}" in
  elf)
    machine=$(readelf -h "${binary}" | sed -n 's/ *Machine: *//p')
    echo "  machine: ${machine}"
    case "${ASSET}" in
      linux-x86_64)  expected_machine="X86-64" ;;
      linux-aarch64) expected_machine="AArch64" ;;
      linux-ppc64le) expected_machine="PowerPC64" ;;
      *)             expected_machine="" ;;
    esac
    case "${machine}" in
      *"${expected_machine}"*) ;;
      *) fail "expected a ${expected_machine} binary, readelf says '${machine}'" ;;
    esac

    needed=$(readelf -d "${binary}" | sed -n 's/.*(NEEDED).*\[\(.*\)\]/\1/p' | tr '\n' ' ')
    echo "  NEEDED: ${needed:-<none>}"

    floor=$(readelf -V "${binary}" | grep -o 'GLIBC_[0-9.]*' | sed 's/GLIBC_//' | sort -V | tail -1)
    echo "  needs glibc >= ${floor:-<none>}"
    if [ -n "${floor}" ] ; then
      highest=$(printf '%s\n%s\n' "${floor}" "${MAX_GLIBC}" | sort -V | tail -1)
      test "${highest}" = "${MAX_GLIBC}" ||
        fail "needs glibc ${floor}, above the ${MAX_GLIBC} floor -- was it built in the pinned container?"
    fi
    ;;
  macho)
    # lipo is only on the macOS runners, where the Mach-O assets are made.
    if command -v lipo > /dev/null 2>&1 ; then
      archs=$(lipo -archs "${binary}")
      echo "  architectures: ${archs}"
      case "${ASSET}" in
        macos-aarch64)   expected_archs="arm64" ;;
        macos-x86_64)    expected_archs="x86_64" ;;
        macos-universal) expected_archs="x86_64 arm64" ;;
        *)               expected_archs="${archs}" ;;
      esac
      # lipo prints the slices in the order they sit in the fat binary, so
      # compare the sorted sets rather than the strings.
      sort_words() { printf '%s' "$1" | tr ' ' '\n' | sort | tr '\n' ' ' ; }
      test "$(sort_words "${archs}")" = "$(sort_words "${expected_archs}")" ||
        fail "expected the architectures '${expected_archs}' and found '${archs}'"
    fi
    ;;
  pe)
    echo "  PE binary, no further checks"
    ;;
  *)
    fail "unrecognised binary format"
    ;;
esac

exit "${status}"
