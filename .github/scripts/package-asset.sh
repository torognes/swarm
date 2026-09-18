#!/bin/sh
# Package one swarm release asset from a source tree that has been built.
#
# Driven entirely by the environment, so the workflow matrix stays data and
# this stays the only copy of the packaging logic. It runs inside the build
# container (Linux), on an ordinary runner (Windows cross-build) and on
# macOS, and it assumes nothing beyond a POSIX shell, tar and zip.
#
#   SRCDIR   built swarm source tree               (default: current directory)
#   ASSET    asset base name, e.g. linux-aarch64   (required)
#   VERSION  version string, e.g. 3.2.0            (required)
#   BINARY   name of the built binary              (default: swarm)
#   FORMAT   tar | zip                             (default: tar)
#   OUTDIR   where to write the archive            (default: current directory)
#
# The layout matches what releases up to 3.2.0 shipped, with two changes
# agreed with the maintainers: the shell completions are now included (they
# are installed by "make install" but had never reached an asset), and the
# generated doc/swarm_manual.pdf is gone, leaving man/swarm.1 as the manual.
#
#   swarm-<version>-<asset>/
#   |-- bin/swarm[.exe]
#   |-- man/swarm.1
#   |-- completion/swarm.bash
#   |-- completion/_swarm
#   |-- README.md
#   |-- CITATION.cff
#   `-- LICENSE
set -eu

SRCDIR="${SRCDIR:-$(pwd)}"
: "${ASSET:?ASSET is required}"
: "${VERSION:?VERSION is required}"
BINARY="${BINARY:-swarm}"
FORMAT="${FORMAT:-tar}"
OUTDIR="${OUTDIR:-$(pwd)}"

dirname="swarm-${VERSION}-${ASSET}"

cd "${SRCDIR}"
test -f "bin/${BINARY}" || {
  echo "::error::bin/${BINARY} is missing -- was the build run?" >&2
  exit 1
}

rm -rf "${dirname:?}"
mkdir -p "${dirname}/bin" "${dirname}/man" "${dirname}/completion"
cp "bin/${BINARY}" "${dirname}/bin/"
cp man/swarm.1 "${dirname}/man/"
cp completion/swarm.bash completion/_swarm "${dirname}/completion/"
cp README.md CITATION.cff LICENSE "${dirname}/"

mkdir -p "${OUTDIR}"
if [ "${FORMAT}" = "zip" ] ; then
  rm -f "${OUTDIR}/${dirname}.zip"
  zip -qr "${OUTDIR}/${dirname}.zip" "${dirname}"
else
  # --format=ustar keeps the archive free of the vendor extensions GNU tar
  # and bsdtar each add by default, and COPYFILE_DISABLE stops bsdtar on
  # macOS slipping AppleDouble "._swarm" members and com.apple.* extended
  # attributes into it -- both of which the hand-packed assets carried, so
  # that GNU tar warned on every extraction of a Linux tarball.
  set -- --format=ustar
  if tar --version 2>/dev/null | grep -q GNU ; then
    set -- "$@" --owner=0 --group=0 --numeric-owner
  fi
  rm -f "${OUTDIR}/${dirname}.tar.gz"
  COPYFILE_DISABLE=1 tar "$@" -czf "${OUTDIR}/${dirname}.tar.gz" "${dirname}"
fi

rm -rf "${dirname:?}"
ls -l "${OUTDIR}/${dirname}".*
