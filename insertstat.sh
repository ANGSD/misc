SAM=/usr/bin/samtools
PRG=./getinsertsize
PID=$BASHPID

TMPDIR="/scratch/"

if [ ! -d "$TMPDIR"  ]; then
    TMPDIR=/tmp/
fi

if [ ! -d "$TMPDIR"  ]; then
    echo "tmpdir: '$TMPDIR' doesnt exists"
fi

if [ "$#" -ne 1 ]
then
  echo "Usage: insertstat.sh afile.bam [onam]"
  exit 1
fi
TMPNAM=${TMPDIR}/${PID}.res
${SAM} view ${1} |${PRG} >$TMPNAM
echo "dumped tmpresults in ${TMPNAM}"
