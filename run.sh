#!/bin/bash

POSITIONAL=()
HELP="NO"

for i in "$@"
    do
    case $i in
        -s=*|--start=*)
        START_DATE="${i#*=}"
        shift
        ;;
        -e=*|--end=*)
        END_DATE="${i#*=}"
        shift
        ;;
        -d=*|--data=*)
        DATADIR="${i#*=}"
        shift
        ;;
        -o=*|--out=*)
        OUTDIR="${i#*=}"
        shift
        ;;
        -h|--help)
        HELP="YES"
        shift
        ;;
        --default)
        shift
        ;;
        *)
        ;;
    esac
    done

DATADIR=${DATADIR:-<default_path>}
OUTDIR=${OUTDIR:-~<default_path>}
WORKDIR=`pwd`

export RADASREF=<path_to_script>
export PATH=$RADASREF:$PATH

START_DATE=${START_DATE:-1995-06-28}
END_DATE=${END_DATE:-2011-07-02}

echo " " > $OUTDIR/std_output_file.txt
i="$START_DATE"
while [ "$i" != "$(date -I -d "$END_DATE + 1 day")" ];
    do
    DIRDATE="$(date +%Y/%m/%d -d "$i" )"
    FILDATE="$(date +%Y%m%d -d "$i" )"
    if [ -d ${DATADIR}/${DIRDATE} ]
        then
        L1files=`ls ${DATADIR}/${DIRDATE}/*.nc`
        if [ ! -f ${OUTDIR}/ER2_RPRO_GOM_L1B_rad_average_HCHO_${FILDATE}.nc ]
            then
            echo "Today is ${i}" >> $OUTDIR/std_output_file.txt
            begin_day=$(date +%s)
            echo $L1files >> $OUTDIR/std_output_file.txt
            RadAsRef.py -i $L1files -o ${OUTDIR}/ER2_RPRO_GOM_L1B_rad_average_HCHO_${FILDATE}.nc --gome --apex --wlen 289.59 405.18 0.01 --lat -15 15 --lon 150 250 --plot > $OUTDIR/average_${FILDATE}.txt
            end_day=$(date +%s)
            echo "It took me $( expr $end_day - $begin_day ) seconds for ${i}" >> $OUTDIR/output_dalibor.txt
        fi
    fi
    i=$(date -I -d "$i + 1 day")
    done
