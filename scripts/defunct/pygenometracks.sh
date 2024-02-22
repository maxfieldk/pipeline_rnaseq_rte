#!/bin/bash

make_tracks_file --trackFiles \
${snakemake_params[l1hs6kbintactbed]} \
${snakemake_input[bw]} \
${snakemake_params[repeatsbed]} \
-o ${snakemake_params[outputdir]}/atracks.ini 2> ${snakemake_log[0]}

#modify tracks ini file to show labels
#I set the max height for all rna tracks to 50
sed 's/labels = false/labels = true/g' ${snakemake_params[outputdir]}/atracks.ini > ${snakemake_params[outputdir]}/atracksMOD1.ini
sed 's/#overlay_previous = yes/overlay_previous = share-y/g' ${snakemake_params[outputdir]}/atracksMOD1.ini > ${snakemake_params[outputdir]}/atracksMOD2.ini
sed 's/overlay_previous = share-y/#overlay_previous = yes/1' ${snakemake_params[outputdir]}/atracksMOD2.ini > ${snakemake_params[outputdir]}/atracksMOD3.ini
sed 's/overlay_previous = share-y/#overlay_previous = yes/1' ${snakemake_params[outputdir]}/atracksMOD3.ini > ${snakemake_params[outputdir]}/atracksMOD4.ini
sed 's/#max_value = auto/max_value = 50/1' ${snakemake_params[outputdir]}/atracksMOD4.ini > ${snakemake_params[outputdir]}/atracksMOD.ini


#manually modify tracks to show labels
for telocaltype in ${snakemake_params[telocaltypes]}
do
echo $telocaltype
for contrast in ${snakemake_params[contrasts]}
do
echo $contrast
cat ${snakemake_input[DETEsbyContrast]} | while read line
do
tetype=$(awk '{print $1}' <<< $line)
echo $tetype
te=$(awk '{print $2}' <<< $line)
chr=$(awk '{print $3}' <<< $line)
start=$(awk '{print $4}' <<< $line)
stop=$(awk '{print $5}' <<< $line)
strand=$(awk '{print $6}' <<< $line)
direction=$(awk '{print $7}' <<< $line)
counttype=$(awk '{print $9}' <<< $line)
contrasttype=$(awk '{print $8}' <<< $line)
echo $chr $start $stop
echo $counttype $contrasttype
if [[ $tetype == AluY ]]
then
flanklength=200
elif [[ $tetype == L1HS ]]
then
flanklength=1000
else
flanklength=1000
fi
echo $flanklength
if [[ $counttype == $telocaltype ]] && [[ $contrasttype == $contrast ]]
then
echo we did it joe
pyGenomeTracks --tracks ${snakemake_params[outputdir]}/atracksMOD.ini --region ${chr}:$((${start}-${flanklength}))-$((${stop}+${flanklength})) \
--dpi 300 --title ${te} -o ${snakemake_params[outputdir]}/${telocaltype}/${contrast}/${te}${chr}${start}${stop}.png 2>> {log}
fi
done

done

done

touch ${snakemake_output[outfile]}