#!/bin/bash
traj="$1/mytraj.lammpstrj"
ref_size_traj=50055009
dsts="$1/dsts.dat"
ref_size_dsts=5000002
rgyr="$1/rgyr.dat"
ref_size_rgyr=5000002

echo "Start analyzing traj file"
fact_num_lines_traj=$(wc -l ${traj} | awk '{print $1}')
if [ "$ref_size_traj" != "$fact_num_lines_traj" ]
then
    steps=( $(grep -A1 "ITEM: TIMESTEP" ${traj} | grep 0000) ) #awk -F ":" '{print $1}') )
    duplicates=( $(printf '%s\n' "${steps[@]}" | awk '!($0 in seen){seen[$0];next} 1') )
    lims=( $(grep -n ${duplicates[0]} mytraj.lammpstrj | awk -F ":" '{print $1}') )
    start_del_traj=${lims[0]}
    let end_del_traj=${lims[1]}-1
    if [ -z ${start_del_traj} ] || [ -z ${end_del_traj} ]
    then
	break
    fi
    sed -i ${start_del_traj},${end_del_traj}d ${traj}
fi

echo "Traj is adjusted"
echo "Start analyzing dsts.dat"

fact_num_lines_dsts=$(wc -l ${dsts} | awk '{print $1}')
if [ "$ref_size_dsts" != "$fact_num_lines_dsts" ]
then
    let start_del_dst=${duplicates[0]}/10+2
    let end_del_dst=${duplicates[0]}/10+$fact_num_lines_dsts-$ref_size_dsts+1
    if [ -z ${start_del_dst} ] || [ -z ${end_del_dst} ]
    then
        break
    fi
    sed -i ${start_del_dst},${end_del_dst}d ${dsts}
fi

echo "dsts.dat is adjusted"
echo "Start analyzing rgyr"

fact_num_lines_rgyr=$(wc -l ${rgyr} | awk '{print $1}')
if [ "$ref_size_rgyr" != "$fact_num_lines_rgyr" ]
then
    let start_del_rgyr=${duplicates[0]}/10+2
    let end_del_rgyr=${duplicates[0]}/10+$fact_num_lines_rgyr-$ref_size_rgyr+1
    if [ -z ${start_del_rgyr} ] || [ -z ${end_del_rgyr} ]
    then
        break
    fi
    sed -i ${start_del_rgyr},${end_del_rgyr}d ${rgyr}
fi
echo "Unnecessary data removed"

