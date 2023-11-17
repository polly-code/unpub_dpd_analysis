#!bin/bash

## Author(s): Pavel Kos
## Contact: pavel.kos@fmi.ch
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

GREEN=$(tput setaf 2)
RED=$(tput setaf 1)
NC=$(tput sgr0)
temp="temp_$(uuidgen).dat"
echo "Temporary file $temp"
echo Start >> $temp
echo $1 >> $temp

if [ "$2" == "replot" ]; then
replot=true
else
replot=false
fi

if [ "$3" == "difflen" ]; then
    difflen=true
else
    difflen=false
fi

if [ "$2" == "difflen" ]; then
    difflen=true
else
    difflen=false
fi

echo -en "\nFolder path to be analyzed: $1\n"
# Check whether files are here
if [ ! -f $1/mytraj.lammpstrj ]; then
echo -en "${RED}[-][-]${NC} No traj file\n"
return 1
fi
if [ ! -f $1/dsts.dat ]; then
echo -en "${RED}[-][-]${NC} No file with distances\n"
return 1
fi
if [ ! -f $1/rgyr.dat ]; then
echo -en "${RED}[-][-]${NC} No file with radius of gyration\n"
return 1
fi

# Check whether python env is the proper one
if [ "$CONDA_PREFIX" != "/tungstenfs/scratch/ggiorget/pavel/code_projects/multiple_trajs/analysis_chrom_dynamics/andyn" ]; then
echo -en "${RED}[-]${NC} Environment is wrong"
conda activate ../multiple_trajs/analysis_chrom_dynamics/andyn/ >> $temp 2>&1
echo -en "\033[2K"
echo -en "\r${GREEN}[+]${NC} Environment is the proper one\n"
else
echo -en "${GREEN}[+]${NC} Environment is the proper one\n"
fi
wait

# Check wether restart generated obsolete data
if [ "$difflen" = false ] ; then
    echo -en "${RED}[-]${NC} Checking/adjusting files because of restart"
    . code/remove_overlapping_due_to_restart.sh $1 >> $temp 2>&1
    if [ $? -eq 0 ]; then
        echo -en "\033[2K"
        echo -en "\r${GREEN}[+]${NC} Data adjusted\n"
    else
        echo -en "\033[2K"
        echo -en "\r${RED}[-][-]${NC} Something went wrong with traj/dsts/rgyr files\n"
        return 1
    fi
fi

# Check whether motion correction has been already done
if [ ! -f $1/mytraj_c.lammpstrj ] || [ "$replot" = true ] ; then
    echo -en "${RED}[-]${NC} Subtracting motion of entire system"
    python3 code/motion_correction/lmp_motion_correction.py $1/mytraj.lammpstrj code/motion_correction/mychain.mol2 >> $temp 2>&1
    if [ -f $1/mytraj_c.lammpstrj ]; then
        echo -en "\033[2K"
        echo -en "\r${GREEN}[+]${NC} Motion corrected\n"
    else
        echo -en "\033[2K"
        echo -en "\r${RED}[-][-]${NC} Something went wrong with substraction\n"
        return 1
    fi
else
    echo -en "${GREEN}[+]${NC} Motion corrected: file exists\n"
fi
wait

# Check whether maps have been already calculated 
if ([ ! -f $1/mytraj_c.lammpstrj_dm.dat ] && [ ! -f $1/mytraj_c.lammpstrj_cm.dat ]) || [ "$replot" = true ] ; then
echo -en "${RED}[-]${NC} Calculating contact/distance maps"
mpiexec -np 10 code/distance_matrix/cmdm -p $1/mytraj_c.lammpstrj >> $temp 2>&1
if [ -f $1/mytraj_c.lammpstrj_dm.dat ] && [ -f $1/mytraj_c.lammpstrj_cm.dat ]; then
echo -en "\033[2K"
echo -en "\r${GREEN}[+]${NC} Contact and distance matrices exist\n"
else
echo -en "\033[2K"
echo -en "\r${RED}[-][-]${NC} Something went wrong with calculating contact/distance maps\n"
return 1
fi
else
echo -en "${GREEN}[+]${NC} Contact and distance matrices exist\n"
fi
wait

# Check whether msd single has been calculated
if [ ! -f $1/mytraj_c.lammpstrj_msd_s.dat ] || [ "$replot" = true ]; then
echo -en "${RED}[-]${NC} Calculating MSD of single beads"
mpiexec -np 10 code/msd_c/msd -m single -lmp -p $1/mytraj_c.lammpstrj >> $temp 2>&1
if [ -f $1/mytraj_c.lammpstrj_msd_s.dat ]; then
echo -en "\033[2K"
echo -en "\r${GREEN}[+]${NC} MSD per bead (single) has been calculated\n"
else
echo -en "\033[2K"
echo -en "\r${RED}[-][-]${NC} Something went wrong with calculating MSD of single beads\n"
return 1
fi
else
echo -en "${GREEN}[+]${NC} MSD per bead (single) has been calculated\n"
fi

# Check whether msd average has been calculated
if [ ! -f $1/mytraj_c.lammpstrj_msd_a.dat ] || [ "$replot" = true ] ; then
echo -en "${RED}[-]${NC} Calculating averaged MSD"
mpiexec -np 10 code/msd_c/msd -m aver -lmp -p $1/mytraj_c.lammpstrj >> $temp 2>&1
if [ -f $1/mytraj_c.lammpstrj_msd_a.dat ]; then
echo -en "\033[2K"
echo -en "\r${GREEN}[+]${NC} Average MSD has been calculated\n"
else
echo -en "\033[2K"
echo -en "\r${RED}[-][-]${NC} Something went wrong with calculating averaged MSD\n"
return 1
fi
else
echo -en "${GREEN}[+]${NC} Average MSD has already been calculated\n"
fi

# Check whether distance map has been already averaged to 1k map
if [ ! -f $1/mytraj_c.lammpstrj_dm.dat_s ] || [ "$replot" = true ] ; then
echo -en "${RED}[-]${NC} Averaging distance map"
python3 code/10kto1k.py $1/mytraj_c.lammpstrj_dm.dat aver >> $temp 2>&1
if [ -f $1/mytraj_c.lammpstrj_dm.dat_s ]; then
echo -en "\033[2K"
echo -en "\r${GREEN}[+]${NC} Distance map averaged\n"
else
echo -en "\033[2K"
echo -en "\r${RED}[-][-]${NC} Something went wrong with averaging distance map\n"
return 1
fi
else
echo -en "${GREEN}[+]${NC} Distance map has already been averaged to 1k map\n"
fi

# Check whether contact map has been already averaged to 1k map
if [ ! -f $1/mytraj_c.lammpstrj_cm.dat_s ] || [ "$replot" = true ] ; then
echo -en "${RED}[-]${NC} Averaging contact map"
python3 code/10kto1k.py $1/mytraj_c.lammpstrj_cm.dat >> $temp 2>&1
if [ -f $1/mytraj_c.lammpstrj_cm.dat_s ]; then
echo -en "\033[2K"
echo -en "\r${GREEN}[+]${NC} Contact map is averaged\n"
else
echo -en "\033[2K"
echo -en "\r${RED}[-][-]${NC} Something went wrong with averaging contact map\n"
return 1
fi
echo -en "${GREEN}[+]${NC} Contact map has already been averaged to 1k map\n"
fi

# Check whether distance map has been already plotted
if [ ! -f $1/mytraj_c.lammpstrj_dm.dat_s.png ] || [ "$replot" = true ] ; then
echo -en "${RED}[-]${NC} Plotting distance map"
python3 code/dmap.py $1/mytraj_c.lammpstrj_dm.dat_s 1000 >> $temp 2>&1
if [ -f $1/mytraj_c.lammpstrj_dm.dat_s.png ]; then
echo -en "\033[2K"
echo -en "\r${GREEN}[+]${NC} Distance map is plotted\n"
else
echo -en "\033[2K"
echo -en "\r${RED}[-][-]${NC} Something went wrong with plotting distance map\n"
return 1
fi
else
echo -en "${GREEN}[+]${NC} Distance map has already been plotted\n"
fi

# Check whether contact map has been already plotted
if [ ! -f $1/mytraj_c.lammpstrj_cm.dat_s.png ] || [ "$replot" = true ] ; then
echo -en "${RED}[-]${NC} Plotting contact map"
python3 code/cmap.py $1/mytraj_c.lammpstrj_cm.dat_s >> $temp 2>&1
if [ -f $1/mytraj_c.lammpstrj_cm.dat_s.png ]; then
echo -en "\033[2K"
echo -en "\r${GREEN}[+]${NC} Contact map is plotted\n"
else
echo -en "\033[2K"
echo -en "\r${RED}[-][-]${NC} Something went wrong with plotting contact map\n"
return 1
fi
else
echo -en "${GREEN}[+]${NC} Contact map has already been plotted\n"
fi

# Check whether plots from dsts.dat are done
if [ ! -f $1/hist_out.png ] || [ ! -f $1/hist_in.png ] || [ ! -f $1/hist_across.png ] || [ ! -f $1/hist_ctcf.png ] || [ ! -f $1/hist_in_out_across.png ] || [ "$replot" = true ] ; then
echo -en "${RED}[-]${NC} Plotting hists of distances"
python3 code/hist_distances.py $1  >> $temp 2>&1
if [ -f $1/hist_out.png ] && [ -f $1/hist_in.png ] && [ -f $1/hist_across.png ] && [ -f $1/hist_ctcf.png ] && [ -f $1/hist_in_out_across.png ]; then
echo -en "\033[2K"
echo -en "\r${GREEN}[+]${NC} Hists are plotted\n"
else
echo -en "\033[2K"
echo -en "\r${RED}[-][-]${NC} Something went wrong with plotting hists of distances\n"
return 1
fi
else
echo -en "${GREEN}[+]${NC} Hists have already been plotted\n"
fi

# Check whether rgyr has been already plotted
if [ ! -f $1/hist_rgyr_in_out_across.png ] || [ ! -f $1/bp_rgyr_in_out_across.png ] || [ "$replot" = true ] ; then
echo -en "${RED}[-]${NC} Plotting rgyr distributions"
python3 code/hist_rgyr.py $1 >> $temp 2>&1
if [ -f $1/hist_rgyr_in_out_across.png ] && [ -f $1/bp_rgyr_in_out_across.png ]; then
echo -en "\033[2K"
echo -en "\r${GREEN}[+]${NC} Rgyr distributions are plotted\n"
else
echo -en "\033[2K"
echo -en "\r${RED}[-][-]${NC} Something went wrong with plotting rgyr distributions\n"
return 1
fi
else
echo -en "${GREEN}[+]${NC} Rgyr distributions have already been plotted\n"
fi

# Check whether pc has been already plotted
if [ ! -f $1/pc_contact_prob.png ] || [ "$replot" = true ] ; then
echo -en "${RED}[-]${NC} Plotting contact probability"
python3 code/contact_probability.py $1 >> $temp 2>&1
if [ -f $1/pc_contact_prob.png ]; then
echo -en "\033[2K"
echo -en "\r${GREEN}[+]${NC} Pc is plotted\n"
else
echo -en "\033[2K"
echo -en "\r${RED}[-][-]${NC} Something went wrong with plotting Pc\n"
return 1
fi
else
echo -en "${GREEN}[+]${NC} Pc has already been plotted\n"
fi

echo -en "${GREEN}[+][+]${NC} Done!\n"
rm $temp
return 0
