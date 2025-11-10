#!/bin/bash

module purge 
module load gromacs/2025.2


cd ../1_PDBs || exit

run_num=3

for file in *.pdb; do
    sequence_name="${file%.pdb}"

    echo "Working on simulation based on $sequence_name"

    CONFIG_FILE="../4_ConfigFiles/${sequence_name}.config"
    if [[ -f "$CONFIG_FILE" ]]; then
        source "$CONFIG_FILE"
    else
        echo "Config for $sequence_name not found at $CONFIG_FILE"
        exit 1
    fi

    for force in "${forces[@]}"; do
        force_label="${sequence_name}_${force}nN"
        scaled_force=$(awk "BEGIN { printf \"%.0f\", 602 * $force }")
        current_folder="../5_Runs/${force_label}"
        mkdir -p "$current_folder"
        for (( run=1; run<=run_num; run++ )); do
            run_dir="${current_folder}/R${run}"
            mkdir -p "$run_dir"
            cd "$run_dir" || continue

            cp "../../../1_PDBs/$file" .
            cp -r ../../../3_InputFiles/* .

            sed -i "s|^#SBATCH --job-name=.*|#SBATCH --job-name=${force_label}_R${run}|" resubmit_phase_full_node.sh
            sed -i "s|^#SBATCH --job-name=.*|#SBATCH --job-name=${force_label}_R${run}|" resubmit_phase_part_node.sh

            sed -i "s|^pull-coord1-k.*|pull-coord1-k            = ${scaled_force} ; ${force}nN|" pull.mdp
            sed -i "s|^pull-coord2-k.*|pull-coord2-k            = -${scaled_force} ; ${force}nN|" pull.mdp

	          echo "Starting pdb2gmx"

            echo -e "1\n1\n" | gmx pdb2gmx -f "${sequence_name}.pdb" \
                                            -o "${sequence_name}.gro" \
                                            -p "${sequence_name}.top" \
                                            -i "${sequence_name}.itp"

            gmx editconf -f "${sequence_name}.gro" \
                         -o "${sequence_name}_newbox.gro" \
                         -rotate "$rot_x" "$rot_y" "$rot_z" \
                         -c \
                         -box "$box_x" "$box_y" "$box_z"

            gmx solvate -cp "${sequence_name}_newbox.gro" \
                        -cs spc216.gro \
                        -o "${sequence_name}_solvated.gro" \
                        -p "${sequence_name}.top"

            gmx grompp -f ions.mdp \
                       -c "${sequence_name}_solvated.gro" \
                       -p "${sequence_name}.top" \
                       -o ions.tpr

            echo -e "3\n" | gmx genion -s ions.tpr \
                                      -o "${sequence_name}_ions.gro" \
                                      -p "${sequence_name}.top" \
                                      -pname NA \
                                      -nname CL \
                                      -neutral

            gmx grompp -f minim.mdp \
                       -c "${sequence_name}_ions.gro" \
                       -p "${sequence_name}.top" \
                       -o "${sequence_name}_EM.tpr" \
                       -maxwarn 1

            gmx mdrun -v -deffnm "${sequence_name}_EM" -nt 72 -pin on

            # Create index file for pulling
            ndx_cmd_file="make_ndx_commands.txt"
            echo "ri $PullLeft1 | ri $PullLeft2" > "$ndx_cmd_file"
            echo "ri $PullRight1 | ri $PullRight2" >> "$ndx_cmd_file"
            echo "q" >> "$ndx_cmd_file"
            gmx make_ndx -f "${sequence_name}_EM.gro" -o "${sequence_name}_EM.ndx" < "$ndx_cmd_file"
            left_pattern="\[ r_${PullLeft1}_r_${PullLeft2} \]"
            right_pattern="\[ r_${PullRight1}_r_${PullRight2} \]"
            sed -i "s/${left_pattern}/[ PullLeft ]/" "${sequence_name}_EM.ndx"
            sed -i "s/${right_pattern}/[ PullRight ]/" "${sequence_name}_EM.ndx"


            mv "${sequence_name}.top" "${sequence_name}_EM.top"
            rm -f "#"* step*
            rm -f "$ndx_cmd_file"

            cd ../../../1_PDBs || exit
        done
    done
done
