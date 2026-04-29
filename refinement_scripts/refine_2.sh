#!/bin/bash

#this script will run a refinement based on the input PDB names you have.

#________________________________________________INPUTS________________________________________________#
base_dir=''
PDB_file=${base_dir}/pxr_ids.txt

#________________________________________________SET PATHS________________________________________________#
source ~/phenix-2.0-5824/phenix_env.sh
export PHENIX_OVERWRITE_ALL=true

#________________________________________________Read from input csv________________________________________________#
PDB=$(cat $PDB_file | head -n $SGE_TASK_ID | tail -n 1)

#This section extracts the required information in a way that can accomodate multiple ligands/ligand types
group_csv() {
    local line="$1"
    local varA="$2"
    local varB="$3"
    local varC="$4"
    local varD="$5"

    # Split into array
    IFS=',' read -ra fields <<< "$line"

    # First entry is pdb_id
    pdb_id="${fields[0]}"

    # Initialize empty strings
    local A=""
    local B=""
    local C=""
    local D=""

    # Walk through rest of fields (start at index 1)
    for ((i=1; i<${#fields[@]}; i++)); do
        pos=$(( (i-1) % 4 ))   # cycle: 0 -> A, 1 -> B, 2 -> C
        echo "${fields[i]}"
        case $pos in
            0) A+="${A:+,}${fields[i]}" ;;
            1) B+="${B:+,}${fields[i]}" ;;
            2) C+="${C:+,}${fields[i]}" ;;
            3) D+="${D:+,}${fields[i]}" ;;
        esac
    done

    # Export into user-named variables
    printf -v "$varA" "%s" "$A"
    printf -v "$varB" "%s" "$B"
    printf -v "$varC" "%s" "$C"
    printf -v "$varD" "%s" "$D"
}


group_csv "${PDB}" ligand_name chain res smiles
echo "pdb name: ${pdb_id}"
echo "ligand(s): ${ligand_name}"
echo "chain(s) each ligand is found in: ${chain}"
echo "residue each ligand is found in: ${res}"
echo "smiles string: ${smiles}"

cd ${base_dir}/${pdb_id}
echo "Working in directory ${base_dir}/${pdb_id}"

#________________________________________________Refinement Preprocessing________________________________________________#

# Convert cif file to mtz
phenix.cif_as_mtz ${pdb_id}-sf.cif --merge
echo "${pdb_id}-sf.cif merging into mtz"
if [ ! -f "${pdb_id}-sf.mtz" ]; then
    echo "merge failed"
fi

# Move the converted file to the new mtz. If conversion fails, script uses original mtz
echo "moving cif file"
mv ${pdb_id}-sf.mtz ${pdb_id}.mtz
echo "cif file moved "

# Determine FOBS VS IOBS VS FP
mtzmetadata=`phenix.mtz.dump "${pdb_id}.mtz"`

obstypes=("FP" "FOBS" "F-obs" "I" "IOBS" "I-obs" "F(+)" "I(+)" "FSIM")
ampfields=`grep -E "amplitude|intensity|F\(\+\)|I\(\+\)" <<< "${mtzmetadata}"`

# Clear xray_data_labels variable
xray_data_labels=""

# Is amplitude an Fo?
for field in ${ampfields}; do
  # Check field in obstypes
  if [[ " ${obstypes[*]} " =~ " ${field} " ]]; then
    # Check SIGFo is in the mtz too!
    if grep -F -q -w "SIG$field" <<< "${mtzmetadata}"; then
      xray_data_labels="${field},SIG${field}";
      break
    fi

    ##This if statement is to catch a specific edge case where the labels are IOBS,SIGI
    if [[ "$field" == "IOBS" ]]; then
      if grep -F -q -w "SIGI" <<< "${mtzmetadata}"; then
        xray_data_labels="${field}";
        break
      fi
    fi
  fi
done
if [ -z "${xray_data_labels}" ]; then
  echo >&2 "Could not determine Fo field name with corresponding SIGFo in .mtz.";
  echo >&2 "Was not among "${obstypes[*]}". Please check .mtz file\!";
  exit 1;
else
  echo "data labels: ${xray_data_labels}"
fi

#________________________________________________Run ready_set and refinement________________________________________________#

# Run ready_set
echo "Running ready set in direcory $PWD"
phenix.ready_set ${pdb_id}.pdb

# Check if the .ligands.cif file exists
if [ -f "${pdb_id}.ligands.cif" ]; then
  echo "ligand restraints not found in distribution"

  #Since base restraints do not work well for complex ligands, run elbow seperately with mopac 
  ###We have to run elbow seperately for each ligand
  IFS=',' read -ra ligs <<< "$ligand_name"
  IFS=',' read -ra sms  <<< "$smiles"

  > ${pdb_id}_ligands.cif

  for i in "${!ligs[@]}"; do
    lig="${ligs[i]}"
    smi="${sms[i]}"
    echo "running elbow for ${lig}-${smi}"

  phenix.elbow --template "${pdb_id}.updated.pdb" \
              --residue "$lig" \
              --smiles "$smi" \
              --mopac

    cat "elbow.${lig}.cif" >> ${pdb_id}_ligands.cif
  done

  if [[ ! -s "${pdb_id}_ligands.cif" ]]; then
    echo "Error: ${pdb_id}_ligands.cif is empty or does not exist"
    echo "elbow failed. exiting"
    exit 1
  fi

  #run refinement with elbow restraints
  phenix.refine ${pdb_id}.mtz ${pdb_id}.updated.pdb \
  refinement.input.monomers.file_name=${pdb_id}_ligands.cif \
  refinement.refine.strategy=individual_sites+individual_adp+occupancies \
  output.prefix=${pdb_id} \
  refinement.main.number_of_macro_cycles=5 \
  refinement.main.nqh_flips=True \
  refinement.output.write_maps=False \
  refinement.hydrogens.refine=riding \
  refinement.main.ordered_solvent=True \
  refinement.target_weights.optimize_xyz_weight=true \
  refinement.target_weights.optimize_adp_weight=true \
  data_manager.fmodel.xray_data.r_free_flags.generate=true \
  data_manager.miller_array.labels.name="${xray_data_labels}"

else
  echo "ligand restraints found in distribution"
  #run with in distribution restraints

  phenix.refine ${pdb_id}.mtz ${pdb_id}.updated.pdb \
  refinement.refine.strategy=individual_sites+individual_adp+occupancies \
  output.prefix=${pdb_id} \
  refinement.main.number_of_macro_cycles=5 \
  refinement.main.nqh_flips=True \
  refinement.output.write_maps=False \
  refinement.hydrogens.refine=riding \
  refinement.main.ordered_solvent=True \
  refinement.target_weights.optimize_xyz_weight=true \
  refinement.target_weights.optimize_adp_weight=true \
  data_manager.fmodel.xray_data.r_free_flags.generate=true \
  data_manager.miller_array.labels.name="${xray_data_labels}"
fi


