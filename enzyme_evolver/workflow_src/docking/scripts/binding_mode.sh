#this script is to combine the poses after docking with protein and save it to a file: binding_poses.pdbqt
#first argument is protein name
#second argument is the ligand name

rec_name=$1
lig_name=$2
cp $rec_name'.pdbqt' $rec_name'_'$lig_name'_binding_poses.pdbqt'
cat out.pdbqt >> $rec_name'_'$lig_name'_binding_poses.pdbqt'