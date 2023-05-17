#remove the END tag in protein pdb file
rec_structure=$1
sed -i '/END/d' $rec_structure
sed -i '/CONECT/d' $rec_structure