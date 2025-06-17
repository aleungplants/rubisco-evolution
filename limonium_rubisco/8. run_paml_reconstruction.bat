ECHO ON

cd %~dp0
copy "clean_aligned_paml.phylip" "paml_reconstruction\clean_aligned_paml.phylip"
copy "koutroumpa_tree_paml.tre" "paml_reconstruction\koutroumpa_tree_paml.tre"
cd paml_reconstruction
"codeml.exe"

PAUSE
