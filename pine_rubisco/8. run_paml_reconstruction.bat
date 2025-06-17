ECHO ON

cd %~dp0
copy "clean_aligned_paml.phylip" "paml_reconstruction\clean_aligned_paml.phylip"
copy "Pinus_time_paml.tre" "paml_reconstruction\Pinus_time_paml.tre"
cd paml_reconstruction
"codeml.exe"

PAUSE
