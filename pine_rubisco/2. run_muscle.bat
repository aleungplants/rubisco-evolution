ECHO ON

pushd %~dp0
muscle5.1.win64.exe -align sequences.fasta -output aligned.fasta -threads 7
python "%~dp0fasta_to_phylip.py"
muscle5.1.win64.exe -align aligned.fasta -stratified -output ensemble.efa -threads 7
popd

PAUSE
