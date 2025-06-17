ECHO ON

pushd %~dp0
muscle5.1.win64.exe -align clean_aligned.fasta -stratified -output ensemble.efa -threads 16
muscle5.1.win64.exe -disperse ensemble.efa -log disperse.log
popd

PAUSE
