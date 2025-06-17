ECHO ON

pushd %~dp0
muscle5.1.win64.exe -align clean_spinach.fasta -output clean_aligned.fasta -threads 16
popd

PAUSE
