ECHO ON

pushd %~dp0
muscle5.1.win64.exe -disperse ensemble.efa -log disperse.log
popd

PAUSE
