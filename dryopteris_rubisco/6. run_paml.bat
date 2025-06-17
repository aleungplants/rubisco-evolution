ECHO ON

cd %~dp0

bash cd "$(dirname "$0")"
bash ./codeml

PAUSE
