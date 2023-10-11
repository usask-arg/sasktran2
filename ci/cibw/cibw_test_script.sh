set -xe

PROJECT_DIR="$1"

python -c "import sasktran2 as sk; sk.appconfig.download_standard_databases()"

pytest $PROJECT_DIR/tests