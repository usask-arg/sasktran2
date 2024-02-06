set -xe

PROJECT_DIR="$1"

pytest -W 'once:The distutils package is deprecated:DeprecationWarning' $PROJECT_DIR/tests
