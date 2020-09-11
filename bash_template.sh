# Author: Marco Kreuzer <marco.kreuzer@bioinformatics.unibe.ch>
#
#/ Usage: SCRIPTNAME [OPTIONS]... [ARGUMENTS]...
#/ 
#/ OPTIONS
#/
#/ EXAMPLES
#/  

#{{{ CL arguments
while getopts ":t:" opt; do
  case ${opt} in
    t )
	  target=$OPTARG
	  ;;
	\? )
      echo "Invalid option: $OPTARG" 1>&2
	  ;;
	: )
	  echo "Invalid option: $OPTARG requires an argument" 1>&2
	  ;;
  esac
done
shift $((OPTIND -1))
#}}}

# Bash settings

set -o errexit # abort on nonzero exitstatus
set -o nounset # abort on unbound variable     
set -o pipefail # dont hide errors within pipes


#{{{ Variables
readonly script_name=$(basename "${0}")
readonly script_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
IFS=$'\t\n'   # Split on newlines and tabs (but not on spaces)
#}}}

main() {
	
}
#{{{ Helper functions

myfunc(){
}

#}}}

main "${@}"
