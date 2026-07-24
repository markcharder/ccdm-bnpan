#!/bin/bash --login

### Common functions for workflow submissions ###

print_nice() {
    for msg in "$@"; do
        printf '%-30s\n' "$msg"
    done
}

parse_args_generic() {

    # Read argument definitions as a string
    # that is a name reference (-n) to an associative array.
    local -n arg_definitions="$1"
    shift

    # Declare 3 associative arrays to hold arg mappings.
    declare -A args_map
    declare -A required_args
    declare -A default_values

    # Build arbitrary usage and help text.
    # Start as minimal strings.
    local usage_text="Usage: $0"
    local help_text="\nOptions:\n"

    # For each key in arg definitions.
    for key in "${!arg_definitions[@]}"; do

        # Read the input key with a ':' separator.
        IFS='|' read -r long_flag short_flag requirement restk <<< "$key"
        IFS='|' read -r description default_val restv <<< "${arg_definitions[$key]}"

        # Make sure a colon has not been used in key or value
        # strings as it will break the parsing.
        if [[ -n "$restk" ]]; then
            print_nice "Error: Invalid key definition: ${key}." \
                       "Make sure that colons (:) are not used in argument flag names or requirements."
            exit 1
        fi

        if [[ -n "$restv" ]]; then
            print_nice "Error: Invalid value definition for argument: ${arg_definitions[$key]}." \
                       "Make sure that colons (:) are not used in argument descriptions or default values."
            exit 1
        fi

        # Add '--' to long flags.
        args_map["--${long_flag}"]="$long_flag"

        # If there is a short flag, add '-' to it.
        if [[ -n "$short_flag" ]]; then
            args_map["-${short_flag}"]="$long_flag"
        fi

        if [[ "$requirement" == "required" ]]; then
            required_args["$long_flag"]=1
            usage_text+=" --${long_flag} <value>"
        else
            usage_text+=" [--${long_flag} <value>]"
            if [[ -n "$default_val" ]]; then
                default_values["$long_flag"]="$default_val"
            fi
        fi

        # Build help text
        local flag_text="  "
        [[ -n "$short_flag" ]] && flag_text+="-${short_flag}, "
        flag_text+="--${long_flag}"
        flag_text=$(print_nice "$flag_text")
        help_text+="${flag_text}${description}"
        [[ -n "$default_val" ]] && help_text+=" (default: ${default_val})"
        help_text+="\n"
    done

    help_text+="$(printf '%-30s' '  -h, --help')Show this help message\n"

    # Initialize with defaults.
    for arg_name in "${!default_values[@]}"; do
        export "${arg_name//-/_}"="${default_values[$arg_name]}"
    done

    if [[ ! -z "$CUSTOM_USAGE" ]]; then
        usage_text="$CUSTOM_USAGE"
        help_text=""
    fi

    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                echo -e "$usage_text"
                echo -e "$help_text"
                exit 0
                ;;
            -*|--*)
                local flag="$1"
                if [[ -z "${args_map[$flag]}" ]]; then
                    echo "Error: Unknown option: $flag"
                    echo -e "$usage_text"
                    exit 1
                fi

                local var_name="${args_map[$flag]}"
                if [[ $# -lt 2 ]]; then
                    echo "Error: Option $flag requires a value."
                    exit 1
                fi

                # Export variable (replace hyphens with underscores for valid bash variables).
                export "${var_name//-/_}"="$2"
                unset required_args["$var_name"]
                shift 2
                ;;
            *)
                echo "Error: Unknown argument: $1"
                echo -e "$usage_text"
                exit 1
                ;;
        esac
    done

    # Check required arguments
    if [[ ${#required_args[@]} -gt 0 ]]; then
        echo "Error: Missing required arguments: ${!required_args[*]}"
        echo -e "$usage_text"
        exit 1
    fi
}

monitor_tmp () {

    local pid="$1"
    local tmp_dir="$2"
    local sleep_interval="$3"
    local nlines="$4"
    shift 4
    local files_monitor=("$@")
    local tmp_size=$(df -h "$tmp_dir" | tail -1 | awk '{print $3}' | sed 's/[A-Za-z]//g')
    local tmp_peak=0

    while kill -0 "$pid" 2> /dev/null; do
        tmp_size=$(df "$tmp_dir" | tail -1 | awk '{print $3}')

        if (( $(echo "$tmp_size > $tmp_peak" | bc -l) )); then
            tmp_peak="$tmp_size"
        fi

       if [[ "${#files_monitor[@]}" -gt 0 ]]; then

            for file in "${files_monitor[@]}"; do

                if [ -f "$file" ]; then
                    tail -n "$nlines" "$file"
                fi

            done

        fi

        sleep "$sleep_interval"

    done

    echo -e "\nPeak $tmp_dir usage: ${tmp_peak}\n"

    wait "$pid"

}

usage_info() {
    echo -e "\nCommon Functions Script"
    echo -e "\nOptions:"
    echo -e "  -n, --nxf-opts            Specify additional Nextflow options"
    echo -e "  -s, --skip-singularity    Skip pulling Singularity containers"
    echo -e "  -w, --workflow-root       Specify the workflow root directory"
    echo -e "  -h, --help                Show this help message and exit\n"
}

parse_args() {

    NXF_OPTS=""
    SINGULARITY_PULL=true

    while [[ $# -gt 0 ]]; do

        case "$1" in
            -n|--nxf-opts)
                NXF_OPTS="$2"
                shift 2
                ;;
            -s|--skip-singularity)
                SINGULARITY_PULL=false
                shift
                ;;
            -w|--workflow-root)
                WF_ROOT="$2"
                shift 2
                ;;
            -h|--help)
                usage_info
                exit 0
                ;;
            *)
                echo "Unknown option: $1"
                exit 1
                ;;
        esac

    done

}

load_modules() {

    modules=("$@")
    for module in "${modules[@]}"; do
        module load "$module"
    done

}

singularity_setup() {

    load_modules singularity/4.1.0-nompi
    unset SBATCH_EXPORT
    export SINGULARITY_CACHEDIR="$MYSCRATCH/.singularity"
    mkdir -p "$SINGULARITY_CACHEDIR"

}

singularity_pull() {

    singularity_setup

    if [ "$SINGULARITY_PULL" = false ] ; then
        echo -e "\nSkipping Singularity container pulls as requested.\n"
        return
    fi

    containers=("$@")
    CURRENT_DIR=`pwd`
    cd "$SINGULARITY_CACHEDIR" || exit 1

    echo -e "\nPulling Singularity containers:\n\n"

    for container in "${containers[@]}"; do
        echo "  - $container"
    done

    echo ""

    for container in "${containers[@]}"; do
        echo "Pulling $container..."
        singularity pull --force "$container"
    done

    cd "$CURRENT_DIR" || exit 1

}

nextflow_setup() {

    load_modules nextflow/25.04.6
    export NXF_SINGULARITY_CACHEDIR="$MYSCRATCH/.singularity"
    export NXF_HOME="$MYSCRATCH/.nextflow"
    mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$NXF_HOME"

}
