#!/bin/bash -l
set -euo pipefail

#################################################
############### Edit these values ###############
#################################################

PIPELINE_NAME=MCWorkflow
TMUX_SESSION_NAME=MCWorkflow
REPO_PATH=/path/to/repository/MCWorkflow
PIXI_COMMAND="pixi run nextflow main.nf -profile uppmax -params-file params.yml --project allocation"
SINGULARITY_CACHE_DIRECTORY="/path/to/shared/singularity/cache"




#################################################
############### No edit below     ###############
#################################################

# Set cluster variable
function get_cluster_name {
    if command -v sacctmgr >/dev/null 2>&1; then
        sacctmgr show cluster -P -n | \
        cut -f1 -d'|' | \
        head -n 1
    else
        echo "unknown"
    fi
}
CLUSTER=$( get_cluster_name )

# Welcome message
function print_welcome_message() {
    local tmux_version="$1"

    cat <<EOF
Welcome to $PIPELINE_NAME

Running on system: $CLUSTER.

Your job is running in the background via tmux version $tmux_version (session name: $TMUX_SESSION_NAME).

The command supplied was:
$PIXI_COMMAND

You can monitor workflow progress either via the '.nextflow.log' file, or with:
  squeue -u $USER

View your active tmux sessions with:
  tmux list-sessions

Cancel this job with:
  tmux kill-session -t $TMUX_SESSION_NAME
EOF
}

# Launch function
function launch_workflow {

    if [ ! -d "$REPO_PATH" ]; then
        echo "Error: '$REPO_PATH' does not exist."
        return 1
    else
        cd "$REPO_PATH"
    fi

    if ! TMUX_VERSION_OUTPUT=$(tmux -V 2>&1); then
        echo "Error: tmux not found"
        return 1
    fi

    if tmux has-session -t "$TMUX_SESSION_NAME" 2>/dev/null; then
        echo "Error: tmux session '$TMUX_SESSION_NAME' already exists."
        echo "Hint: run 'tmux kill-session -t $TMUX_SESSION_NAME' to shut it down."
        return 1
    fi

    # Set up environment variables
    if [ ! -d "$SINGULARITY_CACHE_DIRECTORY" ]; then
        echo "Error: '$SINGULARITY_CACHE_DIRECTORY' does not exist."
        return 1
    else
        export NXF_SINGULARITY_CACHEDIR="$SINGULARITY_CACHE_DIRECTORY"
    fi

    MAJORVERSION=$(echo "$TMUX_VERSION_OUTPUT" | awk '{split($2,v,"."); print v[1]}')

    # Launch run
    if [ "$MAJORVERSION" -eq 1 ]; then
        tmux new-session -s "$TMUX_SESSION_NAME" -d
        tmux send-keys -t "$TMUX_SESSION_NAME" "$PIXI_COMMAND" C-m
    else
        tmux new -s "$TMUX_SESSION_NAME" -d /bin/bash -c "$PIXI_COMMAND"
    fi

    print_welcome_message "$MAJORVERSION"

}

# Cluster specific preamble
if [ "$CLUSTER" == "dardel" ]; then
    module load PDC apptainer
    export APPTAINER_CACHEDIR=$PDC_TMP/apptainer/cache
    export SINGULARITY_CACHEDIR=$PDC_TMP/singularity/cache
    launch_workflow
elif [ "$CLUSTER" == "nac" ]; then
    module load Singularity
    launch_workflow
else
    echo "WARN: this launch script contains no configuration for '$CLUSTER'." >&2
    echo "Hint: you may need to pre-load the requested container platform (e.g. Apptainer/Singularity)" >&2
    launch_workflow
fi
