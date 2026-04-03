export STORAGE1=/storage1/fs1/dinglab/Active/
export LSF_DOCKER_VOLUMES="$HOME:$HOME $STORAGE1:$STORAGE1"
export PASSWORD=igvcp1
export LSF_DOCKER_PRESERVE_ENVIRONMENT=true

# Take note of the node your job lands on (compute1-exec-263 for example)
PATH="/app:/IGV:$PATH" LSF_DOCKER_PORTS='8080:8080' PASSWORD='igvcp1' bsub -G compute-dinglab -Is -R 'select[port8080=1]' -q dinglab-interactive -a 'docker(fernandamrod/novnc_igv_cp1)' supervisord -c /app/supervisord.conf


