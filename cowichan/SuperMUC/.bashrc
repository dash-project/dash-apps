module load cmake/3.6
module swap intel/16.0 intel/17.0.2
module swap mpi.ibm/1.4 mpi.intel/2017.2
module load papi
module load hwloc
export CC=icc
export CXX=icpc


# Print to STDERR.
echoerr()
{
  cat <<< "$@" 1>&2;
}


#------------------------------------------------
#-------   Randmat   ----------------------------
#------------------------------------------------

genrandmatcmdfile()
{
  #Problem InputParameters RunExecutionParameters
  if [ -z "$1" ]
    then echoerr "No command file template given"; return -1;
  fi
  if [ -z "$2" ]
    then echoerr "No matrix size given"; return -2;
  fi
  MATRIX_SIZE="$2"  
  if [ -z "$3" ]
    then echoerr "No run number given"; return -2;
  fi
  RUN_NR="$3"
  if [ -z "$4" ]
    then echoerr "No number of nodes given"; return -2;
  fi
  NUM_NODES="$4"
  
  if [ -z "$5" ]
  then
    NUM_THREADS=1
  else
    NUM_THREADS="$5"
  fi
  if [ -z "$6" ]
  then
    NUM_PROCS=28
  else
    NUM_PROCS="$6"
  fi
  UNITS=$(($NUM_PROCS*$NUM_NODES))
  COMMAND="echo $MATRIX_SIZE $MATRIX_SIZE 2 | mpiexec -n $UNITS ~/cowichan/randmat/randmat --is_bench"
  JOB_NAME_SUFFIX="mx$MATRIX_SIZE.n$NUM_NODES.p$NUM_PROCS.t$NUM_THREADS.r$RUN_NR"
  CMD_FILE_NAME="randmat.$JOB_NAME_SUFFIX"
  cat $1 | sed "s/<<NUM_PROCS>>/$NUM_PROCS/g" | \
           sed "s/<<NUM_THREADS>>/$NUM_THREADS/g" | \
           sed "s/<<NUM_NODES>>/$NUM_NODES/g" | \
           sed "s/<<JOBNAME>>/$CMD_FILE_NAME/g" | \
           sed "s%<<COMMAND>>%$COMMAND%g" \
         > ~/tmp/$CMD_FILE_NAME
  echo ~/tmp/$CMD_FILE_NAME
  return 0
}
export -f genrandmatcmdfile

benchrandmat()
{
  genrandmatcmdfile $1 $2 $3 $4 $5 $6 | xargs llsubmit
}
export -f benchrandmat


#------------------------------------------------
#-------   Thresh    ----------------------------
#------------------------------------------------

genthreshcmdfile()
{
  #Problem InputParameters RunExecutionParameters
  if [ -z "$1" ]
    then echoerr "No command file template given"; return -1;
  fi
  if [ -z "$2" ]
    then echoerr "No matrix size given"; return -2;
  fi
  MATRIX_SIZE="$2"  
	if [ -z "$3" ]
    then echoerr "No thresh value given"; return -2;
  fi
  THRESH="$3"
  if [ -z "$4" ]
    then echoerr "No run number given"; return -2;
  fi
  RUN_NR="$4"
  if [ -z "$5" ]
    then echoerr "No number of nodes given"; return -2;
  fi
  NUM_NODES="$5"
  
  if [ -z "$6" ]
  then
    NUM_THREADS=1
  else
    NUM_THREADS="$6"
  fi
  if [ -z "$7" ]
  then
    NUM_PROCS=28
  else
    NUM_PROCS="$7"
  fi
  UNITS=$(($NUM_PROCS*$NUM_NODES))
  COMMAND="echo $MATRIX_SIZE $MATRIX_SIZE $THRESH | mpiexec -n $UNITS ~/cowichan/thresh/thresh --is_bench"
  JOB_NAME_SUFFIX="mx$MATRIX_SIZE.trsh$THRESH.n$NUM_NODES.p$NUM_PROCS.t$NUM_THREADS.r$RUN_NR"
  CMD_FILE_NAME="thresh.$JOB_NAME_SUFFIX"
  cat $1 | sed "s/<<NUM_PROCS>>/$NUM_PROCS/g" | \
           sed "s/<<NUM_THREADS>>/$NUM_THREADS/g" | \
           sed "s/<<NUM_NODES>>/$NUM_NODES/g" | \
           sed "s/<<JOBNAME>>/$CMD_FILE_NAME/g" | \
           sed "s%<<COMMAND>>%$COMMAND%g" \
         > ~/tmp/$CMD_FILE_NAME
  echo ~/tmp/$CMD_FILE_NAME
  return 0
}
export -f genthreshcmdfile

benchthresh()
{
  genthreshcmdfile $1 $2 $3 $4 $5 $6 $7 | xargs llsubmit
}
export -f benchthresh

#------------------------------------------------
#-------   Winnow    ----------------------------
#------------------------------------------------

genwinnowcmdfile()
{
  #Problem InputParameters RunExecutionParameters
  if [ -z "$1" ]
    then echoerr "No command file template given"; return -1;
  fi
  if [ -z "$2" ]
    then echoerr "No matrix size given"; return -2;
  fi
  MATRIX_SIZE="$2"  
	if [ -z "$3" ]
    then echoerr "No thresh value given"; return -2;
  fi
  THRESH="$3"
	if [ -z "$4" ]
    then echoerr "No nelts value given"; return -2;
  fi
  NELTS="$4"
  if [ -z "$5" ]
    then echoerr "No run number given"; return -2;
  fi
  RUN_NR="$5"
  if [ -z "$6" ]
    then echoerr "No number of nodes given"; return -2;
  fi
  NUM_NODES="$6"
  
  if [ -z "$7" ]
  then
    NUM_THREADS=1
  else
    NUM_THREADS="$7"
  fi
  if [ -z "$8" ]
  then
    NUM_PROCS=28
  else
    NUM_PROCS="$8"
  fi
  UNITS=$(($NUM_PROCS*$NUM_NODES))
  COMMAND="echo $MATRIX_SIZE $MATRIX_SIZE $NELTS $THRESH | mpiexec -n $UNITS ~/cowichan/winnow/winnow --is_bench"
  JOB_NAME_SUFFIX="mx$MATRIX_SIZE.trsh$THRESH.nelts$NELTS.n$NUM_NODES.p$NUM_PROCS.t$NUM_THREADS.r$RUN_NR"
  CMD_FILE_NAME="winnow.$JOB_NAME_SUFFIX"
  cat $1 | sed "s/<<NUM_PROCS>>/$NUM_PROCS/g" | \
           sed "s/<<NUM_THREADS>>/$NUM_THREADS/g" | \
           sed "s/<<NUM_NODES>>/$NUM_NODES/g" | \
           sed "s/<<JOBNAME>>/$CMD_FILE_NAME/g" | \
           sed "s%<<COMMAND>>%$COMMAND%g" \
         > ~/tmp/$CMD_FILE_NAME
  echo ~/tmp/$CMD_FILE_NAME
  return 0
}
export -f genwinnowcmdfile

benchwinnow()
{
  genwinnowcmdfile $1 $2 $3 $4 $5 $6 $7 $8 | xargs llsubmit
}
export -f benchwinnow



#------------------------------------------------
#-------   Outer     ----------------------------
#------------------------------------------------

genoutercmdfile()
{
  #Problem InputParameters RunExecutionParameters
  if [ -z "$1" ]
    then echoerr "No command file template given"; return -1;
  fi
  if [ -z "$2" ]
    then echoerr "No number of elements given"; return -2;
  fi
  NELTS="$2"  
  if [ -z "$3" ]
    then echoerr "No run number given"; return -2;
  fi
  RUN_NR="$3"
  if [ -z "$4" ]
    then echoerr "No number of nodes given"; return -2;
  fi
  NUM_NODES="$4"
  
  if [ -z "$5" ]
  then
    NUM_THREADS=1
  else
    NUM_THREADS="$5"
  fi
  if [ -z "$6" ]
  then
    NUM_PROCS=28
  else
    NUM_PROCS="$6"
  fi
  UNITS=$(($NUM_PROCS*$NUM_NODES))
  COMMAND="echo $NELTS | mpiexec -n $UNITS ~/cowichan/outer/outer --is_bench"
  JOB_NAME_SUFFIX="nelts$NELTS.n$NUM_NODES.p$NUM_PROCS.t$NUM_THREADS.r$RUN_NR"
  CMD_FILE_NAME="outer.$JOB_NAME_SUFFIX"
  cat $1 | sed "s/<<NUM_PROCS>>/$NUM_PROCS/g" | \
           sed "s/<<NUM_THREADS>>/$NUM_THREADS/g" | \
           sed "s/<<NUM_NODES>>/$NUM_NODES/g" | \
           sed "s/<<JOBNAME>>/$CMD_FILE_NAME/g" | \
           sed "s%<<COMMAND>>%$COMMAND%g" \
         > ~/tmp/$CMD_FILE_NAME
  echo ~/tmp/$CMD_FILE_NAME
  return 0
}
export -f genoutercmdfile

benchouter()
{
  genoutercmdfile $1 $2 $3 $4 $5 $6 | xargs llsubmit
}
export -f benchouter


#------------------------------------------------
#-------   Product   ----------------------------
#------------------------------------------------

genproductcmdfile()
{
  #Problem InputParameters RunExecutionParameters
  if [ -z "$1" ]
    then echoerr "No command file template given"; return -1;
  fi
  if [ -z "$2" ]
    then echoerr "No number of elements given"; return -2;
  fi
  NELTS="$2"  
  if [ -z "$3" ]
    then echoerr "No run number given"; return -2;
  fi
  RUN_NR="$3"
  if [ -z "$4" ]
    then echoerr "No number of nodes given"; return -2;
  fi
  NUM_NODES="$4"
  
  if [ -z "$5" ]
  then
    NUM_THREADS=1
  else
    NUM_THREADS="$5"
  fi
  if [ -z "$6" ]
  then
    NUM_PROCS=28
  else
    NUM_PROCS="$6"
  fi
  UNITS=$(($NUM_PROCS*$NUM_NODES))
  COMMAND="echo $NELTS | mpiexec -n $UNITS ~/cowichan/product/product --is_bench"
  JOB_NAME_SUFFIX="nelts$NELTS.n$NUM_NODES.p$NUM_PROCS.t$NUM_THREADS.r$RUN_NR"
  CMD_FILE_NAME="product.$JOB_NAME_SUFFIX"
  cat $1 | sed "s/<<NUM_PROCS>>/$NUM_PROCS/g" | \
           sed "s/<<NUM_THREADS>>/$NUM_THREADS/g" | \
           sed "s/<<NUM_NODES>>/$NUM_NODES/g" | \
           sed "s/<<JOBNAME>>/$CMD_FILE_NAME/g" | \
           sed "s%<<COMMAND>>%$COMMAND%g" \
         > ~/tmp/$CMD_FILE_NAME
  echo ~/tmp/$CMD_FILE_NAME
  return 0
}
export -f genproductcmdfile

benchproduct()
{
  genproductcmdfile $1 $2 $3 $4 $5 $6 | xargs llsubmit
}
export -f benchproduct



#------------------------------------------------
#-------   Chain     ----------------------------
#------------------------------------------------

genchaincmdfile()
{
  #Problem InputParameters RunExecutionParameters
  if [ -z "$1" ]
    then echoerr "No command file template given"; return -1;
  fi
  if [ -z "$2" ]
    then echoerr "No matrix size given"; return -2;
  fi
  MATRIX_SIZE="$2"  
	if [ -z "$3" ]
    then echoerr "No thresh value given"; return -2;
  fi
  THRESH="$3"
	if [ -z "$4" ]
    then echoerr "No nelts value given"; return -2;
  fi
  NELTS="$4"
  if [ -z "$5" ]
    then echoerr "No run number given"; return -2;
  fi
  RUN_NR="$5"
  if [ -z "$6" ]
    then echoerr "No number of nodes given"; return -2;
  fi
  NUM_NODES="$6"
  
  if [ -z "$7" ]
  then
    NUM_THREADS=1
  else
    NUM_THREADS="$7"
  fi
  if [ -z "$8" ]
  then
    NUM_PROCS=28
  else
    NUM_PROCS="$8"
  fi
  UNITS=$(($NUM_PROCS*$NUM_NODES))
  COMMAND="echo $MATRIX_SIZE 2 $THRESH $NELTS | mpiexec -n $UNITS ~/cowichan/chain/chain --is_bench"
  JOB_NAME_SUFFIX="mx$MATRIX_SIZE.trsh$THRESH.nelts$NELTS.n$NUM_NODES.p$NUM_PROCS.t$NUM_THREADS.r$RUN_NR"
  CMD_FILE_NAME="chain.$JOB_NAME_SUFFIX"
  cat $1 | sed "s/<<NUM_PROCS>>/$NUM_PROCS/g" | \
           sed "s/<<NUM_THREADS>>/$NUM_THREADS/g" | \
           sed "s/<<NUM_NODES>>/$NUM_NODES/g" | \
           sed "s/<<JOBNAME>>/$CMD_FILE_NAME/g" | \
           sed "s%<<COMMAND>>%$COMMAND%g" \
         > ~/tmp/$CMD_FILE_NAME
  echo ~/tmp/$CMD_FILE_NAME
  return 0
}
export -f genchaincmdfile

benchchain()
{
  genchaincmdfile $1 $2 $3 $4 $5 $6 $7 $8 | xargs llsubmit
}
export -f benchchain
