LocalPath='/home/zzhu/tmp'

usage_str="Usage: $0 [ -i ID ] [ -t DNAnexus access token ] [ -p Project ID ] [-lp Local Path] [-o Output Dir] [-h]"
usage() {                                 # Function: Print a help message.
  echo "${usage_str}"
}

exit_abnormal() {                         # Function: Exit with error.
  usage
  exit 1
}

Help()
{
  # Display Help
   echo "${usage_str}"
   echo -e "\nOptions:"
   echo "-t     Specify a DNAnexus access token for easy login."
   echo "-i     Provide a Vidium Accession ID for sample."
   echo "-p     Project ID for the sample to be downloaded."
   echo "-lp    Path of local folder for the files to be downloaded."
   echo -e "-h     Print this Help.\n"
}

if [[ $# -eq 0 ]] ; then
    Help
    exit 0
fi

while getopts ":i:ht:p:lp:" options; do         # Loop: Get the next option;
                                          # use silent error checking;
                                          # options n and t take arguments.
  case "${options}" in                    # 
    i)                                    # If the option is n,
      ID=${OPTARG}                      # set $NAME to specified value.  VDM-21001170
      ;;
    t)                                    # If the option is t,
      dx_apitoken=${OPTARG} 
      ;;
    p)
      projectID=${OPTARG}
      exit;;
    lp)
      LocalPath=${OPTARG}
      exit;;
    h)
      Help
      exit;;                    
    :)                                    # If expected argument omitted:
      echo "Error: -${OPTARG} requires an argument."
      exit_abnormal                       # Exit abnormally.
      ;;
    *)                                    # If unknown (any other) option:
      exit_abnormal                       # Exit abnormally.
      ;;
  esac
done

#
# Downloading data from DNAnexus
# only download when there are only one file that matches the given name
#
downloadFileFromDNAnexus() { 
  local errCode=0
  rFile=$1
  lFile=$2
  if [[ -e "$lFile}" ]]
  then
    echo "Warning - ${lFile} exists"
    return 1;
  else
    # list file at DNANexus
    nf=`dx ls -l ${rFile} > ${LocalPath}/.tmp.txt`
    if [[ "$?" -ne 0 ]]
    then
      echo "Error in accessing DNANexus file - ${rFile}"
      echo "      no such file"
      return 1
    else
      nf=`wc -l ${LocalPath}/.tmp.txt`
      if [[ "$nf" -ne 1 ]]
      then
        # more than one file than matches the given name
        cat ${LocalPath}/tmp.txt
        echo "Error: more than one file for ${rFile}"
        return 2
      else
        # we find one and only one more. So download ....
        cmd=`dx download -o  ${lFile} ${rFile}`
        if [[ "$?" -eq 0 ]]
        then
          echo "${lFile}" >> "${LocalPath}/downloaded.txt"
          return 0
        else
          echo "Download failed!"
          return 3
        fi
      fi
    fi
  fi
}

projectID='projectX'

dx login --token $dx_apitoken --noprojects
dx select $projectID

remote_json="/Output/reports/vcfanno.json"
local_json="${LocalPath}/${ID}.json"
echo -e "\n##################################"
echo "# Downloading data from DNAnexus "
echo -e "##################################\n"

echo "#### (1) Downloading JSON ####"
ret_code=`downloadFileFromDNAnexus ${remote_json}  ${local_json}`
if [[ "$?" -ne 0 ]]
then 
  echo "JSON download failed!"
  echo "$ret_code"
  exit
fi
