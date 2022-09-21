#!/bin/bash


# You may source this bashrc-extension, if you wish to have some shortcuts which can be nice for development and deployment.



# Neat aliases
alias mconda='mamba'
alias survey='watch "sensors; nvidia-smi"'
alias citament='git add -u && git commit -m "amend" &&  git pull && git push && echo OK'



#Human Monkeypox Virus 2022
start_viraligner_mpx() {

	example="\n\texample:\n\t start_viraligner_mpx path/to/my_samplesheet.xlsx path/to/my_rundir\n"

if [ -z "$1" ]; then
    echo "Input error: The samplesheet argument is empty. Please specify a samplesheet."
    echo -e $example
     
elif [ -z "$2" ]; then
    echo "Input error: The rundir argument is empty. Please specify a rundir."
    echo -e $example
else

    #clear
    cd ~/viraligner && touch CoverMon.flag && conda activate viraligner && snakemake --profile default --use-conda --conda-frontend mamba --config scheme_directory="primer_schemes" scheme_name="hMPX-2022" scheme_version="1" samplesheet="${1}" rundir="${2}" ${3} && echo && cowsay "viraligner_mpx finished successfully, you may now close this window."
fi
}



