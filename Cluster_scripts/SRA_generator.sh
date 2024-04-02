SRAIDS=("SRR26383972" "SRR26383973" "SRR26383974" "SRR26383975" "SRR26383982" "SRR26383983" "SRR26383984" "SRR26383985")
mkdir "./slurmScripts/"
# loop through each SRA ID and create a new slurm file for each
for SRAID in ${SRAIDS[@]}
do
    sed "s/%%SRAID%%/$SRAID/g" SRA_template.sh > "./slurmScripts/SRA_$SRAID.sh"
done
