# Manually input Index directory that contains human-2020-A_splici3, and which SRA files to use
SRAIDS=("SRR26383972" "SRR26383973" "SRR26383974" "SRR26383975" "SRR26383982" "SRR26383983" "SRR26383984" "SRR26383985")
# Must follow $AF_SAMPLE_DIR/Index_X_X/human-2020-A_splici
INDEX="Index150_150"
# loop through each SRA ID and create a new slurm file for each
for SRAID in ${SRAIDS[@]}
do
    sed "s/%%SRAID%%/$SRAID/g" SAF_template.slurm > "./slurmScripts/SAF_$SRAID.slurm"
    sed -i "s/%%INDEX%%/$INDEX/g" "./slurmScripts/SAF_$SRAID.slurm"
done
