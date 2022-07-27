# # Create conda environment with Bismark installed

# See my lock file for env content

# Note: i used a patched version of Bismark, this is no longer necessary, just install with conda

# Note: I also used a patched version of FASTQ-screen. not sure if that would still be necessary. By default, FastQ-screen is anyway not run

# # Clone the smk_wgbs repo

# private repo -> send me your github name

# git clone https://github.com/stephenkraemer/smk_wgbs

# # Install python package

# pip -e /path/to/smk_wgbs

# # Copy and adapt the example yaml config file

# doc/smk_wgbs_default_config.yaml
#open screen
#open modul
#module load anaconda3/2019.07
#open conda 
#conda activate smk_wgbs
# # Obtain cluster submission scripts for your scheduler

# ie scripts: jobscript, submission, status

# # In python script

# ## Import snakemake, smk_wgbs

import snakemake
import smk_wgbs

# ## Setup a list of all samples you want to align

# Format either ['entity1', 'entity2', ...]
# Or [['entity1', 'sample1'], ...]
# Or a mix of both
# In ODCF terms, entity=Pid ("Patient1", 'HSCs'), sample=sample_type ("tumor", "normal", "cell01", ...)

plate_1_entities = [
    "JMMLC_D123_pr-sctm_p1_c0",
    "JMMLC_D123_pr-sctm_p1_c1",
    "JMMLC_D123_pr-sctm_p1_c25",
    "JMMLC_D124_pr-sctm_p1_c0",
    "JMMLC_D124_pr-sctm_p1_c1",
    "JMMLC_D124_pr-sctm_p1_c25",
    "JMMLC_D129_pr-sctm_p1_c0",
    "JMMLC_D129_pr-sctm_p1_c1",
    "JMMLC_D129_pr-sctm_p1_c25"
]

# ## Run workflow
snakemake.snakemake(
        snakefile=smk_wgbs.get_snakefile_fp(),
        latency_wait=60,
        configfiles=['/omics/groups/OE0219/internal/jmmlc_pbat/sc_methyl/210107_completeRun/config.yaml'],
        # Three config values are not usually specified in the YAML
        # Because they change for each run
        config=dict(
                # name for a given batch, eg one plate
                experiment_name="plate-1",
                # sctm or scnmt
                protocol="sctm",
                entities=plate_1_entities,
        ),
        nodes=5000,
        restart_times=4,
        keepgoing=True,
        jobscript="/omics/groups/OE0219/internal/scMeth_Workflow_Stephen/lsf-jobscript.sh",
        cluster="/omics/groups/OE0219/internal/scMeth_Workflow_Stephen/lsf-submit.py",
        cluster_status="/omics/groups/OE0219/internal/scMeth_Workflow_Stephen/lsf-status.py",
        max_jobs_per_second=10,
        force_incomplete=True,
        dryrun=False,
)
