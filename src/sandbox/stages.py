"""
This file exists to define all the Stages for the workflow.
The logic for each stage can be contained here (if it is not too complex),
or can be delegated to a separate file in jobs.

Naming conventions for Stages are not enforced, but a series of recommendations have been made here:

https://cpg-populationanalysis.atlassian.net/wiki/spaces/ST/pages/185597962/Pipeline+Naming+Convention+Specification

A suggested naming convention for a stages is:
  - PascalCase (each word capitalized, no hyphens or underscores)
  - If the phrase contains an initialism (e.g. VCF), only the first character should be capitalised
  - Verb + Subject (noun) + Preposition + Direct Object (noun)  TODO(anyone): please correct my grammar is this is false
  e.g. AlignShortReadsWithBowtie2, or MakeSitesOnlyVcfWithBcftools
  - This becomes self-explanatory when reading the code and output folders

Each Stage should be a Class, and should inherit from one of
  - SequencingGroupStage
  - DatasetStage
  - CohortStage
  - MultiCohortStage
"""

from typing import TYPE_CHECKING

import cpg_utils
from cpg_flow.stage import CohortStage, StageInput, StageOutput, stage
from cpg_flow.targets import Cohort
from cpg_utils.config import config_retrieve
from loguru import logger

from sandbox.jobs.generate_sites_table import generate_sites_table
from sandbox.jobs.vds_to_vcf import vds_to_vcf

if TYPE_CHECKING:
    # Path is a classic return type for a Stage, and is a shortcut for [CloudPath | pathlib.Path]
    from hailtop.batch.job import PythonJob


@stage()
class GenerateSitesTable(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> cpg_utils.Path:  # noqa: ARG002
        outpath: cpg_utils.Path = cpg_utils.to_path(config_retrieve(['generate_sites_table', 'sites_table_outpath']))
        logger.info(f'Outpath is : {outpath}')
        return outpath

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:  # noqa: ARG002
        outputs: cpg_utils.Path = self.expected_outputs(cohort=cohort)
        sites_table_outpath: str = config_retrieve(['generate_sites_table', 'sites_table_outpath'])
        j: PythonJob = generate_sites_table(cohort=cohort, sites_table_outpath=sites_table_outpath)

        return self.make_outputs(target=cohort, data=outputs, jobs=j)


@stage()
class ExportVdsToVcf(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> cpg_utils.Path:  # noqa: ARG002
        outdir: str = config_retrieve(['export_vds_to_vcf', 'vcf_out_dir'])
        logger.info(f'Outdir is : {outdir}')
        chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
        return {f'{chrom}': cpg_utils.to_path(outdir) / f'{chrom}.vcf.bgz' for chrom in chroms}

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:  # noqa: ARG002
        vcf_out_dict: dict[str, cpg_utils.Path] = self.expected_outputs(cohort=cohort)
        vds_path: str = config_retrieve(['export_vds_to_vcf', 'vds_path'])

        jobs = []
        for chrom, outpath in vcf_out_dict.items():
            j: PythonJob = vds_to_vcf(cohort=cohort, vds_path=vds_path, vcf_outpath=str(outpath), chrom=chrom)
            jobs.append(j)

        return self.make_outputs(target=cohort, data=vcf_out_dict, jobs=jobs)
