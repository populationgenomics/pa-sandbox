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
from cpg_flow.stage import CohortStage, SequencingGroupStage, StageInput, StageOutput, stage
from cpg_flow.targets import Cohort
from cpg_utils.config import config_retrieve

from sandbox.jobs.generate_sites_table import generate_sites_table
from sandbox.jobs.run_deepvariant import run

if TYPE_CHECKING:
    # Path is a classic return type for a Stage, and is a shortcut for [CloudPath | pathlib.Path]
    from cpg_flow.targets import SequencingGroup
    from cpg_utils import Path
    from hailtop.batch.job import PythonJob


@stage()
class RunPangenomeAwareDeepVariant(SequencingGroupStage):
    def expected_outputs(self, sequencing_group: 'SequencingGroup') -> 'Path':
        # self.prefix() is a more concise shortcut for multicohort.analysis_dataset_bucket/ StageName / Hash
        return {
            'gvcf': sequencing_group.dataset.prefix(category='tmp') / self.name / f'{sequencing_group.id}.g.vcf.gz',
            'vcf': sequencing_group.dataset.prefix(category='tmp') / self.name / f'{sequencing_group.id}.vcf.gz',
        }

    def queue_jobs(self, sequencing_group: 'SequencingGroup', inputs: 'StageInput') -> 'StageOutput':  # noqa: ARG002
        """
        This is where we generate jobs for this stage.
        """
        # locate the intended output path
        outputs = self.expected_outputs(sequencing_group)

        # generate the output
        j = run(
            output_vcf=str(outputs['vcf']),
            output_gvcf=str(outputs['gvcf']),
            sequencing_group_name=str(sequencing_group.id),
            cram_path=sequencing_group.cram or sequencing_group.make_cram_path(),
        )

        # return the jobs and outputs
        return self.make_outputs(sequencing_group, data=outputs, jobs=j)


@stage()
class GenerateSitesTable(CohortStage):
    def expected_outputs(self, cohort: Cohort) -> cpg_utils.Path:  # noqa: ARG002
        return cpg_utils.to_path(config_retrieve(['generate_sites_table', 'sites_table_outpath']))

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:  # noqa: ARG002
        outputs: cpg_utils.Path = self.expected_outputs(cohort=cohort)
        sites_table_outpath: str = config_retrieve(['generate_sites_table', 'sites_table_outpath'])
        j: PythonJob = generate_sites_table(cohort=cohort, sites_table_outpath=sites_table_outpath)

        return self.make_outputs(target=cohort, data=outputs, jobs=j)
