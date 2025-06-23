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

from cpg_flow.stage import SequencingGroupStage, stage

from deep_variant.jobs.run_deepvariant import run

if TYPE_CHECKING:
    # Path is a classic return type for a Stage, and is a shortcut for [CloudPath | pathlib.Path]
    from cpg_flow.stage import StageInput, StageOutput
    from cpg_flow.targets import SequencingGroup
    from cpg_utils import Path


@stage()
class RunDeepVariant(SequencingGroupStage):
    def expected_outputs(self, sequencing_group: 'SequencingGroup') -> 'Path':
        # self.prefix() is a more concise shortcut for multicohort.analysis_dataset_bucket/ StageName / Hash
        return sequencing_group.dataset.prefix(category='tmp') / self.name / f'{sequencing_group.id}.g.vcf.gz'

    def queue_jobs(self, sequencing_group: 'SequencingGroup', inputs: 'StageInput') -> 'StageOutput':  # noqa: ARG002
        """
        This is where we generate jobs for this stage.
        """
        # locate the intended output path
        outputs = self.expected_outputs(sequencing_group)

        # generate the output
        j = run(str(outputs))

        # return the jobs and outputs
        return self.make_outputs(sequencing_group, data=outputs, jobs=j)


# @stage(required_stages=[RunDeepVariant])
# class PrintPreviousJobOutputInAPythonJob(MultiCohortStage):
#     """
#     This is a stage that cats the output of a previous stage to the logs.
#     This uses a method imported from a job file, run as a PythonJob.
#     """

#     def expected_outputs(self, multicohort: 'MultiCohort') -> 'Path':
#         return multicohort.analysis_dataset.prefix(category='tmp') / self.name / 'cat.txt'

#     def queue_jobs(self, multicohort: 'MultiCohort', inputs: 'StageInput') -> 'StageOutput':
#         # get the previous stage's output
#         previous_stage = inputs.as_str(multicohort, DoSomethingGenericWithBash)

#         # localise the file
#         local_input = get_batch().read_input(previous_stage)

#         # generate the expected output path
#         outputs = self.expected_outputs(multicohort)

#         # run the PythonJob
#         job = get_batch().new_python_job(f'Read {previous_stage}')
#         job.image(config_retrieve(['workflow', 'driver_image']))
#         pyjob_output = job.call(
#             print_file_contents,
#             local_input,
#         )
#         get_batch().write_output(pyjob_output.as_str(), str(outputs))

#         return self.make_outputs(multicohort, data=outputs, jobs=job)
