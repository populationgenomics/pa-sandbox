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

import dataclasses
import logging
from collections.abc import Callable

from cpg_flow.filetypes import CramPath
from cpg_flow.stage import (
    CohortStage,
    DatasetStage,
    SequencingGroupStage,
    StageInput,
    StageInputNotFoundError,
    StageOutput,
    stage,
)
from cpg_flow.targets import Cohort, Dataset, SequencingGroup
from cpg_flow.utils import exists
from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from sandbox.jobs import somalier, verifybamid
from sandbox.jobs.multiqc import multiqc
from sandbox.jobs.picard import vcf_qc


# There are the above QC functionality, but I wonder how many of these metrics
# that they calculate and output are actually already provided by DRAGEN?
# If DRAGEN already provides them, then we can skip running these tools again
# and just use the DRAGEN outputs directly in MultiQC.
# TODO: DRAGEN provides all of these metrics in its output except VerifyBamID and Somalier.
@dataclasses.dataclass
class QcOut:
    """QC output file"""

    suf: str
    multiqc_key: str


@dataclasses.dataclass
class Qc:
    """QC function definition, and corresponding outputs"""

    func: Callable | None
    outs: dict[str, QcOut | None]
    optional: bool = False


def qc_functions() -> list[Qc]:
    """
    QC functions and their outputs for MultiQC aggregation
    """
    if config_retrieve(['workflow', 'skip_qc'], False):
        return []

    return [
        Qc(func=somalier.extract, outs={'somalier': None}),
        Qc(
            func=verifybamid,
            outs={'verify_bamid': QcOut('.verify-bamid.selfSM', 'verifybamid/selfsm')},
        ),
    ]

@stage()
class DragenCramQC(SequencingGroupStage):
    """
    Calling tools that process CRAM for QC purposes.
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        outs = {}
        for qc in qc_functions():
            for key, out in qc.outs.items():
                if key == 'somalier':
                    # Somalier outputs will be written to self.dataset.prefix() / 'cram' / f'{self.id}.cram.somalier' regardless of input cram path.
                    outs[key] = sequencing_group.make_cram_path().somalier_path
                elif out:
                    outs[key] = sequencing_group.dataset.prefix() / 'dragen_qc' / key / f'{sequencing_group.id}{out.suf}'
        return outs

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        cram_path = inputs.as_path(sequencing_group, Align, 'cram')
        crai_path = inputs.as_path(sequencing_group, Align, 'crai')

        jobs = []
        # This should run if either the stage or the sequencing group is being forced.
        forced = self.forced or sequencing_group.forced
        for qc in qc_functions():
            out_path_kwargs = {
                f'out_{key}_path': self.expected_outputs(sequencing_group)[key] for key in qc.outs.keys()
            }
            if qc.func:
                j = qc.func(  # type: ignore
                    get_batch(),
                    CramPath(cram_path, crai_path),
                    job_attrs=self.get_job_attrs(sequencing_group),
                    overwrite=forced,
                    **out_path_kwargs,
                )
                if j:
                    jobs.append(j)

        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=jobs)

@stage(required_stages=[DragenCramQC])
class SomalierPedigree(DatasetStage):
    """
    Checks pedigree from CRAM fingerprints.
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        """
        Return the report for MultiQC, plus putting an HTML into the web bucket.
        MultiQC expects the following patterns:
        * *.samples.tsv
        * *.pairs.tsv
        https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns
        .yaml#L472-L481
        """
        if config_retrieve(['workflow', 'skip_qc'], False):
            return {}

        prefix = dataset.prefix() / 'somalier' / 'cram' / dataset.get_alignment_inputs_hash()
        web_prefix = dataset.web_prefix() / 'somalier' / 'cram' / dataset.get_alignment_inputs_hash()
        return {
            'samples': prefix / f'{dataset.name}.samples.tsv',
            'expected_ped': prefix / f'{dataset.name}.expected.ped',
            'pairs': prefix / f'{dataset.name}.pairs.tsv',
            'html': web_prefix / 'cram-somalier-pedigree.html',
            'checks': prefix / f'{dataset.name}-checks.done',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput | None:
        """
        Checks calls job from the pedigree module
        """
        verifybamid_by_sgid = {}
        somalier_path_by_sgid = {}
        for sequencing_group in dataset.get_sequencing_groups():
            if config_retrieve(['somalier', 'exclude_high_contamination'], False):
                verify_bamid_path = inputs.as_path(stage=DragenCramQC, target=sequencing_group, key='verify_bamid')
                if not exists(verify_bamid_path):
                    logging.warning(
                        f'VerifyBAMID results {verify_bamid_path} do not exist for '
                        f'{sequencing_group}, somalier pedigree estimations might be affected',
                    )
                else:
                    verifybamid_by_sgid[sequencing_group.id] = verify_bamid_path
            somalier_path = inputs.as_path(stage=DragenCramQC, target=sequencing_group, key='somalier')
            somalier_path_by_sgid[sequencing_group.id] = somalier_path

        html_path = self.expected_outputs(dataset)['html']
        if base_url := dataset.web_url():
            html_url = str(html_path).replace(str(dataset.web_prefix()), base_url)
        else:
            html_url = None

        if any(sg.pedigree.dad or sg.pedigree.mom for sg in dataset.get_sequencing_groups()):
            expected_ped_path = dataset.write_ped_file(self.expected_outputs(dataset)['expected_ped'])
            jobs = somalier.pedigree(
                b=get_batch(),
                dataset=dataset,
                expected_ped_path=expected_ped_path,
                somalier_path_by_sgid=somalier_path_by_sgid,
                verifybamid_by_sgid=verifybamid_by_sgid,
                out_samples_path=self.expected_outputs(dataset)['samples'],
                out_pairs_path=self.expected_outputs(dataset)['pairs'],
                out_html_path=html_path,
                out_html_url=html_url,
                out_checks_path=self.expected_outputs(dataset)['checks'],
                job_attrs=self.get_job_attrs(dataset),
                send_to_slack=config_retrieve(['workflow', 'somalier_pedigree', 'send_to_slack'], default=True),
            )
            return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=jobs)
        return self.make_outputs(dataset, skipped=True)



# @stage(required_stages=[DragenCramQC, SomalierPedigree], analysis_type='qc', analysis_keys=['json'])
@stage(analysis_type='qc', analysis_keys=['json'])
class DragenCramMultiQC(CohortStage):
    """
    Run MultiQC to aggregate CRAM QC stats across a Cohort, rather than a Dataset.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Expected to produce an HTML and a corresponding JSON file.
        """
        if config_retrieve(['workflow', 'skip_qc'], False):
            return {}

        # get the unique hash for these Sequencing Groups
        sg_hash = cohort.get_alignment_inputs_hash()
        return {
            'html': cohort.dataset.web_prefix() / 'qc' / 'cram' / sg_hash / 'cohort_multiqc.html',
            'json': cohort.dataset.prefix() / 'qc' / 'cram' / sg_hash / 'cohort_multiqc_data.json',
            'checks': cohort.dataset.prefix() / 'qc' / 'cram' / sg_hash / '.cohort_checks',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Call a function from the `jobs` module using inputs from `cramqc`
        and `somalier` stages.
        """
        if config_retrieve(['workflow', 'skip_qc'], False):
            return self.make_outputs(cohort)

        json_path = self.expected_outputs(cohort)['json']
        html_path = self.expected_outputs(cohort)['html']
        checks_path = self.expected_outputs(cohort)['checks']
        if base_url := cohort.dataset.web_url():
            html_url = str(html_path).replace(str(cohort.dataset.web_prefix()), base_url)
        else:
            html_url = None

        cohort_sgs = cohort.get_sequencing_groups()

        dragen_metrics_paths: list[Path] = [
            to_path(f'gs://cpg-bioheart-test/ica/dragen_3_7_8/output/dragen_metrics/{sg.id}')
            for sg in cohort_sgs
            ]
        paths = []

        try:
            somalier_samples = inputs.as_path(cohort, SomalierPedigree, key='samples')
            somalier_pairs = inputs.as_path(cohort, SomalierPedigree, key='pairs')
        except StageInputNotFoundError:
            pass
        else:
            paths = [
                somalier_samples,
                somalier_pairs,
            ]

        ending_to_trim = set()  # endings to trim to get sample names
        modules_to_trim_endings = set()

        # NOTE: Commenting out because we need a way to represent default DRAGEN CRAM paths.
        # for sequencing_group in cohort.get_sequencing_groups():
        #     for qc in qc_functions():
        #         for key, out in qc.outs.items():
        #             if not out:
        #                 continue
        #             try:
        #                 path = inputs.as_path(sequencing_group, CramQC, key)
        #             except StageInputNotFoundError:  # allow missing inputs
        #                 logging.warning(
        #                     f'Output CramQc/"{key}" not found for {sequencing_group}, '
        #                     f'it will be silently excluded from MultiQC',
        #                 )
        #                 continue
        #             modules_to_trim_endings.add(out.multiqc_key)
        #             paths.append(path)
        #             ending_to_trim.add(path.name.replace(sequencing_group.id, ''))

        if not paths:
            logging.warning('No CRAM QC found to aggregate with MultiQC')
            return self.make_outputs(cohort)

        send_to_slack = config_retrieve(['workflow', 'cram_multiqc', 'send_to_slack'], default=True)
        extra_config = config_retrieve(['workflow', 'cram_multiqc', 'extra_config'], default={})
        extra_config['table_columns_visible'] = {'FastQC': False}

        paths += dragen_metrics_paths
        jobs = multiqc(
            get_batch(),
            tmp_prefix=cohort.dataset.tmp_prefix() / 'multiqc' / 'cram',
            paths=paths,
            ending_to_trim=ending_to_trim,
            modules_to_trim_endings=modules_to_trim_endings,
            dataset=cohort.dataset,
            out_json_path=json_path,
            out_html_path=html_path,
            out_html_url=html_url,
            out_checks_path=checks_path,
            job_attrs=self.get_job_attrs(cohort),
            sequencing_group_id_map=cohort.dataset.rich_id_map(),
            label='CRAM',
            send_to_slack=send_to_slack,
            extra_config=extra_config,
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)

@stage()
class DragenGvcfQC(SequencingGroupStage):
    """
    Calling tools that process GVCF for QC purposes.
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        """
        Generate a GVCF and corresponding TBI index, as well as QC.
        """
        outs: dict[str, Path] = {}
        if not config_retrieve(['workflow', 'skip_qc'], False):
            qc_prefix = sequencing_group.dataset.prefix() / 'qc' / sequencing_group.id
            outs |= {
                'qc_summary': to_path(f'{qc_prefix}.variant_calling_summary_metrics'),
                'qc_detail': to_path(f'{qc_prefix}.variant_calling_detail_metrics'),
            }
        return outs

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        """
        Use function from the jobs module
        """
        gvcf_path = inputs.as_path(sequencing_group, Genotype, 'gvcf')

        j = vcf_qc(
            b=get_batch(),
            vcf_or_gvcf=GvcfPath(gvcf_path).resource_group(get_batch()),
            is_gvcf=True,
            job_attrs=self.get_job_attrs(sequencing_group),
            output_summary_path=self.expected_outputs(sequencing_group)['qc_summary'],
            output_detail_path=self.expected_outputs(sequencing_group)['qc_detail'],
            overwrite=sequencing_group.forced,
        )
        return self.make_outputs(sequencing_group, data=self.expected_outputs(sequencing_group), jobs=[j])

@stage(
    required_stages=[DragenGvcfQC],
    analysis_type='qc',
    analysis_keys=['json'],
    update_analysis_meta=_update_meta,
)
class GvcfMultiQC(CohortStage):
    """
    Run MultiQC to summarise all GVCF QC.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Expected to produce an HTML and a corresponding JSON file.
        """
        if config_retrieve(['workflow', 'skip_qc'], False):
            return {}

        # get the unique hash for these Sequencing Groups
        sg_hash = cohort.get_alignment_inputs_hash()
        return {
            'html': cohort.web_prefix() / 'qc' / 'gvcf' / sg_hash / 'multiqc.html',
            'json': cohort.prefix() / 'qc' / 'gvcf' / sg_hash / 'multiqc_data.json',
            'checks': cohort.prefix() / 'qc' / 'gvcf' / sg_hash / '.checks',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:
        """
        Collect QC.
        """
        if config_retrieve(['workflow', 'skip_qc'], False):
            return self.make_outputs(cohort)

        json_path = self.expected_outputs(cohort)['json']
        html_path = self.expected_outputs(cohort)['html']
        checks_path = self.expected_outputs(cohort)['checks']
        if base_url := cohort.web_url():
            html_url = str(html_path).replace(str(cohort.web_prefix()), base_url)
        else:
            html_url = None

        paths = []
        ending_to_trim = set()  # endings to trim to get sample names

        for sequencing_group in cohort.get_sequencing_groups():
            for _stage, key in [(DragenGvcfQC, 'qc_detail')]:
                try:
                    path = inputs.as_path(sequencing_group, _stage, key)
                except StageInputNotFoundError:  # allow missing inputs
                        logging.warning(
                            f'Output {_stage.__name__}/"{key}" not found for {sequencing_group}, '
                            f'it will be silently excluded from MultiQC',
                        )
                else:
                    paths.append(path)
                    ending_to_trim.add(path.name.replace(sequencing_group.id, ''))

        if not paths:
            logging.warning('No GVCF QC found to aggregate with MultiQC')
            return self.make_outputs(cohort)

        modules_to_trim_endings = {'picard/variant_calling_metrics'}

        send_to_slack = config_retrieve(['workflow', 'gvcf_multiqc', 'send_to_slack'], default=True)
        extra_config = config_retrieve(['workflow', 'gvcf_multiqc', 'extra_config'], default={})
        extra_config['table_columns_visible'] = {'Picard': True}

        jobs = multiqc(
            get_batch(),
            tmp_prefix=cohort.tmp_prefix() / 'multiqc' / 'gvcf',
            paths=paths,
            ending_to_trim=ending_to_trim,
            modules_to_trim_endings=modules_to_trim_endings,
            cohort=cohort,
            out_json_path=json_path,
            out_html_path=html_path,
            out_html_url=html_url,
            out_checks_path=checks_path,
            job_attrs=self.get_job_attrs(cohort),
            sequencing_group_id_map=cohort.rich_id_map(),
            extra_config=extra_config,
            send_to_slack=send_to_slack,
            label='GVCF',
        )
        return self.make_outputs(cohort, data=self.expected_outputs(cohort), jobs=jobs)
@stage()
class CheckSampleMetrics(CohortStage):
    pass
