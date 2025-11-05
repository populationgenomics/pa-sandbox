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
import json
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
from cpg_utils.config import config_retrieve, dataset_path
from cpg_utils.hail_batch import get_batch

from sandbox.jobs import somalier, verifybamid
from sandbox.jobs.multiqc import multiqc


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
        Qc(
            func=verifybamid.verifybamid,
            outs={'verify_bamid': QcOut('.verify-bamid.selfSM', 'verifybamid/selfsm')},
        ),
    ]

@stage()
class DragenCramQC(SequencingGroupStage):
    """
    Calling tools that process CRAM for QC purposes.
    """

    def expected_outputs(self, sequencing_group: SequencingGroup) -> dict[str, Path]:
        dragen_prefix = 'ica/dragen_3_7_8/qc'
        outs = {}
        for qc in qc_functions():
            for key, out in qc.outs.items():
                outs[key] = sequencing_group.dataset.prefix() / dragen_prefix / key / f'{sequencing_group.id}{out.suf}'
        return outs

    def queue_jobs(self, sequencing_group: SequencingGroup, inputs: StageInput) -> StageOutput | None:
        cram_path = sequencing_group.cram
        crai_path = sequencing_group.cram.index_path

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
            somalier_path = dataset_path(f'ica/dragen_3_7_8/output/somalier/{sequencing_group.id}.somalier')
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


@stage(required_stages=[DragenCramQC, SomalierPedigree], analysis_type='qc', analysis_keys=['json'])
class DragenCramMultiQC(CohortStage):
    """
    Run MultiQC to aggregate CRAM QC stats across a Cohort, rather than a Dataset.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        """
        Expected to produce an HTML and a corresponding JSON file.
        """
        dragen_prefix = 'ica/dragen_3_7_8/qc'
        return {
            'html': cohort.dataset.web_prefix() / dragen_prefix / cohort.id / 'multiqc' / 'cohort_multiqc.html',
            'json': cohort.dataset.prefix() / dragen_prefix / cohort.id / 'multiqc' / 'cohort_multiqc_data.json',
            'checks': cohort.dataset.prefix() / dragen_prefix / cohort.id / 'multiqc' / '.cohort_checks',
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
            to_path(dataset_path(f'ica/dragen_3_7_8/output/dragen_metrics/{sg.id}'))
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

        for sequencing_group in cohort.get_sequencing_groups():
            for qc in qc_functions():
                for key, out in qc.outs.items():
                    if not out:
                        continue
                    try:
                        path = inputs.as_path(sequencing_group, DragenCramQC, key)
                    except StageInputNotFoundError:  # allow missing inputs
                        logging.warning(
                            f'Output DragenCramQC/"{key}" not found for {sequencing_group}, '
                            f'it will be silently excluded from MultiQC',
                        )
                        continue
                    modules_to_trim_endings.add(out.multiqc_key)
                    paths.append(path)
                    ending_to_trim.add(path.name.replace(sequencing_group.id, ''))

        paths += dragen_metrics_paths

        if not paths:
            logging.warning('No CRAM QC found to aggregate with MultiQC')
            return self.make_outputs(cohort)

        send_to_slack = config_retrieve(['workflow', 'cram_multiqc', 'send_to_slack'], default=True)
        extra_config = config_retrieve(['workflow', 'cram_multiqc', 'extra_config'], default={})
        extra_config['table_columns_visible'] = {'FastQC': False}

        jobs = multiqc(
            get_batch(),
            tmp_prefix=cohort.dataset.tmp_prefix() / 'multiqc' / 'cram',
            paths=paths,
            ending_to_trim=ending_to_trim,
            modules_to_trim_endings=modules_to_trim_endings,
            cohort=cohort,
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
@stage(required_stages=[DragenCramMultiQC], analysis_type='qc', analysis_keys=['json'])
class RegisterDragenSampleFailures(CohortStage):
    """
    Check MultiQC report against defined thresholds.
    """

    def expected_outputs(self, cohort: Cohort) -> dict[str, Path]:
        return {
            'checks': cohort.dataset.prefix() / 'qc' / 'cram' / cohort.id / '.cohort_checks_registered',
        }

    def queue_jobs(self, cohort: Cohort, inputs: StageInput) -> StageOutput | None:

        cohort_sgs: list[SequencingGroup] = cohort.get_sequencing_groups()
        mqc_checks = json.loads(inputs.as_path(cohort, DragenCramMultiQC, key='checks'))

        return


