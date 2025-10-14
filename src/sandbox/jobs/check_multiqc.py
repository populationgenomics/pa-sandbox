#!/usr/bin/env python3

"""
Checks metrics in MultiQC output, based on thresholds in the qc_thresholds
config section.

Script can send a report to a Slack channel. To enable that, set SLACK_TOKEN
and SLACK_CHANNEL environment variables, and add "Seqr Loader" app into
a channel with:

/invite @Seqr Loader
"""
import json
import logging
import pprint
from collections import defaultdict

import click
from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.slack import send_message

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)


@click.command()
@click.option(
    '--multiqc-json',
    'multiqc_json_path',
    required=True,
    help='Path to MultiQC JSON output',
)
@click.option(
    '--html-url',
    'html_url',
    help='MultiQC HTML URL',
)
@click.option('--dataset', 'dataset', help='Dataset name')
@click.option('--title', 'title', help='Report title')
@click.option(
    '--send-to-slack/--no-send-to-slack',
    'send_to_slack',
    help='Send log to Slack message, according to environment variables SLACK_CHANNEL and SLACK_TOKEN',
)
@click.option(
    '--failed-samples-path',
    'failed_samples_path',
    help='Path to write JSON file with failed samples and their failed metrics',
)
def main(
    multiqc_json_path: str,
    html_url: str | None = None,
    dataset: str | None = None,
    title: str | None = None,
    send_to_slack: bool = True,
    failed_samples_path: str | None = None,
):
    """
    Check metrics in MultiQC json and send info about failed samples
    as a Slack message.
    """
    run(
        multiqc_json_path=multiqc_json_path,
        html_url=html_url,
        dataset=dataset,
        title=title,
        send_to_slack=send_to_slack,
        failed_samples_path=failed_samples_path,
    )

QC_MAPPING = {
    'mean_coverage': {
        'multiqc_report_name': 'Average sequenced coverage over genome',
        'display_name': 'Mean Coverage',
    },
    'pct_genome_gt_20x': {
        'multiqc_report_name': 'wgs pct of genome with coverage [20x:inf)',
        'display_name': 'Pct Genome @ >20x',
    },
    'q30_bases': {
        'multiqc_report_name': 'Q30 bases',
        'display_name': 'Q30 Bases',
    },
    'contamination_verifybamid': {
        'multiqc_report_name': 'FREEMIX',
        'display_name': 'Contamination (VerifyBamID)',
    },
    'contamination_dragen': {
        'multiqc_report_name': 'Estimated sample contamination',
        'display_name': 'Contamination (DRAGEN)',
    },
    'mapping_rate_pct': {
        'multiqc_report_name': 'Mapped reads pct',
        'display_name': 'Mapping Rate (%)',
    },
    'duplication_rate_pct': {
        'multiqc_report_name': 'Number of duplicate marked reads pct',
        'display_name': 'Duplication Rate (%)',
    },
    'chimera_rate': {
        'calculator': lambda d: d.get('Supplementary (chimeric) alignments', 0) / d.get('Total alignments', 1),
        'display_name': 'Chimera Rate',
    },
    'mean_insert_size': {
        'multiqc_report_name': 'Insert length: mean',
        'display_name': 'Mean Insert Size',
    },
    'insert_size_sd': {
        'multiqc_report_name': 'Insert length: standard deviation',
        'display_name': 'Insert Size SD',
    },
    'ti_tv_ratio': {
        'multiqc_report_name': 'Ti/Tv ratio',
        'display_name': 'Ti/Tv Ratio (SNPs)',
    },
    'het_hom_ratio': {
        'multiqc_report_name': 'Het/Hom ratio',
        'display_name': 'Het/Hom Ratio',
    },
}

def build_qc_thresholds(seq_type: str, config_key: str) -> dict[str, dict]:
    """
    Build a dictionary of desired QC thresholds from config.
    Example config structure:
        [qc_thresholds.genome.min]
        mean_coverage = 30
        q30_bases = 8e10
        [qc_thresholds.genome.max]
        contamination_verifybamid = 0.05
        contamination_dragen = 0.03
        chimera_rate = 0.03
    """
    threshold_d = get_config()['qc_thresholds'].get(seq_type, {}).get(config_key, {})
    qc_thresholds = {}
    for metric, threshold in threshold_d.items():
        if metric in QC_MAPPING:
            qc_thresholds[metric] = {
                'threshold': threshold,
                **QC_MAPPING[metric],
            }
        else:
            logging.warning(
                f"Metric '{metric}' has a threshold but is not defined in QC_MAPPING. "
                f"Using default names."
            )
            qc_thresholds[metric] = {
                'threshold': threshold,
                'multiqc_report_name': metric,
                'display_name': metric,
            }
    return qc_thresholds


def run(
    multiqc_json_path: str,
    html_url: str | None = None,
    dataset: str | None = None,
    title: str | None = None,
    send_to_slack: bool = True,
    failed_samples_path: str | None = None,
):
    seq_type = get_config()['workflow']['sequencing_type']

    with to_path(multiqc_json_path).open() as f:
        d = json.load(f)
        sections = d['report_general_stats_data']

    bad_lines_by_sample = defaultdict(list)
    for min_or_max, fail_sign, good_sign, is_fail in [
        ('min', '<', '≥', lambda val_, thresh_: val_ < thresh_),
        ('max', '>', '≤', lambda val_, thresh_: val_ > thresh_),
    ]:
        threshold_d = build_qc_thresholds(seq_type, min_or_max)
        logging.info(f'{min_or_max} thresholds: {pprint.pformat(threshold_d)}')
        for section_data in sections.values():
            for sg_id, val_by_metric in section_data.items():
                for metric_config in threshold_d.values():
                    val = None
                    # DRAGEN does not provide pct chimeras directly, so we calculate it
                    if 'calculator' in metric_config:
                        try:
                            val = metric_config['calculator'](val_by_metric)
                        except (KeyError, ZeroDivisionError):
                            continue
                    elif 'multiqc_report_name' in metric_config:
                        val = val_by_metric.get(metric_config['multiqc_report_name'])

                    if val is None:
                        continue

                    threshold = metric_config['threshold']
                    display_name = metric_config['display_name']

                    if is_fail(val, threshold):
                        line = f'{display_name}={val:.4f} {fail_sign} {threshold:.4f}'
                        bad_lines_by_sample[sg_id].append(line)
                        logging.warning(f'❗ {sg_id}: {line}')
                    else:
                        line = f'{display_name}={val:.4f} {good_sign} {threshold:.4f}'
                        logging.info(f'✅ {sg_id}: {line}')
    logging.info('')

    if bad_lines_by_sample and failed_samples_path:
        logging.info(f'Writing {len(bad_lines_by_sample)} failed sample(s) to {failed_samples_path}')
        with to_path(failed_samples_path).open('w') as f:
            json.dump(bad_lines_by_sample, f, indent=2)

    # Constructing Slack message
    if dataset and html_url:
        title = f'*[{dataset}]* <{html_url}|{title or "MultiQC report"}>'
    elif not title:
        title = 'MultiQC report'
    messages = []
    if bad_lines_by_sample:
        messages.append(f'{title}. {len(bad_lines_by_sample)} samples are flagged:')
        for sample, bad_lines in bad_lines_by_sample.items():
            messages.append(f'❗ {sample}: ' + ', '.join(bad_lines))
    else:
        messages.append(f'✅ {title}')
    text = '\n'.join(messages)
    logging.info(text)

    if send_to_slack:
        send_message(text)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
