"""Configuration of the report module."""

import os
from datetime import datetime
from functools import cached_property
from pathlib import Path

import yaml

from .utils import path_safe

ROOT_DIR = Path(__file__).parent
GITHUB_URL = 'https://github.com/maelyg/ont_amplicon'


class Config:

    # METADATA_FILE = 'index.csv'
    TIMESTAMP_FILE = '*_start_timestamp.txt'
    VERSIONS_PATH = ROOT_DIR.parents[1] / 'versions.yml'
    DEFAULT_PARAMS_PATH = ROOT_DIR.parents[1] / 'params/default_params.yml'

    class SCHEMA:
        BLAST_HITS_FIELD_CSV = ROOT_DIR / 'schema/blast_hits_columns.csv'

    class OUTPUTS:
        REPORT_FILE_TEMPLATE = '{sample_id}_report.html'
        BAM_HTML_FILE_TEMPLATE = '{sample_id}_bam-alignment.html'

    class REPORT:
        TITLE = "Amplicon sequencing assembly report"
        SUBTITLE = (
            "Results generated through the"
            f' <a href="{GITHUB_URL}" target="_blank">'
            'ONT-amplicon NextFlow workflow</a>.')

    class CRITERIA:
        MIN_RAW_READS = 2500
        MIN_FILTERED_READS = 800

    @property
    def default_params(self) -> dict[str, str]:
        """Return dict of default workflow parameters."""
        with self.DEFAULT_PARAMS_PATH.open() as f:
            return yaml.safe_load(f)

    @property
    def result_dir(self) -> Path:
        """Ensure that result dir is propagated between Config instances."""
        return Path(os.environ['RESULT_DIR'])

    @property
    def report_path(self) -> Path:
        path = self.result_dir / self.OUTPUTS.REPORT_FILE_TEMPLATE.format(
            sample_id=path_safe(self.sample_id),
        )
        return path.absolute()

    @property
    def bam_html_path(self) -> Path:
        path = self.result_dir / self.OUTPUTS.BAM_HTML_FILE_TEMPLATE.format(
            sample_id=path_safe(self.sample_id),
        )
        return path.absolute()

    @property
    def metadata_path(self) -> Path:
        return self.result_dir / self.METADATA_FILE

    @property
    def run_qc_path(self) -> Path:
        return self._get_file_by_pattern('run_qc_report_*.txt')

    @property
    def run_qc_html_file(self) -> str:
        return self._get_file_by_pattern('run_qc_report.html').name

    @property
    def nanoplot_raw_html_path(self) -> Path:
        return self._get_file_by_pattern('*raw_nanoplot-report.html')

    @property
    def nanoplot_filtered_html_path(self) -> Path:
        return self._get_file_by_pattern('*filtered_nanoplot-report.html')

    @property
    def bam_path(self) -> Path:
        return self._get_file_by_pattern("*.bam")

    @property
    def bai_path(self) -> Path:
        return self._get_file_by_pattern("*.bai")

    @property
    def consensus_fasta_path(self) -> Path:
        return self._get_file_by_pattern("*polished_consensus.fasta")
    
    @property
    def consensus_match_fasta_path(self) -> Path:
        return self._get_file_by_pattern("*polished_consensus_match.fasta")

    @property
    def blast_hits_path(self) -> Path:
        return self._get_file_by_pattern("*top_blast_with_cov_stats.txt")

    @property
    def blast_hits_polished_path(self) -> Path:
        return self._get_file_by_pattern(
            "*polished_consensus_rc_megablast_top_10_hits.txt")

    @cached_property
    def sample_id(self) -> str:
        """Return the sample ID from the result directory."""
        bam_path = self._get_file_by_pattern('*.bam')
        return bam_path.name.split('_aln.')[0]

    @cached_property
    def start_time(self) -> str:
        """Return the timestamp of the start of the workflow."""
        timestamp_path = self._get_file_by_pattern(self.TIMESTAMP_FILE)
        if timestamp_path.exists():
            return datetime.strptime(
                timestamp_path.read_text().strip(),
                '%Y%m%d%H%M%S',
            )
        return None

    def load(self, result_dir: Path):
        """Read the results from the result directory."""
        os.environ['RESULT_DIR'] = str(result_dir)

    def _get_file_by_pattern(self, file_pattern: str) -> Path:
        paths = list(self.result_dir.glob(file_pattern, case_sensitive=False))
        if paths:
            return paths[0]
        raise FileNotFoundError(
            f'No file matching pattern: {self.result_dir / file_pattern}'
        )
