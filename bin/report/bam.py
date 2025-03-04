import logging
from jinja2 import Environment, FileSystemLoader
from pathlib import Path

from .config import Config

config = Config()
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

TEMPLATE_DIR = Path(__file__).parent / 'templates'
STATIC_DIR = Path(__file__).parent / 'static'
TEMPLATE_NAME = "bam-viewer.html"


def file_to_js_array(path):
    byte_array = path.read_bytes()
    return "[" + ",".join(str(b) for b in byte_array) + "]"


def render_bam_html():
    j2 = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
    template = j2.get_template(TEMPLATE_NAME)
    context = {
        k: file_to_js_array(v)
        for k, v in [
            ('bam_binary_arr', config.bam_path),
            ('bai_binary_arr', config.bai_path),
            ('fasta_binary_arr', config.consensus_fasta_path),
        ]
    }
    context.update({
        'sample_id': config.sample_id,
        'loading_svg': (STATIC_DIR / 'img/spinner.svg').read_text(),
        'igv_js': (STATIC_DIR / 'js/igv.min.js').read_text(),
        'bootstrap_css': (STATIC_DIR / 'css/bootstrap.min.css').read_text(),
    })
    rendered_html = template.render(**context)
    path = config.bam_html_output_path
    path.write_text(rendered_html)
    logger.info(f"BAM Viewer HTML generated: {path}")
