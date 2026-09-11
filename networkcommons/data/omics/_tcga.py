#!/usr/bin/env python

#
# This file is part of the `networkcommons` Python module
#
# Copyright 2024
# Heidelberg University Hospital
#
# File author(s): Saez Lab (omnipathdb@gmail.com)
#
# Distributed under the GPLv3 license
# See the file `LICENSE` or read a copy at
# https://www.gnu.org/licenses/gpl-3.0.txt
#

"""
Access to open TCGA data through the Genomic Data Commons API.
"""

from __future__ import annotations

__all__ = [
    'tcga_projects',
    'tcga_cases',
    'tcga_files',
    'tcga_datatypes',
    'tcga_table',
    'tcga_rppa_files',
    'tcga_download',
    'tcga_query',
]

import os
import re
from collections.abc import Sequence

import pandas as pd

from . import _common
from networkcommons import _conf
from networkcommons._session import _log


GDC_API = 'https://api.gdc.cancer.gov'
_ENDPOINTS = {'projects', 'cases', 'files', 'annotations'}


def _values(value: str | Sequence[str] | None) -> list[str] | None:

    if value is None:

        return None

    return [value] if isinstance(value, str) else list(value)


def _filter(field: str, value: str | Sequence[str] | None) -> dict | None:

    values = _values(value)

    if not values:

        return None

    op = '=' if len(values) == 1 else 'in'

    return {'op': op, 'content': {'field': field, 'value': values}}


def _and(*filters: dict | None) -> dict | None:

    filters = [f for f in filters if f]

    if not filters:

        return None

    return filters[0] if len(filters) == 1 else {'op': 'and', 'content': filters}


def _fields(fields: str | Sequence[str] | None) -> str | None:

    if fields is None or isinstance(fields, str):

        return fields

    return ','.join(fields)


def tcga_query(
        endpoint: str,
        filters: dict | None = None,
        fields: str | Sequence[str] | None = None,
        size: int = 100,
        from_: int = 0,
    ) -> pd.DataFrame:
    """
    Query a GDC search endpoint and return normalized hits.
    """

    endpoint = endpoint.strip('/')

    if endpoint not in _ENDPOINTS:

        raise ValueError(f'Unsupported GDC endpoint: {endpoint}.')

    payload = {
        'format': 'JSON',
        'size': size,
        'from': from_,
    }

    if filters:

        payload['filters'] = filters

    if fields := _fields(fields):

        payload['fields'] = fields

    _log(f'DATA: Querying GDC {endpoint} endpoint...')

    resp = _common._requests_session().post(f'{GDC_API}/{endpoint}', json=payload)
    resp.raise_for_status()

    return pd.json_normalize(resp.json().get('data', {}).get('hits', []))


def tcga_projects(size: int = 1000) -> pd.DataFrame:
    """
    List TCGA projects from GDC.
    """

    return tcga_query(
        'projects',
        filters=_filter('program.name', 'TCGA'),
        size=size,
    )


def tcga_cases(
        project_id: str | Sequence[str] | None = None,
        fields: str | Sequence[str] | None = None,
        size: int = 100,
        from_: int = 0,
    ) -> pd.DataFrame:
    """
    List TCGA cases from GDC.
    """

    return tcga_query(
        'cases',
        filters=_and(
            _filter('cases.project.program.name', 'TCGA'),
            _filter('project.project_id', project_id),
        ),
        fields=fields,
        size=size,
        from_=from_,
    )


def tcga_files(
        project_id: str | Sequence[str] | None = None,
        data_category: str | Sequence[str] | None = None,
        data_type: str | Sequence[str] | None = None,
        experimental_strategy: str | Sequence[str] | None = None,
        workflow_type: str | Sequence[str] | None = None,
        fields: str | Sequence[str] | None = None,
        size: int = 100,
        from_: int = 0,
    ) -> pd.DataFrame:
    """
    List open-access TCGA files from GDC.
    """

    return tcga_query(
        'files',
        filters=_and(
            _filter('cases.project.program.name', 'TCGA'),
            _filter('cases.project.project_id', project_id),
            _filter('access', 'open'),
            _filter('data_category', data_category),
            _filter('data_type', data_type),
            _filter('experimental_strategy', experimental_strategy),
            _filter('analysis.workflow_type', workflow_type),
        ),
        fields=fields,
        size=size,
        from_=from_,
    )


def tcga_datatypes() -> pd.DataFrame:
    """
    Common open TCGA data types available through GDC.
    """

    return pd.DataFrame({
        'type': ['rppa', 'rnaseq', 'mirna', 'mutation', 'methylation', 'cnv'],
        'data_category': [
            'Proteome Profiling',
            'Transcriptome Profiling',
            'Transcriptome Profiling',
            'Simple Nucleotide Variation',
            'DNA Methylation',
            'Copy Number Variation',
        ],
        'data_type': [
            'Protein Expression Quantification',
            'Gene Expression Quantification',
            'miRNA Expression Quantification',
            'Masked Somatic Mutation',
            'Methylation Beta Value',
            'Gene Level Copy Number',
        ],
        'experimental_strategy': [
            'Reverse Phase Protein Array',
            'RNA-Seq',
            'miRNA-Seq',
            'WXS',
            'Methylation Array',
            'Genotyping Array',
        ],
        'description': [
            'Reverse phase protein array protein expression',
            'RNA-seq gene expression quantification',
            'miRNA expression quantification',
            'Masked somatic mutation calls',
            'DNA methylation beta values',
            'Gene-level copy number estimates',
        ],
    })


def tcga_table(file_id: str, sep: str = '\t') -> pd.DataFrame:
    """
    One TCGA data table from GDC.
    """

    _log(f'DATA: Retrieving TCGA table {file_id}...')

    return pd.read_csv(
        tcga_download(file_id, path=_conf.get('cachedir')),
        sep=sep,
    )


def tcga_rppa_files(
        project_id: str | Sequence[str] | None = None,
        fields: str | Sequence[str] | None = None,
        size: int = 100,
        from_: int = 0,
    ) -> pd.DataFrame:
    """
    List open-access TCGA RPPA files from GDC.
    """

    return tcga_files(
        project_id=project_id,
        data_category='Proteome Profiling',
        data_type='Protein Expression Quantification',
        experimental_strategy='Reverse Phase Protein Array',
        fields=fields,
        size=size,
        from_=from_,
    )


def tcga_download(file_id: str, path: str = '.') -> str:
    """
    Download one open-access GDC file.
    """

    _log(f'DATA: Downloading GDC file {file_id}...')

    resp = _common._requests_session().get(f'{GDC_API}/data/{file_id}', stream=True)
    resp.raise_for_status()

    match = re.search(r'filename="?([^";]+)"?', resp.headers.get('Content-Disposition', ''))
    filename = match.group(1) if match else file_id
    outpath = os.path.join(path, filename)

    os.makedirs(path, exist_ok=True)

    with open(outpath, 'wb') as fp:

        for chunk in resp.iter_content(chunk_size=8192):

            if chunk:

                fp.write(chunk)

    return outpath
