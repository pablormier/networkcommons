import pytest

import pandas as pd
import anndata as ad

from networkcommons.data.omics import _common
from networkcommons.data import omics

from unittest.mock import patch, MagicMock, mock_open
import zipfile
import bs4

import responses
import contextlib


# FILE: omics/_common.py
def test_datasets():

    dsets = _common._datasets()

    assert 'baseurl' in dsets
    assert isinstance(dsets['datasets'], dict)
    assert 'decryptm' in dsets['datasets']


def test_datasets_2():

    dsets = _common.datasets()

    assert isinstance(dsets, pd.DataFrame)
    assert dsets.columns.tolist() == ['name', 'description', 'publication_link', 'detailed_description']
    assert 'decryptm' in dsets.index
    assert 'CPTAC' in dsets.index


def test_commons_url():

    url = _common._commons_url('test', table = 'meta')

    assert 'metadata' in url


@pytest.mark.slow
def test_open():

    url = _common._commons_url('test', table = 'meta')

    with _common._open(url) as fp:

        line = next(fp)

    assert line.startswith('sample_ID\t')


@pytest.mark.slow
def test_open_df():

    url = _common._commons_url('test', table = 'meta')
    df = _common._open(url, df = {'sep': '\t'})

    assert isinstance(df, pd.DataFrame)
    assert df.shape == (4, 2)


@patch('networkcommons.data.omics._common._maybe_download')
@patch('pandas.read_csv')
def test_open_with_pandas_readers(mock_csv, mock_download):
    mock_download.return_value = 'test.csv'
    ftype = 'csv'
    _common._open('http://example.com/test.csv', ftype, df=True)

    mock_download.assert_called_once_with('http://example.com/test.csv')

    mock_csv.assert_called_once_with('test.csv')


def test_open_tsv():
    url = "http://example.com/test.tsv"
    with patch('networkcommons.data.omics._common._maybe_download', return_value='path/to/test.tsv'), \
         patch('builtins.open', mock_open(read_data="col1\tcol2\nval1\tval2")):
        with _common._open(url, ftype='tsv') as f:
            content = f.read()
            assert "col1\tcol2\nval1\tval2" in content


def test_open_html():
    url = "http://example.com/test.html"
    with patch('networkcommons.data.omics._common._maybe_download', return_value='path/to/test.html'), \
         patch('builtins.open', mock_open(read_data="<html><body>Test</body></html>")):
        result = _common._open(url, ftype='html')
        assert isinstance(result, bs4.BeautifulSoup)
        assert result.body.text == "Test"


@patch('networkcommons.data.omics._common._maybe_download')
@patch('contextlib.closing')
@patch('zipfile.ZipFile')
def test_open_zip(mock_zip, contextlib_mock, mock_maybe_download):
    url = "http://example.com/test.zip"
    mock_maybe_download.return_value = 'path/to/test.zip'
    mock_zip.return_value = MagicMock()

    result = _common._open(url, ftype='zip')
    mock_zip.assert_called_once_with('path/to/test.zip', 'r')
    contextlib_mock.assert_called_once_with(mock_zip.return_value)


@patch('networkcommons.data.omics._common._download')
@patch('networkcommons.data.omics._common._log')
@patch('networkcommons.data.omics._common._conf.get')
@patch('os.path.exists')
@patch('hashlib.md5')
def test_maybe_download_exists(mock_md5, mock_exists, mock_conf_get, mock_log, mock_download):
    # Setup mock values
    url = 'http://example.com/file.txt'
    md5_hash = MagicMock()
    md5_hash.hexdigest.return_value = 'dummyhash'
    mock_md5.return_value = md5_hash
    mock_conf_get.return_value = '/mock/cache/dir'
    mock_exists.return_value = True

    # Call the function
    path = _common._maybe_download(url)

    # Assertions
    mock_md5.assert_called_once_with(url.encode())
    mock_conf_get.assert_called_once_with('cachedir')
    mock_exists.assert_called_once_with('/mock/cache/dir/dummyhash-file.txt')
    mock_log.assert_called_once_with('Utils: Looking up in cache: `http://example.com/file.txt` -> `/mock/cache/dir/dummyhash-file.txt`.')
    mock_download.assert_not_called()
    assert path == '/mock/cache/dir/dummyhash-file.txt'


@patch('networkcommons.data.omics._common._download')
@patch('networkcommons.data.omics._common._log')
@patch('networkcommons.data.omics._common._conf.get')
@patch('os.path.exists')
@patch('hashlib.md5')
def test_maybe_download_not_exists(mock_md5, mock_exists, mock_conf_get, mock_log, mock_download):
    # Setup mock values
    url = 'http://example.com/file.txt'
    md5_hash = MagicMock()
    md5_hash.hexdigest.return_value = 'dummyhash'
    mock_md5.return_value = md5_hash
    mock_conf_get.return_value = '/mock/cache/dir'
    mock_exists.return_value = False

    # Call the function
    path = _common._maybe_download(url)

    # Assertions
    mock_md5.assert_called_once_with(url.encode())
    mock_conf_get.assert_called_once_with('cachedir')
    mock_exists.assert_called_once_with('/mock/cache/dir/dummyhash-file.txt')
    mock_log.assert_any_call('Utils: Looking up in cache: `http://example.com/file.txt` -> `/mock/cache/dir/dummyhash-file.txt`.')
    mock_log.assert_any_call('Utils: Not found in cache, initiating download: `http://example.com/file.txt`.')
    mock_download.assert_called_once_with(url, '/mock/cache/dir/dummyhash-file.txt')
    assert path == '/mock/cache/dir/dummyhash-file.txt'


@patch('networkcommons.data.omics._common._requests_session')
@patch('networkcommons.data.omics._common._log')
@patch('networkcommons.data.omics._common._conf.get')
def test_download(mock_conf_get, mock_log, mock_requests_session, tmp_path):
    # Setup mock values
    url = 'http://example.com/file.txt'
    path = tmp_path / 'file.txt'
    timeouts = (5, 5)
    mock_conf_get.side_effect = lambda k: 5 if k in ('http_read_timout', 'http_connect_timout') else None
    mock_session = MagicMock()
    mock_requests_session.return_value = mock_session
    mock_response = MagicMock()
    mock_response.iter_content.return_value = [b'test content']
    mock_session.get.return_value.__enter__.return_value = mock_response

    # Call the function
    _common._download(url, str(path))

    # Assertions
    mock_conf_get.assert_any_call('http_read_timout')
    mock_conf_get.assert_any_call('http_connect_timout')
    mock_log.assert_any_call(f'Utils: Downloading `{url}` to `{path}`.')
    mock_log.assert_any_call(f'Utils: Finished downloading `{url}` to `{path}`.')
    mock_requests_session.assert_called_once()
    mock_session.get.assert_called_once_with(url, timeout=(5, 5), stream=True)
    mock_response.raise_for_status.assert_called_once()
    mock_response.iter_content.assert_called_once_with(chunk_size=8192)

    # Check that the file was written correctly
    with open(path, 'rb') as f:
        content = f.read()
    assert content == b'test content'


def test_ls_success():
    url = "http://example.com/dir/"
    html_content = '''
    <html>
        <body>
            <a href="file1.txt">file1.txt</a>
            <a href="file2.txt">file2.txt</a>
            <a href="..">parent</a>
        </body>
    </html>
    '''

    with responses.RequestsMock() as rsps:
        rsps.add(responses.GET, url, body=html_content, status=200)
        result = _common._ls(url)
        assert result == ["file1.txt", "file2.txt"]


def test_ls_not_found():
    url = "http://example.com/dir/"

    with responses.RequestsMock() as rsps:
        rsps.add(responses.GET, url, status=404)
        with pytest.raises(FileNotFoundError, match="URL http://example.com/dir/ returned status code 404"):
            _common._ls(url)


@patch('networkcommons.data.omics._common._maybe_download')
def test_open_unknown_file_type(mock_maybe_download):
    url = 'http://example.com/file.unknown'
    mock_maybe_download.return_value = 'file.unknown'
    with pytest.raises(NotImplementedError, match='Can not open file type `unknown`.'):
        _common._open(url, 'unknown')


@patch('networkcommons.data.omics._common._maybe_download')
def test_open_no_extension(mock_maybe_download):
    url = 'http://example.com/file'
    mock_maybe_download.return_value = 'file'
    with pytest.raises(RuntimeError, match='Cannot determine file type for http://example.com/file.'):
        _common._open(url)


# FILE: omics/_decryptm.py
@pytest.fixture
def decryptm_args():
    return 'KDAC_Inhibitors', 'Acetylome', 'curves_CUDC101.txt'


@patch('networkcommons.data.omics._decryptm._common._ls')
@patch('networkcommons.data.omics._decryptm._common._baseurl', return_value='http://example.com')
@patch('pandas.read_pickle')
@patch('os.path.exists', return_value=False)
@patch('pandas.DataFrame.to_pickle')
def test_decryptm_datasets_update(mock_to_pickle, mock_path_exists, mock_read_pickle, mock_baseurl, mock_ls):
    # Mock the directory listing
    mock_ls.side_effect = [
        ['experiment1', 'experiment2'],  # First call, list experiments
        ['data_type1', 'data_type2'],  # Second call, list data types for experiment1
        ['curves_file1.txt', 'curves_file2.txt'],  # Third call, list files for experiment1/data_type1
        ['curves_file3.txt', 'curves_file4.txt'],  # Fourth call, list files for experiment1/data_type2
        ['data_type1', 'data_type2'],  # Fifth call, list data types for experiment2
        ['curves_file5.txt', 'curves_file6.txt'],  # Sixth call, list files for experiment2/data_type1
        ['curves_file7.txt', 'curves_file8.txt']   # Seventh call, list files for experiment2/data_type2
    ]

    dsets = omics.decryptm_datasets(update=True)

    assert isinstance(dsets, pd.DataFrame)
    assert dsets.shape == (8, 3)  # 4 experiments * 2 data types = 8 files
    assert dsets.columns.tolist() == ['experiment', 'data_type', 'fname']
    mock_to_pickle.assert_called_once()


@patch('pandas.read_pickle')
@patch('os.path.exists', return_value=True)
def test_decryptm_datasets_cached(mock_path_exists, mock_read_pickle):
    # Mock the cached DataFrame
    mock_df = pd.DataFrame({
        'experiment': ['experiment1', 'experiment2'],
        'data_type': ['data_type1', 'data_type2'],
        'fname': ['curves_file1.txt', 'curves_file2.txt']
    })
    mock_read_pickle.return_value = mock_df

    dsets = omics.decryptm_datasets(update=False)

    assert isinstance(dsets, pd.DataFrame)
    assert dsets.shape == (2, 3)
    assert dsets.columns.tolist() == ['experiment', 'data_type', 'fname']
    mock_read_pickle.assert_called_once()


@patch('networkcommons.data.omics._decryptm._common._open')
def test_decryptm_table(mock_open, decryptm_args):
    mock_df = pd.DataFrame({'EC50': [0.5, 1.0, 1.5]})
    mock_open.return_value = mock_df

    df = omics.decryptm_table(*decryptm_args)

    assert isinstance(df, pd.DataFrame)
    assert df.shape == (3, 1)
    assert df.EC50.dtype == 'float64'
    mock_open.assert_called_once()


@patch('networkcommons.data.omics._decryptm.decryptm_datasets')
@patch('networkcommons.data.omics._decryptm.decryptm_table')
def test_decryptm_experiment(mock_decryptm_table, mock_decryptm_datasets, decryptm_args):
    mock_decryptm_datasets.return_value = pd.DataFrame({
        'experiment': ['KDAC_Inhibitors', 'KDAC_Inhibitors'],
        'data_type': ['Acetylome', 'Acetylome'],
        'fname': ['curves_CUDC101.txt', 'curves_other.txt']
    })
    mock_df = pd.DataFrame({'EC50': [0.5, 1.0, 1.5]})
    mock_decryptm_table.return_value = mock_df

    dfs = omics.decryptm_experiment(decryptm_args[0], decryptm_args[1])

    assert isinstance(dfs, list)
    assert len(dfs) == 2
    assert all(isinstance(df, pd.DataFrame) for df in dfs)
    assert dfs[0].shape == (3, 1)
    assert dfs[0].EC50.dtype == 'float64'
    mock_decryptm_table.assert_called()


@patch('networkcommons.data.omics._decryptm.decryptm_datasets')
def test_decryptm_experiment_no_dataset(mock_decryptm_datasets):
    mock_decryptm_datasets.return_value = pd.DataFrame({
        'experiment': ['KDAC_Inhibitors'],
        'data_type': ['Acetylome'],
        'fname': ['curves_CUDC101.txt']
    })

    with pytest.raises(ValueError, match='No such dataset in DecryptM: `Invalid_Experiment/Invalid_Type`.'):
        omics.decryptm_experiment('Invalid_Experiment', 'Invalid_Type')


# FILE: omics/_panacea.py

@patch('urllib.request.urlopen')
@patch('pandas.read_csv')
@patch('os.path.exists', return_value=False)
@patch('pandas.DataFrame.to_pickle')
def test_panacea_experiments(mock_to_pickle, mock_path_exists, mock_read_csv, mock_urlopen):
    # Mock the metadata file (panacea__metadata.tsv)
    mock_metadata = pd.DataFrame({
        'group': ['ASPC_AEE788', 'ASPC_AFATINIB', 'DU145_CRIZOTINIB', 'ASPC_CEDIRANIB'],
        'sample_ID': ['ID1', 'ID2', 'ID3', 'ID4']
    })
    mock_read_csv.return_value = mock_metadata

    # Mock the HTML content returned by the remote server, 
    # with some lines containing the "__TF_scores.tsv" pattern
    mock_html_content = """
        <html>
        <body>
            <a href="ASPC_AEE788__TF_scores.tsv">ASPC_AEE788__TF_scores.tsv</a><br>
            <a href="ASPC_AFATINIB__TF_scores.tsv">ASPC_AFATINIB__TF_scores.tsv</a><br>
            <a href="DU145_CRIZOTINIB__TF_scores.tsv">DU145_CRIZOTINIB__TF_scores.tsv</a><br>
        </body>
        </html>
    """
    
    # Mock the HTTP response for the TF scores directory
    mock_response = MagicMock()
    mock_response.read.return_value = mock_html_content.encode('utf-8')
    mock_urlopen.return_value = mock_response

    # Call the function under test
    result_df = omics.panacea_experiments(update=True)

    # Check that the result is a DataFrame
    assert isinstance(result_df, pd.DataFrame)

    # Check that the expected columns are in the DataFrame
    assert 'group' in result_df.columns
    assert 'cell' in result_df.columns
    assert 'drug' in result_df.columns
    assert 'tf_scores' in result_df.columns

    # Check that the 'tf_scores' column has the correct True/False values
    expected_tf_scores = [True, True, True, False]  # CEDIRANIB should not have TF scores
    assert result_df['tf_scores'].tolist() == expected_tf_scores

    # Verify that the function attempts to save the result as a pickle file
    mock_to_pickle.assert_called_once()


@patch('networkcommons.data.omics._panacea._common._baseurl', return_value='http://example.com')
@patch('pandas.read_pickle')
@patch('os.path.exists', return_value=True)
def test_panacea_experiments_cached(mock_path_exists, mock_read_pickle, mock_baseurl):
    # Mock the cached data
    mock_df = pd.DataFrame({'cell': ['A', 'C'], 'drug': ['B', 'D']})
    mock_read_pickle.return_value = mock_df

    result_df = omics.panacea_experiments(update=False)

    mock_read_pickle.assert_called_once()
    assert result_df.equals(mock_df)


def test_panacea_datatypes():
    dtypes = omics.panacea_datatypes()

    expected_df = pd.DataFrame({
        'type': ['raw', 'diffexp', 'TF_scores'],
        'description': [
            'RNA-Seq raw counts and metadata containing sample, name, and group',
            'Differential expression analysis with filterbyExpr+DESeq2',
            'Transcription factor activity scores with CollecTRI + T-values'
        ]
    })

    pd.testing.assert_frame_equal(dtypes, expected_df)


@patch('pandas.read_csv')
@patch('networkcommons.data.omics._panacea._common._baseurl', return_value='http://example.com')
def test_panacea_tables_diffexp(mock_baseurl, mock_read_csv):
    # Mock the data
    mock_df = pd.DataFrame({
        'gene': ['gene1', 'gene2'],
        'log2FoldChange': [1.5, -2.3],
        'pvalue': [0.01, 0.05]
    })
    mock_read_csv.return_value = mock_df

    result_df = omics.panacea_tables(cell_line='cell1', drug='drug1', type='diffexp')

    assert isinstance(result_df, pd.DataFrame)
    assert 'gene' in result_df.columns
    assert 'log2FoldChange' in result_df.columns
    assert 'pvalue' in result_df.columns


@patch('networkcommons.data.omics._panacea._common._open')
@patch('networkcommons.data.omics._panacea._common._baseurl', return_value='http://example.com')
def test_panacea_tables_convert_to_list(mock_baseurl, mock_open):
    # Mock the metadata
    mock_meta = pd.DataFrame({
        'sample_ID': ['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6'],
        'group': ['cell1_drug1', 'cell1_drug2', 'cell2_drug1', 'cell2_drug2', 'cell1_drug1', 'cell1_drug2']
    })
    # Mock the count data
    mock_count = pd.DataFrame({
        'gene_symbol': ['gene1', 'gene2'],
        'sample1': [100, 200],
        'sample2': [150, 250],
        'sample3': [100, 200],
        'sample4': [150, 250],
        'sample5': [100, 200],
        'sample6': [150, 250]
    })
    mock_open.side_effect = [mock_meta, mock_count] * 5


    # Test with cell_line and drug as strings
    df_count, df_meta = omics.panacea_tables(cell_line='cell1', drug='drug1', type='raw')
    assert isinstance(df_count, pd.DataFrame)
    assert isinstance(df_meta, pd.DataFrame)
    assert df_count.shape == (2, 3)
    assert df_meta.shape == (2, 4)

    # Test with cell_line and drug as lists
    df_count, df_meta = omics.panacea_tables(cell_line=['cell1'], drug=['drug1'], type='raw')
    assert isinstance(df_count, pd.DataFrame)
    assert isinstance(df_meta, pd.DataFrame)
    assert df_count.shape == (2, 3)
    assert df_meta.shape == (2, 4)

    # Test with cell_line and drug both None
    df_count, df_meta = omics.panacea_tables(type='raw')
    assert isinstance(df_count, pd.DataFrame)
    assert isinstance(df_meta, pd.DataFrame)
    assert df_count.shape == (2, 7)
    assert df_meta.shape == (6, 4)

    # Test with cell_line as None and drug as string
    df_count, df_meta = omics.panacea_tables(cell_line=None, drug='drug1', type='raw')
    assert isinstance(df_count, pd.DataFrame)
    assert isinstance(df_meta, pd.DataFrame)
    assert df_count.shape == (2, 4)
    assert df_meta.shape == (3, 4)

    # Test with cell_line as string and drug as None
    df_count, df_meta = omics.panacea_tables(cell_line='cell1', drug=None, type='raw')
    assert isinstance(df_count, pd.DataFrame)
    assert isinstance(df_meta, pd.DataFrame)
    assert df_count.shape == (2, 5)
    assert df_meta.shape == (4, 4)

    # Test with unknown type to trigger the ValueError
    with pytest.raises(ValueError, match='Unknown data type: unknown_type'):
        omics.panacea_tables(cell_line='cell1', drug='drug1', type='unknown_type')


@patch('pandas.read_csv')
@patch('networkcommons.data.omics._panacea._common._baseurl', return_value='http://example.com')
def test_panacea_tables_diffexp(mock_baseurl, mock_read_csv):
    # Mock the data
    mock_df = pd.DataFrame({
        'gene': ['gene1', 'gene2'],
        'log2FoldChange': [1.5, -2.3],
        'pvalue': [0.01, 0.05]
    })
    mock_read_csv.return_value = mock_df

    result_df = omics.panacea_tables(cell_line='cell1', drug='drug1', type='diffexp')

    assert isinstance(result_df, pd.DataFrame)
    assert 'gene' in result_df.columns
    assert 'log2FoldChange' in result_df.columns
    assert 'pvalue' in result_df.columns


@patch('pandas.read_csv')
@patch('networkcommons.data.omics._panacea._common._baseurl', return_value='http://example.com')
def test_panacea_tables_tf_scores(mock_baseurl, mock_read_csv):
    # Mock the data
    mock_df = pd.DataFrame({
        'TF': ['TF1', 'TF2'],
        'score': [2.5, -1.3],
        'pvalue': [0.02, 0.07]
    })
    mock_read_csv.return_value = mock_df

    result_df = omics.panacea_tables(cell_line='cell1', drug='drug1', type='TF_scores')

    assert isinstance(result_df, pd.DataFrame)
    assert 'TF' in result_df.columns
    assert 'score' in result_df.columns
    assert 'pvalue' in result_df.columns


def test_panacea_tables_value_error():
    with pytest.raises(ValueError, match='Please specify cell line and drug.'):
        omics.panacea_tables(type='diffexp')


@patch('networkcommons.data.omics._panacea._common._open')
@patch('networkcommons.data.omics._panacea._common._baseurl', return_value='http://example.com')
def test_panacea_tables_raw(mock_baseurl, mock_open):
    cell_line = 'CellLine1'
    drug = 'Drug1'
    data_type = 'raw'

    # Mock the DataFrames returned by _common._open
    mock_meta_df = pd.DataFrame({'group': ['CellLine1_Drug1', 'CellLine2_Drug2'], 'sample_ID': ['ID1', 'ID2']})
    mock_count_df = pd.DataFrame({'gene_symbol': ['Gene1', 'Gene2'], 'ID1': [10, 20], 'ID2': [30, 40]})
    mock_open.side_effect = [mock_meta_df, mock_count_df]

    result_count_df, result_meta_df = omics.panacea_tables(cell_line=cell_line, drug=drug, type=data_type)

    assert isinstance(result_count_df, pd.DataFrame)
    assert 'gene_symbol' in result_count_df.columns
    assert 'ID1' in result_count_df.columns
    assert isinstance(result_meta_df, pd.DataFrame)
    assert 'group' in result_meta_df.columns
    assert 'sample_ID' in result_meta_df.columns
    mock_open.assert_called()


def test_panacea_tables_no_cell_line_drug():
    with pytest.raises(ValueError, match='Please specify cell line and drug.'):
        omics.panacea_tables(type='diffexp')


@patch('networkcommons.data.omics._panacea._common._open')
@patch('networkcommons.data.omics._panacea._common._baseurl', return_value='http://example.com')
def test_panacea_tables_unknown_type(mock_baseurl, mock_open):
    with pytest.raises(ValueError, match='Unknown data type: unknown.'):
        omics.panacea_tables(cell_line='CellLine1', drug='Drug1', type='unknown')


@patch('networkcommons.data.omics._panacea._log')
@patch('networkcommons.data.omics._panacea._common._open')
@patch('networkcommons.data.omics._panacea._common._baseurl', return_value='http://example.com')
@patch('pandas.read_pickle')
@patch('os.path.exists')
@patch('os.makedirs')
@patch('pandas.DataFrame.to_pickle')
@patch('urllib.request.urlopen')
def test_panacea_gold_standard_update(mock_urllib, mock_pandas_topickle, mock_makedirs, mock_os, mock_pickle, mock_baseurl, mock_open, mock_log):
    
    # Mocking the HTTP response to return CSV-like content
    mock_response = MagicMock()
    mock_response.read.return_value = b"cmpd,cmpd_id,target,rank\nGene1,GeneA,Target1,1\nGene2,GeneB,Target2,t2"
    mock_response.__enter__.return_value = mock_response
    mock_urllib.return_value = mock_response
    mock_os.return_value = False

    # Run the function with `update=True`
    assert omics.panacea_gold_standard(update=True).shape == (2, 4)
    
    # Check if logs are being called correctly
    mock_log.assert_any_call('DATA: Retrieving Panacea offtarget gold standard...')
    mock_log.assert_any_call('DATA: not found in cache, downloading from server...')

    # Run the function with `update=False` to load from the cache
    mock_pickle.return_value = pd.DataFrame({
        'cmpd': ['Gene1', 'Gene2'],
        'cmpd_id': ['GeneA', 'GeneB'],
        'target': ['Target1', 'Target2'],
        'rank': [1, 2]
    })

    mock_os.return_value = True

    assert omics.panacea_gold_standard(update=False).shape == (2, 4)
    
    mock_log.assert_any_call('DATA: Retrieving Panacea offtarget gold standard...')
    mock_log.assert_any_call('DATA: found in cache, loading...')



# FILE: omics/_scperturb.py
import pytest
from unittest.mock import patch, MagicMock
import json

from networkcommons.data.omics import _scperturb

@pytest.fixture
def mock_metadata():
    return {
        'files': {
            'entries': {
                'dataset1.h5ad': {'links': {'content': 'https://example.com/dataset1.h5ad'}},
                'dataset2.h5ad': {'links': {'content': 'https://example.com/dataset2.h5ad'}}
            }
        }
    }


@pytest.fixture
def mock_ann_data():
    return MagicMock(spec=ad.AnnData)


@patch('networkcommons.data.omics._scperturb._common._open')
@patch('networkcommons.data.omics._scperturb.json.loads')
def test_scperturb_metadata(mock_json_loads, mock_open, mock_metadata):
    mock_open.return_value = MagicMock()
    mock_json_loads.return_value = mock_metadata

    metadata = _scperturb.scperturb_metadata()
    assert metadata == mock_metadata
    mock_open.assert_called_once_with('https://zenodo.org/record/10044268', ftype='html')
    mock_json_loads.assert_called_once()


@patch('networkcommons.data.omics._scperturb.scperturb_metadata')
def test_scperturb_datasets(mock_scperturb_metadata, mock_metadata):
    mock_scperturb_metadata.return_value = mock_metadata

    datasets = _scperturb.scperturb_datasets()
    expected_datasets = {
        'dataset1.h5ad': 'https://example.com/dataset1.h5ad',
        'dataset2.h5ad': 'https://example.com/dataset2.h5ad'
    }
    assert datasets == expected_datasets
    mock_scperturb_metadata.assert_called_once()


@patch('networkcommons.data.omics._scperturb.scperturb_datasets')
@patch('networkcommons.data.omics._scperturb._common._maybe_download')
@patch('anndata.read_h5ad')
def test_scperturb(mock_read_h5ad, mock_maybe_download, mock_scperturb_datasets, mock_ann_data):
    mock_scperturb_datasets.return_value = {
        'dataset1.h5ad': 'https://example.com/dataset1.h5ad'
    }
    mock_maybe_download.return_value = 'path/to/dataset1.h5ad'
    mock_read_h5ad.return_value = mock_ann_data

    result = _scperturb.scperturb('dataset1.h5ad')
    assert result is mock_ann_data
    mock_scperturb_datasets.assert_called_once()
    mock_maybe_download.assert_called_once_with('https://example.com/dataset1.h5ad')
    mock_read_h5ad.assert_called_once_with('path/to/dataset1.h5ad')


@pytest.mark.slow
def test_scperturb_metadata_slow():

    m = omics.scperturb_metadata()

    assert isinstance(m, dict)
    assert len(m['files']['entries']) == 50
    assert m['versions'] == {'index': 4, 'is_latest': True}


@pytest.mark.slow
def test_scperturb_datasets_slow():

    example_url = (
        'https://zenodo.org/api/records/10044268/files/'
        'XuCao2023.h5ad/content'
    )
    dsets = omics.scperturb_datasets()

    assert isinstance(dsets, dict)
    assert len(dsets) == 50
    assert dsets['XuCao2023.h5ad'] == example_url


@pytest.mark.slow
def test_scperturb_slow():

    var_cols = ('ensembl_id', 'ncounts', 'ncells')
    adata = omics.scperturb('AdamsonWeissman2016_GSM2406675_10X001.h5ad')

    assert isinstance(adata, ad.AnnData)
    assert tuple(adata.var.columns) == var_cols
    assert 'UMI count' in adata.obs.columns
    assert adata.shape == (5768, 35635)


@patch('networkcommons.data.omics._cptac._conf.get')
@patch('os.path.exists', return_value=True)
@patch('pandas.read_pickle')
def test_cptac_cohortsize_cached(mock_read_pickle, mock_path_exists, mock_conf_get):
    # Mock configuration and data
    mock_conf_get.return_value = '/mock/path'
    mock_df = pd.DataFrame({
        "Cancer_type": ["BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC"],
        "Tumor": [122, 103, 110, 99, 108, 108, 110, 83, 105, 95],
        "Normal": [0, 80, 100, 0, 62, 99, 101, 20, 44, 18]
    })
    mock_read_pickle.return_value = mock_df

    # Run the function with the condition that the pickle file exists
    result_df = omics.cptac_cohortsize()

    # Check that the result is as expected
    mock_read_pickle.assert_called_once_with('/mock/path/cptac_cohort.pickle')
    pd.testing.assert_frame_equal(result_df, mock_df)


@patch('networkcommons.data.omics._cptac._conf.get')
@patch('os.makedirs')  # Patch os.makedirs to prevent FileNotFoundError
@patch('os.path.exists', return_value=False)
@patch('pandas.read_excel')
@patch('pandas.DataFrame.to_pickle')
def test_cptac_cohortsize_download(mock_to_pickle, mock_read_excel, mock_makedirs, mock_conf_get, mock_path_exists):
    # Mock configuration and data
    mock_conf_get.return_value = '/mock/path'
    mock_df = pd.DataFrame({
        "Cancer_type": ["BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC"],
        "Tumor": [122, 103, 110, 99, 108, 108, 110, 83, 105, 95],
        "Normal": [0, 80, 100, 0, 62, 99, 101, 20, 44, 18]
    })
    mock_read_excel.return_value = mock_df

    # Run the function with the condition that the pickle file does not exist
    result_df = omics.cptac_cohortsize(update=True)

    # Check that the result is as expected
    mock_read_excel.assert_called_once()
    mock_to_pickle.assert_called_once()
    pd.testing.assert_frame_equal(result_df, mock_df)


@patch('networkcommons.data.omics._cptac._conf.get')
@patch('os.path.exists', return_value=True)
@patch('pandas.read_pickle')
def test_cptac_fileinfo_cached(mock_read_pickle, mock_path_exists, mock_conf_get):
    # Mock configuration and data
    mock_conf_get.return_value = '/mock/path'
    mock_df = pd.DataFrame({
        "File name": ["file1.txt", "file2.txt"],
        "Description": ["Description1", "Description2"]
    })
    mock_read_pickle.return_value = mock_df

    # Run the function with the condition that the pickle file exists
    result_df = omics.cptac_fileinfo()

    # Check that the result is as expected
    mock_read_pickle.assert_called_once_with('/mock/path/cptac_info.pickle')
    pd.testing.assert_frame_equal(result_df, mock_df)


@patch('networkcommons.data.omics._cptac._conf.get')
@patch('os.makedirs')  # Patch os.makedirs to prevent FileNotFoundError
@patch('os.path.exists', return_value=False)
@patch('pandas.read_excel')
@patch('pandas.DataFrame.to_pickle')
def test_cptac_fileinfo_download(mock_to_pickle, mock_read_excel, mock_makedirs, mock_conf_get, mock_path_exists):
    # Mock configuration and data
    mock_conf_get.return_value = '/mock/path'
    mock_df = pd.DataFrame({
        "File name": ["file1.txt", "file2.txt"],
        "Description": ["Description1", "Description2"]
    })
    mock_read_excel.return_value = mock_df

    # Run the function with the condition that the pickle file does not exist
    result_df = omics.cptac_fileinfo(update=True)

    # Check that the result is as expected
    mock_read_excel.assert_called_once()
    mock_to_pickle.assert_called_once()
    pd.testing.assert_frame_equal(result_df, mock_df)


@patch('networkcommons.data.omics._cptac._common._ls')
@patch('networkcommons.data.omics._cptac._common._baseurl', return_value='http://example.com/')
def test_cptac_datatypes(mock_baseurl, mock_ls):
    # Mock the return value of _ls to simulate the directory listing
    mock_ls.return_value = [
        'directory1',
        'directory2',
        'CPTAC_pancancer_data_freeze_cohort_size.xlsx',
        'CPTAC_pancancer_data_freeze_file_description.xlsx'
    ]

    expected_directories = ['directory1', 'directory2']

    # Call the function
    directories = omics.cptac_datatypes()

    # Check if the returned directories match the expected directories
    assert directories == expected_directories


@patch('networkcommons.data.omics._common._open')
def test_cptac_table(mock_open):
    mock_df = pd.DataFrame({
        "sample_ID": ["sample1", "sample2"],
        "value": [123, 456]
    })
    mock_open.return_value = mock_df

    df = omics.cptac_table('proteomics', 'BRCA', 'file.tsv')

    assert isinstance(df, pd.DataFrame)
    assert df.shape == (2, 2)
    mock_open.assert_called_once_with(
        _common._commons_url('CPTAC', data_type='proteomics', cancer_type='BRCA', fname='file.tsv'),
        df={'sep': '\t'}
    )


def test_cptac_extend_dataframe():
    df = pd.DataFrame({
        "idx": ["sample1", "sample2", "sample3"],
        "Tumor": ["Yes", "No", "Yes"],
        "Normal": ["No", "Yes", "No"]
    })

    extended_df = omics.cptac_extend_dataframe(df)

    print(extended_df)

    expected_df = pd.DataFrame({
        "sample_ID": ["sample1_tumor", "sample3_tumor", "sample2_ctrl"]
    })

    pd.testing.assert_frame_equal(extended_df, expected_df)


@patch('networkcommons.data.omics._common._conf.get')
@patch('pandas.read_pickle')
@patch('os.path.exists', return_value=True)
def test_get_ensembl_mappings_cached(mock_path_exists, mock_read_pickle, mock_conf_get):
    # Mock configuration and data
    mock_conf_get.return_value = '/path/to/pickle/dir'
    mock_df = pd.DataFrame({
        'gene_symbol': ['BRCA2', 'BRCA1'],
        'ensembl_id': ['ENSG00000139618', 'ENSG00000012048']
    })
    mock_read_pickle.return_value = mock_df

    # Run the function with the condition that the pickle file exists
    result_df = _common.get_ensembl_mappings()

    # Check that the result is as expected
    mock_read_pickle.assert_called_once_with('/path/to/pickle/dir/ensembl_map.pickle')


@patch('networkcommons.data.omics._common._conf.get')
@patch('os.path.exists', return_value=False)
@patch('biomart.BiomartServer')
def test_get_ensembl_mappings_download(mock_biomart_server, mock_path_exists, mock_conf_get):
    # Mock configuration and data
    mock_conf_get.return_value = '/path/to/pickle/dir'

    # Mock the biomart server and dataset
    mock_server_instance = MagicMock()
    mock_biomart_server.return_value = mock_server_instance
    mock_dataset = mock_server_instance.datasets['hsapiens_gene_ensembl']
    mock_response = MagicMock()
    mock_dataset.search.return_value = mock_response
    mock_response.raw.data.decode.return_value = (
        'ENST00000361390\tBRCA2\tENSG00000139618\tENSP00000354687\n'
        'ENST00000361453\tBRCA2\tENSG00000139618\tENSP00000354687\n'
        'ENST00000361453\tBRCA1\tENSG00000012048\tENSP00000354688\n'
    )

    with patch('pandas.DataFrame.to_pickle') as mock_to_pickle:
        result_df = _common.get_ensembl_mappings()

        expected_data = {
            'gene_symbol': ['BRCA2', 'BRCA2', 'BRCA1', 'BRCA2', 'BRCA1', 'BRCA2', 'BRCA1'],
            'ensembl_id': ['ENST00000361390', 'ENST00000361453', 'ENST00000361453',
                        'ENSG00000139618', 'ENSG00000012048', 'ENSP00000354687',
                        'ENSP00000354688']
        }
        expected_df = pd.DataFrame(expected_data)

        pd.testing.assert_frame_equal(result_df.reset_index(drop=True), expected_df)
        mock_to_pickle.assert_called_once_with('/path/to/pickle/dir/ensembl_map.pickle')


def test_convert_ensembl_to_gene_symbol_max():
    dataframe = pd.DataFrame({
        'idx': ['ENSG000001.23', 'ENSG000002', 'ENSG000001.19'],
        'value': [10, 20, 15]
    })
    equivalence_df = pd.DataFrame({
        'ensembl_id': ['ENSG000001', 'ENSG000002'],
        'gene_symbol': ['GeneA', 'GeneB']
    })
    result_df = omics.convert_ensembl_to_gene_symbol(dataframe, equivalence_df, summarisation='max')
    expected_df = pd.DataFrame({
        'gene_symbol': ['GeneA', 'GeneB'],
        'value': [15, 20]
    })
    pd.testing.assert_frame_equal(result_df, expected_df)


def test_convert_ensembl_to_gene_symbol_min():
    dataframe = pd.DataFrame({
        'idx': ['ENSG000001.28', 'ENSG000002', 'ENSG000001.23'],
        'value': [10, 20, 15]
    })
    equivalence_df = pd.DataFrame({
        'ensembl_id': ['ENSG000001', 'ENSG000002'],
        'gene_symbol': ['GeneA', 'GeneB']
    })
    result_df = omics.convert_ensembl_to_gene_symbol(dataframe, equivalence_df, summarisation='min')
    expected_df = pd.DataFrame({
        'gene_symbol': ['GeneA', 'GeneB'],
        'value': [10, 20]
    })
    pd.testing.assert_frame_equal(result_df, expected_df)


def test_convert_ensembl_to_gene_symbol_mean():
    dataframe = pd.DataFrame({
        'idx': ['ENSG000001.29', 'ENSG000002', 'ENSG000001.48'],
        'value': [10, 20, 15]
    })
    equivalence_df = pd.DataFrame({
        'ensembl_id': ['ENSG000001', 'ENSG000002'],
        'gene_symbol': ['GeneA', 'GeneB']
    })
    result_df = omics.convert_ensembl_to_gene_symbol(dataframe, equivalence_df, summarisation='mean')
    expected_df = pd.DataFrame({
        'gene_symbol': ['GeneA', 'GeneB'],
        'value': [12.5, 20]
    })
    pd.testing.assert_frame_equal(result_df, expected_df)


def test_convert_ensembl_to_gene_symbol_median():
    dataframe = pd.DataFrame({
        'idx': ['ENSG000001.10', 'ENSG000002', 'ENSG000001.2'],
        'value': [10, 20, 15]
    }).set_index('idx')
    equivalence_df = pd.DataFrame({
        'ensembl_id': ['ENSG000001', 'ENSG000002'],
        'gene_symbol': ['GeneA', 'GeneB']
    })
    result_df = omics.convert_ensembl_to_gene_symbol(dataframe, equivalence_df, summarisation='median')
    expected_df = pd.DataFrame({
        'gene_symbol': ['GeneA', 'GeneB'],
        'value': [12.5, 20]
    })
    pd.testing.assert_frame_equal(result_df, expected_df)

@patch('networkcommons.data.omics._common._log')
def test_convert_ensembl_to_gene_symbol_no_match(mock_log):
    dataframe = pd.DataFrame({
        'idx': ['ENSG000001.1', 'ENSG000003', 'ENSG000001.02'],
        'value': [10, 20, 15]
    })
    equivalence_df = pd.DataFrame({
        'ensembl_id': ['ENSG000001', 'ENSG000002'],
        'gene_symbol': ['GeneA', 'GeneB']
    })
    with patch('builtins.print') as mocked_print:
        result_df = omics.convert_ensembl_to_gene_symbol(dataframe, equivalence_df, summarisation='mean')
        expected_df = pd.DataFrame({
            'gene_symbol': ['GeneA'],
            'value': [12.5]
        })
        pd.testing.assert_frame_equal(result_df, expected_df)
        mock_log.assert_any_call("Utils: Number of non-matched Ensembl IDs: 1 (33.33%)")


@patch('biomart.BiomartServer')
def test_get_ensembl_mappings(mock_biomart_server):
    # Mock the biomart server and dataset
    mock_server_instance = MagicMock()
    mock_biomart_server.return_value = mock_server_instance
    mock_dataset = mock_server_instance.datasets['hsapiens_gene_ensembl']
    mock_response = MagicMock()
    mock_dataset.search.return_value = mock_response
    mock_response.raw.data.decode.return_value = (
        'ENST00000361390\tBRCA2\tENSG00000139618\tENSP00000354687\n'
        'ENST00000361453\tBRCA2\tENSG00000139618\tENSP00000354687\n'
        'ENST00000361453\tBRCA1\tENSG00000012048\tENSP00000354688\n'
    )

    result_df = omics.get_ensembl_mappings()
    print(result_df)

    expected_data = {
        'gene_symbol': ['BRCA2', 'BRCA2', 'BRCA1', 'BRCA2', 'BRCA1', 'BRCA2', 'BRCA1'],
        'ensembl_id': ['ENST00000361390', 'ENST00000361453', 'ENST00000361453',
                       'ENSG00000139618', 'ENSG00000012048', 'ENSP00000354687',
                       'ENSP00000354688']
    }
    expected_df = pd.DataFrame(expected_data)

    pd.testing.assert_frame_equal(result_df.reset_index(drop=True), expected_df)


# FILE: omics/_nci60.py

@patch('networkcommons.data.omics._nci60._common._ls')
@patch('networkcommons.data.omics._nci60._common._baseurl', return_value='http://example.com')
@patch('pandas.read_pickle')
@patch('os.path.exists', return_value=False)
@patch('pandas.DataFrame.to_pickle')
def test_nci60_datasets(mock_to_pickle, mock_path_exists, mock_read_pickle, mock_baseurl, mock_ls):
    # Mock the directory listing
    mock_ls.return_value = ['cell_line1', 'cell_line2', 'cell_line3']

    dsets = omics.nci60_datasets(update=True)

    expected_df = pd.DataFrame({
        'cell_line': ['cell_line1', 'cell_line2', 'cell_line3']
    })
    pd.testing.assert_frame_equal(dsets, expected_df)
    mock_to_pickle.assert_called_once()
    mock_read_pickle.assert_not_called()


@patch('pandas.read_pickle')
@patch('os.path.exists', return_value=True)
def test_nci60_datasets_cached(mock_path_exists, mock_read_pickle):
    mock_df = pd.DataFrame({
        'cell_line': ['cell_line1', 'cell_line2', 'cell_line3']
    })
    mock_read_pickle.return_value = mock_df

    dsets = omics.nci60_datasets()

    pd.testing.assert_frame_equal(dsets, mock_df)
    mock_read_pickle.assert_called_once()
    

def test_nci60_datatypes():
    dtypes = omics.nci60_datatypes()

    expected_df = pd.DataFrame({
        'data_type': ['TF_scores', 'RNA', 'metabolomic'],
        'description': ['TF scores', 'RNA expression', 'metabolomic data']
    })

    pd.testing.assert_frame_equal(dtypes, expected_df)


@patch('networkcommons.data.omics._nci60._common._open')
def test_nci60_table(mock_open):
    cell_line = 'cell_line1'
    data_type = 'RNA'
    mock_df = pd.DataFrame({
        'gene': ['Gene1', 'Gene2'],
        'expression': [100, 200]
    })
    mock_open.return_value = mock_df

    result = omics.nci60_table(cell_line, data_type)

    pd.testing.assert_frame_equal(result, mock_df)
    mock_open.assert_called_once()