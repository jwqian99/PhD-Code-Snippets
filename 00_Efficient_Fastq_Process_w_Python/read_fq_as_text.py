from typing import Union


def get_file_suffix(file_name):
    return file_name.split('.')[-1]


def load_fastq_as_text(fastq_file, return_list=True):
    import mmap

    if get_file_suffix(fastq_file) in ['gz']:
        import gzip

        with open(fastq_file, mode='rb') as f_binary:
            with mmap.mmap(f_binary.fileno(), length=0, access=mmap.ACCESS_READ) as mmap_obj:
                with gzip.open(mmap_obj, mode='r') as g_read:
                    reads = g_read.read().decode('utf-8').strip()  # add strip to exclude the trailing \n
    elif get_file_suffix(fastq_file) in ['fastq', 'fq']:
        with open(fastq_file, mode='rt', encoding='utf-8') as f_text:
            with mmap.mmap(f_text.fileno(), length=0, access=mmap.ACCESS_READ) as mmap_obj:
                reads = mmap_obj.read().decode('utf-8').strip()  # add strip to exclude the trailing \n

    if return_list:
        reads = reads.split('\n')
        return reads
    else:
        return reads


def load_fastq_list_as_dataframe(fastq_list: Union[list, list[list]], mode='single'):
    import pandas as pd
    from tqdm import tqdm

    # logic: generate a list -> load the list directly to pd.DataFrame

    if mode == 'single':
        total_reads = int(len(fastq_list) / 4)  # this list contains lines from the single fastq file

        return pd.DataFrame(
            [
                (
                    fastq_list[i * 4],
                    fastq_list[i * 4 + 1],
                    fastq_list[i * 4 + 3]
                ) for i in tqdm(range(total_reads), total=total_reads)
            ],
            columns=[
                'read_1_header',  # add 1 index to make compatible with other code
                'read_1_sequence',
                'read_1_score'
            ]
        )
    elif mode == 'paired':
        total_reads = int(len(fastq_list[0]) / 4)  # the first one is read 1, and the second one is read 2.

        return pd.DataFrame(
            [
                (
                    fastq_list[0][i * 4], fastq_list[1][i * 4],
                    fastq_list[0][i * 4 + 1], fastq_list[1][i * 4 + 1],
                    fastq_list[0][i * 4 + 3], fastq_list[1][i * 4 + 3]
                ) for i in tqdm(range(total_reads), total=total_reads)
            ],
            columns=[
                'read_1_header', 'read_2_header',
                'read_1_sequence', 'read_2_sequence',
                'read_1_score', 'read_2_score'
            ]
        )


def save_fastq_df_to_np_array(header_as_list, sequence_as_list, score_as_list):
    import numpy as np

    return np.array(
        list(
            zip(
                header_as_list,
                sequence_as_list,
                ['+'] * len(header_as_list),
                score_as_list,
            )
        )
    ).reshape(-1)


def convert_fastq_np_array_to_string(fq_array):
    return '\n'.join(fq_array) + '\n'


def save_with_pgzip(to_save_string, filepath, thread=1, blocksize=2 * 10 ** 8):
    import pgzip

    with pgzip.open(filepath, mode='wt', thread=thread, blocksize=blocksize) as output_file_stream:
        output_file_stream.write(to_save_string)
