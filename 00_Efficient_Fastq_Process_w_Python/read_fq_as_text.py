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
