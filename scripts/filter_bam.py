import subprocess
import os
import sys
import contextlib
import pandas as pd
import click
import os.path as op
import os
import shutil
import six
import tempfile
from collections import defaultdict
import pysam

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option('0.4.5')
@click.pass_context
def cli(obj):
    """read recruitment."""
    pass


def remove_files(fnames):
    for x in fnames:
        if x and os.path.exists(x):
            if os.path.isfile(x):
                os.remove(x)
            elif os.path.isdir(x):
                shutil.rmtree(x, ignore_errors=True)

def remove_tmpdirs(fnames):
    for x in fnames:
        xdir = os.path.dirname(os.path.abspath(x))
        if xdir and os.path.exists(xdir):
            shutil.rmtree(xdir, ignore_errors=True)

def _flatten_plus_safe(rollback_files):
    """
    Flatten names of files and create temporary file names.
    """
    tx_files, orig_files = [], []
    for fnames in rollback_files:
        if isinstance(fnames, six.string_types):
            fnames = [fnames]
        for fname in fnames:
            basedir = safe_makedir(os.path.dirname(fname))
            tmpdir = safe_makedir(tempfile.mkdtemp(dir=basedir))
            tx_file = os.path.join(tmpdir, os.path.basename(fname))
            tx_files.append(tx_file)
            orig_files.append(fname)
    return tx_files, orig_files

def safe_makedir(dname):
    """
    Make a directory if it doesn't exist, handling concurrent race conditions.
    """
    if not dname:
        return dname
    num_tries = 0
    max_tries = 5
    while not os.path.exists(dname):
        try:
            os.makedirs(dname)
        except OSError:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dname


def file_transaction(*rollback_files):
    """
    Wrap file generation in a transaction, moving to output if finishes.
    """
    exts = {".vcf": ".idx", ".bam": ".bai", "vcf.gz": ".tbi", ".fastq.gz": ".count"}

    safe_names, orig_names = _flatten_plus_safe(rollback_files)
    # remove any half-finished transactions
    remove_files(safe_names)
    try:
        if len(safe_names) == 1:
            yield safe_names[0]
        else:
            yield tuple(safe_names)
    # failure -- delete any temporary files
    except:
        remove_files(safe_names)
        remove_tmpdirs(safe_names)
        raise
    # worked -- move the temporary files to permanent location
    else:
        for safe, orig in zip(safe_names, orig_names):
            if os.path.exists(safe):
                shutil.move(safe, orig)
                for check_ext, check_idx in six.iteritems(exts):
                    if safe.endswith(check_ext):
                        safe_idx = safe + check_idx
                        if os.path.exists(safe_idx):
                            shutil.move(safe_idx, orig + check_idx)
        remove_tmpdirs(safe_names)

def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]
        
def read_overlap_pctid(l, pctid, min_len, overlap=0):
    real_len = l.infer_query_length()
    aln_len = l.query_alignment_length
    mismatch = l.get_tag("NM")

    aln_overlap = (aln_len / real_len) * 100
    aln_pctid = ((aln_len - mismatch) / aln_len) * 100
    if aln_overlap >= overlap and aln_pctid >= pctid and aln_len >= min_len:
        return True
    else:
        return False

def filter_bam(bam, outbam, pctid=95, minlen=150, overlap=0,):
    with pysam.AlignmentFile(bam, "rb", check_sq=False) as ih, pysam.AlignmentFile(outbam, "wb", template=ih) as oh:
        good = 0
        good_bp = 0
        outfile = ".".join(outbam.split(".")[:-1]) + ".aln_count"
        
        for i, l in enumerate(ih):
            if l.is_duplicate:
                continue
            
            if read_overlap_pctid(l, pctid, minlen, overlap):
                good += 1
                good_bp += l.query_alignment_length
                oh.write(l)

    with open(outfile, "w") as oh2:
        print(good, good_bp, sep="\n", file=oh2)
    return outbam

def count_bam(bam, pctid=95, minlen=150, overlap=0,):
    good = 0
    good_bp = 0
    
    with pysam.AlignmentFile(bam, "rb", check_sq=False) as ih:
        
        for i, l in enumerate(ih):
            if l.is_duplicate:
                continue
            
            if read_overlap_pctid(l, pctid, minlen, overlap):
                good += 1
                good_bp += l.query_alignment_length
                
    return good, good_bp

def filter_paired_bam(bam, outbam, pctid=95, minlen=150, overlap=0):
    good = 0
    good_bp = 0
    outfile = ".".join(outbam.split(".")[:-1]) + ".aln_count"
    with pysam.AlignmentFile(bamfile, "rb") as ih, pysam.AlignmentFile(outbam, "wb", template=ih) as oh:
        for read1, read2 in read_pair_generator(ih):
            if read1.is_duplicate or read2.is_duplicate:
                continue
                
            if read_overlap_pctid(read1, pctid, minlen, overlap) and read_overlap_pctid(read2, pctid, minlen, overlap):
                good += 1
                good_bp += read1.query_alignment_length + read2.query_alignment_length
                oh.write(read1)
                oh.write(read2)

    with open(outfile, 'w') as oh2:
        print(good, good_bp, sep = "\n", file = oh2)
    return outbam

def count_paired_bam(bam,  pctid=95, minlen=0, overlap=0):
    good = 0
    good_bp = 0
    with pysam.AlignmentFile(bam, "rb") as ih:
        for read1, read2 in read_pair_generator(ih):
            if read1.is_duplicate or read2.is_duplicate:
                continue
                
            if read_overlap_pctid(read1, pctid, minlen, overlap) and read_overlap_pctid(read2, pctid, minlen, overlap):
                good += 1
                good_bp += read1.query_alignment_length + read2.query_alignment_length

    return good, good_bp


def get_coverage(bam_file, bedout=None):
    '''
    create per base coverage patterns from sorted bam
    '''
    bedgraph = ""
    filename, ext = op.splitext(bam_file)
    if bedout is None:
        bedout = filename + ".genomecoverage"

    if op.exists(bedout):
        return bedout

    with file_transaction(bedout) as tx_oh:
        cmd = ("bedtools genomecov -d -ibam {bam} > {out}").format(bam=bam_file, out=tx_oh)
        subprocess.check_call(cmd, shell=True)
    return bedout


def get_recruit_info(gcov):
    '''calculate information on recruited reads based on bedtools genomecoverage table

    Args:
        gcov (str): path to genome coverage file with recruitment pipeline naming convention of:
            metagenome_vs_sag.genomecoverage
            metagenome_vs_sag.aln_count file must also exists within the same directory

    Returns:
        pandas dataframe of genome coverage statistics
    '''
    countfile = gcov.replace("genomecoverage", "aln_count")
    with open(countfile) as infile:
        recruit_count, recruit_bp = infile.read().split("\n")[0:2]
        #recruit_count = infile.read().split()[1].strip()

    metagenome = op.basename(gcov).split("_vs_")[0]
    sag = op.basename(gcov).split("_vs_")[1].split(".")[0]
    cols = ['sag',
            'metagenome',
            'Percent_scaffolds_with_any_coverage',
            'Percent_of_reference_bases_covered',
            'Average_coverage',
            'total_reads_recruited',
            'total_bp_recruited']

    try:
        coverage = pd.read_csv(gcov, sep="\t", header=None)
    except:
        logger.warning("no genome coverage data for {sag}-{metagenome} recruitment".format(sag=sag, metagenome=metagenome))
        data = [sag, metagenome, 0, 0, 0, 0, 0]
        return pd.DataFrame(data, index=cols).transpose()

    mean_per_contig = coverage.groupby([0])[2].mean()
    sum_per_contig = coverage.groupby([0])[2].sum()
    contig_size = coverage.groupby([0])[1].max() + 1
    mean_sag_coverage = mean_per_contig.mean()
    totalbp = contig_size.sum()

    uncovered_bp = len(coverage[coverage[2] == 0])
    pct_covered = ((totalbp - uncovered_bp) / totalbp) * 100
    total_scaffold = len(sum_per_contig)
    uncovered_contig = len(sum_per_contig[sum_per_contig == 0])
    pct_scaffolds_covered = ((total_scaffold - uncovered_contig) / total_scaffold) * 100

    data = [sag,
           metagenome,
           pct_scaffolds_covered,
           pct_covered,
           mean_sag_coverage,
           recruit_count,
           recruit_bp]
    return pd.DataFrame(data, index=cols).transpose()



@cli.command('count-alignments', context_settings=dict(help_option_names=['-h', '--help']), 
             short_help='filter input bam based on pct id, minimum read length and overlap')
@click.argument('inbam')
@click.option('--pctid', 
              default = 95, 
              help='percent identity of read match')
@click.option('--minlen', 
              default = 150, 
              help='minimum read alignment length')
@click.option('--overlap', 
              default = 0, 
              help="minimum overlap")
@click.option('--unpaired',
             is_flag = True,
             help='use flag if input metagenome is unpaired')
def aln_count(inbam, pctid, minlen, overlap, keep_coverage):
    outfile = inbam.replace('.bam','_p-{pctid}_ml-{minlen}_ol-{overlap}.aln_count'.format(pctid=pctid,
                                                                                  minlen=minlen,
                                                                                      overlap=overlap))
    if unpaired:
        count_bam(inbam, outbam, pctid, minlen, overlap)
    else:
        count_paired_bam(inbam, outbam, pctid, minlen, overlap)
                            
    with open(outfile, 'w') as oh:
        print(good, good_bp, sep = "\n", file = oh)
                            
        
@cli.command('count-pctid-range', context_settings=dict(help_option_names=['-h', '--help']), 
             short_help='filter input bam based on pct id, minimum read length and overlap')
@click.argument('inbam')
@click.option('--minlen', 
              default = 0, 
              help='minimum read alignment length')
@click.option('--overlap', 
              default = 0, 
              help="minimum overlap")
def aln_count_pctid_range(inbam, minlen, overlap):
    pid_range = [100, 98, 95, 92, 90, 85, 80, 70]
    #pid_range = [100, 98]
    pids = []
    paired_unpaired = []
    readslist = []
    read_bp_list = []
    outfile = inbam.replace('.bam','_recruit_count_by_pid_range_ml-{minlen}_ol-{overlap}.csv'.format(minlen=minlen,overlap=overlap))
    for pctid in pid_range:
        
        paired_unpaired.append('unpaired')
        pids.append(pctid)
        reads, read_bp = count_bam(inbam, pctid, minlen, overlap)
        readslist.append(reads)
        read_bp_list.append(read_bp)
        
        paired_unpaired.append('paired')
        pids.append(pctid)
        reads, read_bp = count_paired_bam(inbam, pctid, minlen, overlap)
        readslist.append(reads)
        read_bp_list.append(read_bp)
        
    outdf =  pd.DataFrame(data = {'pid':pids, 'paired_unpaired':paired_unpaired, "recruited_read_count": readslist, "bp_recrutied":read_bp_list})
    print(outdf)
    outdf.to_csv(outfile, index = False)
    return outdf
    
@cli.command('filter-and-count', context_settings=dict(help_option_names=['-h', '--help']), 
             short_help='filter input bam based on pct id, minimum read length and overlap')
@click.argument('inbam')
@click.option('--pctid', 
              default = 95, 
              help='percent identity of read match')
@click.option('--minlen', 
              default = 150, 
              help='minimum read alignment length')
@click.option('--overlap', 
              default = 0, 
              help="minimum overlap")
@click.option('--keep-coverage',
             is_flag = True,
             help='delete large genomcoverage and filtered bam file after analysis')
@click.option('--unpaired',
             is_flag = True,
             help='use flag if input metagenome is unpaired')
def filter_and_count(inbam, pctid, minlen, overlap, keep_coverage):
    outbam = inbam.replace('.bam','_p-{pctid}_ml-{minlen}_ol-{overlap}.bam'.format(pctid=pctid,
                                                                                  minlen=minlen,
                                                                                      overlap=overlap))
    if unpaired:
        outbam = filter_bam(inbam, outbam, pctid, minlen, overlap)
    else:
        outbam = filter_paired_bam(inbam, outbam, pctid, minlen, overlap)
    
    gcov = get_coverage(outbam)
    outdf = get_recruit_info(gcov)
    
    outdf_path = outbam.replace('.bam','.csv')
    outdf.to_csv(outdf_path, index = False)
    if not keep_coverage:
        os.remove(gcov)
        os.remove(outbam)
        
if __name__ == '__main__':
    cli()
