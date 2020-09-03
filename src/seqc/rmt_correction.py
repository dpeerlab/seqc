import os
from seqc.sequence.encodings import DNA3Bit
import numpy as np
from seqc import log
from seqc.read_array import ReadArray
import time
import pandas as pd
import multiprocessing as multi
from tqdm import tqdm
import dask
from distributed import Client, LocalCluster
from dask.distributed import wait, performance_report
from tlz import partition_all

log.logging.getLogger("asyncio").setLevel(log.logging.WARNING)

# todo document me
def generate_close_seq(seq):
    """ Return a list of all sequences that are up to 2 hamm distance from seq
    :param seq:
    """
    res = []
    l = DNA3Bit.seq_len(seq)

    # generate all sequences that are dist 1
    for i in range(l):
        mask = 0b111 << (i * 3)
        cur_chr = (seq & mask) >> (i * 3)
        res += [
            seq & (~mask) | (new_chr << (i * 3))
            for new_chr in DNA3Bit.bin2strdict.keys()
            if new_chr != cur_chr
        ]
    # generate all sequences that are dist 2
    for i in range(l):
        mask_i = 0b111 << (i * 3)
        chr_i = (seq & mask_i) >> (i * 3)
        for j in range(i + 1, l):
            mask_j = 0b111 << (j * 3)
            chr_j = (seq & mask_j) >> (j * 3)
            mask = mask_i | mask_j
            res += [
                seq & (~mask) | (new_chr_i << (i * 3)) | (new_chr_j << (j * 3))
                for new_chr_i in DNA3Bit.bin2strdict.keys()
                if new_chr_i != chr_i
                for new_chr_j in DNA3Bit.bin2strdict.keys()
                if new_chr_j != chr_j
            ]

    return res


# todo document me
def probability_for_convert_d_to_r(d_seq, r_seq, err_rate):
    """
    Return the probability of d_seq turning into r_seq based on the err_rate table
    (all binary)

    :param err_rate:
    :param r_seq:
    :param d_seq:
    """

    if DNA3Bit.seq_len(d_seq) != DNA3Bit.seq_len(r_seq):
        return 1

    p = 1.0
    while d_seq > 0:
        if d_seq & 0b111 != r_seq & 0b111:
            if isinstance(err_rate, float):
                p *= err_rate
            else:
                p *= err_rate[(d_seq & 0b111, r_seq & 0b111)]
        d_seq >>= 3
        r_seq >>= 3
    return p


def in_drop(read_array, error_rate, alpha=0.05):
    """ Tag any RMT errors

    :param read_array: Read array
    :param error_rate: Sequencing error rate determined during barcode correction
    :param alpha: Tolerance for errors
    """

    return _correct_errors(read_array, error_rate, alpha)


# a method called by each process to correct RMT for each cell
def _correct_errors_by_cell_group(ra, cell_group, err_rate, p_value):

    from scipy.special import gammainc

    # cell_group = indices_grouped_by_cells[cell_index]
    # Breaks for each gene
    gene_inds = cell_group[np.argsort(ra.genes[cell_group])]
    breaks = np.where(np.diff(ra.genes[gene_inds]))[0] + 1
    splits = np.split(gene_inds, breaks)
    rmt_groups = {}
    res = []

    for inds in splits:
        # RMT groups
        for ind in inds:
            rmt = ra.data["rmt"][ind]
            try:
                rmt_groups[rmt].append(ind)
            except KeyError:
                rmt_groups[rmt] = [ind]

        if len(rmt_groups) == 1:
            continue

        # This logic retains RMTs with N if no donor is found and contributes to the
        # molecule count
        for rmt in rmt_groups.keys():

            # Enumerate all possible RMTs with hamming distances 1 and/or 2
            # to build a probablitiy that this particular RMT was not an error
            # Simulatenously, check if Jaitin error correction can be applied
            jaitin_corrected = False
            expected_errors = 0
            for donor_rmt in generate_close_seq(rmt):

                # Check if donor is detected
                try:
                    donor_count = len(rmt_groups[donor_rmt])
                except KeyError:
                    continue

                # Build likelihood
                # Probability of converting donor to target
                p_dtr = probability_for_convert_d_to_r(donor_rmt, rmt, err_rate)
                # Number of occurrences
                expected_errors += donor_count * p_dtr

                # Check if jaitin correction is feasible
                if not jaitin_corrected:
                    ref_positions = ra.positions[rmt_groups[rmt]]
                    donor_positions = ra.positions[rmt_groups[donor_rmt]]

                    # Is reference a subset of the donor ?
                    if (set(ref_positions)).issubset(donor_positions):
                        jaitin_corrected = True
                        jaitin_donor = donor_rmt

            # Probability that the RMT is an error
            p_val_err = gammainc(len(rmt_groups[rmt]), expected_errors)

            # Remove Jaitin corrected reads if probability of RMT == error is high
            if p_val_err > p_value and jaitin_corrected:
                # Save the RMT donor
                # save the index of the read and index of donor rmt read
                for i in rmt_groups[rmt]:
                    res.append((i, rmt_groups[jaitin_donor][0]))

        rmt_groups.clear()

    return res


def _correct_errors_by_cell_group_chunks(ra, cell_group_chunks, err_rate, p_value):
    return [
        _correct_errors_by_cell_group(ra, cell_group, err_rate, p_value)
        for cell_group in cell_group_chunks
    ]


def _correct_errors(ra, err_rate, p_value=0.05):

    n_workers = int(os.environ.get("SEQC_MAX_NUM_CPU", max(multi.cpu_count() - 1, 1)))

    # configure dask.distributed
    # use total number of available CPUs - 1
    # memory_terminate_fraction doesn't work for some reason
    # https://github.com/dask/distributed/issues/3519
    # https://docs.dask.org/en/latest/setup/single-distributed.html#localcluster
    # https://docs.dask.org/en/latest/scheduling.html#local-threads
    worker_kwargs = {
        "n_workers": n_workers,
        "threads_per_worker": 1,
        "processes": True,
        "memory_limit": "64G",
        "memory_target_fraction": 0.95,
        "memory_spill_fraction": 0.99,
        "memory_pause_fraction": False,
        # "memory_terminate_fraction": False,
    }

    # do not kill worker at 95% memory level
    dask.config.set({"distributed.worker.memory.terminate": False})
    dask.config.set({"distributed.scheduler.allowed-failures": 50})

    # setup Dask distributed client
    cluster = LocalCluster(**worker_kwargs)
    client = Client(cluster)

    # debug message
    log.debug(
        "Dask processes={} threads={}".format(
            len(client.nthreads().values()), np.sum(list(client.nthreads().values()))
        ),
        module_name="rmt_correction",
    )
    log.debug(
        "Dask worker_kwargs "
        + " ".join([f"{k}={v}" for k, v in worker_kwargs.items()]),
        module_name="rmt_correction",
    )
    log.debug("Dask Dashboard=" + client.dashboard_link, module_name="rmt_correction")

    # group by cells (same cell barcodes as one group)
    log.debug("Grouping...", module_name="rmt_correction")
    indices_grouped_by_cells = ra.group_indices_by_cell()

    # send readarray in advance to all workers (i.e. broadcast=True)
    # this way, we reduce the serialization time
    log.debug("Scattering ReadArray...", module_name="rmt_correction")
    [future_ra] = client.scatter([ra], broadcast=True)

    # correct errors per cell group in parallel
    log.debug("Submitting jobs to Dask...", module_name="rmt_correction")
    with performance_report(filename="dask-report.html"):
        futures = []

        # 50 chunks at a time
        chunks = partition_all(50, indices_grouped_by_cells)

        for chunk in tqdm(chunks):
            future = client.submit(
                _correct_errors_by_cell_group_chunks,
                future_ra,
                chunk,
                err_rate,
                p_value,
            )
            futures.append(future)

        # wait until all done
        log.debug("Waiting untill all tasks complete...", module_name="rmt_correction")
        completed, not_completed = wait(futures)

    if len(not_completed) > 1:
        raise Exception("There are uncompleted tasks!")

    # gather the resutls and release
    log.debug(
        "Collecting the task results from the workers...", module_name="rmt_correction"
    )
    results = []
    for future in tqdm(completed):
        # this returns a list of a list
        # len(result) should be the number of chunks e.g. 50
        result = future.result()

        # remove empty lists
        result = list(filter(lambda x: len(x) > 0, result))

        # aggregate and release
        results.extend(result)
        future.release()

    # clean up
    del futures
    del completed
    del not_completed

    client.shutdown()
    client.close()

    # iterate through the list of returned read indices and donor rmts
    # create a mapping tble of pre-/post-correction
    mapping = []
    for result in results:
        for idx, idx_corrected_rmt in result:

            # record pre-/post-correction
            mapping.append(
                (
                    ra.data["cell"][idx],
                    ra.data["rmt"][idx],
                    ra.data["rmt"][idx_corrected_rmt],
                )
            )

    # iterate through the list of returned read indices and donor rmts
    # actually, update the read array object with corrected UMI
    for result in results:
        for idx, idx_corrected_rmt in result:

            # correct
            ra.data["rmt"][idx] = ra.data["rmt"][idx_corrected_rmt]

            # report error
            ra.data["status"][idx] |= ra.filter_codes["rmt_error"]

    return pd.DataFrame(mapping, columns=["CB", "UR", "UB"])
