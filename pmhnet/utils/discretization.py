import numpy as np
from pmhnet.utils.dataset_creation import execute_and_fetch

def percentile_breaks(n_intervals, max_time):
    """Percentile breaks based on max_time and n_intervals"""

    query = """
        SELECT time FROM surv_master
        WHERE time <= %s AND event;
    """
    times = np.concatenate(execute_and_fetch(query, (max_time, )))
    times = np.clip(times, None, max_time)
    percentile_breaks = np.linspace(0, 100, n_intervals+1)

    return np.percentile(times, percentile_breaks)


def compute_breaks(cursor, max_time, n_intervals):
    """Percentile breaks based on observed event times on the training data.

    Arguments:
    * cursor: sql connection cursor
    * max_time: cutoff of event times
    * n_intervals: number of intervals
    """
    cursor.execute(
        """
        SELECT time FROM surv_master
        WHERE time <= %s AND event AND train
        """,
        (max_time, )
    )
    event_times = np.concatenate(cursor.fetchall())
    pbreaks = np.linspace(0, 100, n_intervals+1)
    return np.percentile(event_times, pbreaks)
