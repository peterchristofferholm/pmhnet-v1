import numpy as np
import psycopg2
from psycopg2.extras import NamedTupleCursor


def make_survarray(time, event, breaks):
    """Creates label transformation.

    Args:
        time (integer): Time of event
        event (integer): Failure = 1, Censoring = 0.
        breaks (np.array): Discretization grid. Includes both 0 and max_time.
    """

    surv = np.zeros((2, len(breaks)-1), dtype=np.bool)

    # with failure event
    if event == 1 and time < breaks[-1]:
        surv[1, np.argmax(time <= breaks[1:])] = 1
        surv[0, :] = (time >= breaks[1:])

    # or if censored
    else:
        mids = breaks[:-1] + np.diff(breaks)/2
        surv[0, :] = (time >= mids)  # give credit to half

    return surv


def execute_and_fetch(query, params=()):
    """Executes a postgreSQL query and returns the results

    Args:
        query (string): sql query
        params (tuple): parameters being passed to the query
    """

    # TODO: Replace hardcoded vars with something better
    DBHOST = "dbhost"
    DBNAME = "dbname"
    DBUSER = "dbuser"
    SCHEMA = "pmhnet"

    connection = psycopg2.connect(
        f"host={DBHOST} dbname={DBNAME} user={DBUSER}",
        options = f"-c search_path={SCHEMA}",
        cursor_factory = NamedTupleCursor
    )

    with connection as con:

        cur = con.cursor()
        cur.execute(query, params)

        results = cur.fetchall()

    connection.close()
    return results
