'''
Generate random identifiers for anonymisation.

Previously used random identifiers are stored in an SQLite database.
'''

import random
import sys
import sqlite3
from error import print_error, ERROR_RANDOM_ID_ITERATIONS

MAX_RANDOM_ID_ITERATIONS = 1000000
DEFAULT_USED_IDS_DATABASE = "used_random_sample_ids.db"
# Start the minimum random sample ID at some non-low number.
# Users might think there is something special about low number IDs
MIN_RANDOM_ID = 1000
# XXX This may not be portable: we might want to pick a high
# maximum value ourselves
#MAX_RANDOM_ID = sys.maxint
MAX_RANDOM_ID = sys.maxsize

def make_one_random_id():
    return random.randint(MIN_RANDOM_ID, MAX_RANDOM_ID)

# We keep an sqlite database of previously used random IDs. This ensures that
# we never re-use the same random ID. It may seem like overkill to use a
# database, but the benefit is that sqlite will handle the locking of
# database access, so that we do not have race conditions if multiple instances
# of this program run at the same time.
# If the database does not exist we will create a new empty one.
def make_random_ids(used_ids_database, sample_ids):
    conn = sqlite3.connect(used_ids_database)
    cursor = conn.cursor()
    cursor.execute('CREATE TABLE IF NOT EXISTS unique_ids (id integer)')
    # The set of IDs that have already been used
    used_ids = set([])
    # Grab all the IDs from the database
    for (next_used_id,) in cursor.execute('SELECT * from unique_ids'):
        used_ids.add(next_used_id)
    # The newly generated IDs by the call to this function
    new_ids = []
    # A list of pairs containing the original ID and its new randomised ID
    result = {}
    for old_sample in sample_ids:
        # Make sure the newly generated ID has not been seen before
        new_id = make_one_random_id()
        iter_count = 0
        while new_id in used_ids:
            if iter_count >= MAX_RANDOM_ID_ITERATIONS:
                print_error("Could not make a new random ID, iteration count exceeded")
                exit(ERROR_RANDOM_ID_ITERATIONS)
            new_id = make_one_random_id()
            iter_count += 1
        result[old_sample] = new_id
        # Record this new ID in the set of previously used IDs so we don't
        # use it again
        used_ids.add(new_id)
        new_ids.append(new_id)
    # Write the newly created IDs out to the database
    # XXX should be able to do this as a single INSERT statement
    for new_id in new_ids:
        cursor.execute('INSERT into unique_ids (id) VALUES ({})'.format(new_id))
    conn.commit()
    return result
