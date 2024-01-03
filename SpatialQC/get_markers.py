import pandas as pd
import sqlite3
import os


def get_markers_from_db(species, tissue_class, tissue_type, cancer_type):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    db_path = os.path.join(script_dir, 'data', 'cell_marker.db')
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(
        f"SELECT * FROM markers WHERE species = '{species}' AND tissue_class = '{tissue_class}' And tissue_type = '{tissue_type}' And cancer_type = '{cancer_type}'",
        conn)
    conn.close()

    return df['marker'].tolist()


def marker_decorator(func):
    def wrapper(adata, markers=None, species=None, tissue_class=None, tissue_type=None, cancer_type=None, *args,
                **kwargs):
        if markers is None:
            if species is None or tissue_class is None or tissue_type is None:
                raise ValueError(
                    "If markers are not provided, both species,tissue_class,"
                    "tissue_type and cancer_type must be specified.")
            markers = get_markers_from_db(species, tissue_class, tissue_type, cancer_type)

        return func(adata, markers, *args, **kwargs)

    return wrapper
