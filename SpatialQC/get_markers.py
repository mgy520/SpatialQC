import pandas as pd
import sqlite3
import os


script_dir = os.path.dirname(os.path.realpath(__file__))
db_path = os.path.join(script_dir, 'data', 'cell_marker.db')


def get_available_options_from_db(option_type):
    conn = sqlite3.connect(db_path)
    options = pd.read_sql_query(f"SELECT DISTINCT {option_type} FROM markers", conn)[option_type].tolist()
    conn.close()
    return options


def get_markers_from_db(species, tissue_class, tissue_type, cancer_type):
    available_species = get_available_options_from_db('species')
    available_tissue_class = get_available_options_from_db('tissue_class')
    available_tissue_type = get_available_options_from_db('tissue_type')
    available_cancer_type = get_available_options_from_db('cancer_type')
    if species not in available_species:
        raise ValueError(f"Invalid species: '{species}'. Available species are: {', '.join(available_species)}")
    if tissue_class not in available_tissue_class:
        raise ValueError(
            f"Invalid tissue class: '{tissue_class}'. Available tissue class are: {', '.join(available_tissue_class)}")
    if tissue_type not in available_tissue_type:
        raise ValueError(
            f"Invalid tissue type: '{tissue_type}'. Available tissue type are: {', '.join(available_tissue_type)}")
    if cancer_type not in available_cancer_type:
        raise ValueError(
            f"Invalid cancer type: '{cancer_type}'. Available cancer type are: {', '.join(available_cancer_type)}")

    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(
        f"SELECT * FROM markers WHERE species = '{species}' AND tissue_class = '{tissue_class}' And tissue_type = '{tissue_type}' And cancer_type = '{cancer_type}'",
        conn)
    conn.close()

    if df.empty:
        raise ValueError(
            f"No marker genes found for the specified parameters: species='{species}', "
            f"tissue_class='{tissue_class}', tissue_type='{tissue_type}', cancer_type='{cancer_type}'. "
            f"Please provide valid parameters.")

    return df['marker'].tolist()


def marker_decorator(func):
    def wrapper(adata, markers=None, species=None, tissue_class=None, tissue_type=None, cancer_type=None, *args,
                **kwargs):
        if markers is None:
            if species is None or tissue_class is None or tissue_type is None:
                raise ValueError(
                    "If markers are not provided, both species,tissue_class,"
                    "tissue_type must be specified.")
            markers = get_markers_from_db(species, tissue_class, tissue_type, cancer_type)

        return func(adata, markers, *args, **kwargs)

    return wrapper
