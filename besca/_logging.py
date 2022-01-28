import logging
import os
import datetime


def initialize_logger(log_file):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # create standard log file handler and set level to info
    handler = logging.FileHandler(log_file, "a", encoding=None, delay="true")
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter("%(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # create console handler and set level to info
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter("LOG MESSAGE: %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def initialize_log_file(
    analysis_name,
    root_path,
    species,
    batch_to_correct,
    standard_min_genes,
    standard_min_cells,
    standard_min_counts,
    standard_n_genes,
    standard_percent_mito,
    standard_max_counts,
    version,
):

    logging.info("Standard Pipeline Version " + version + " used")
    logging.info(datetime.datetime.today().strftime("%Y-%m-%d"))
    logging.info(
        "Analysis '" + analysis_name + "' on data located in'" + root_path + "'"
    )
    logging.info("species: " + species)
    logging.info("Batch effect to correct: " + batch_to_correct)
    logging.info("Parameters:")
    logging.info("\tstandard_min_genes = " + str(standard_min_genes))
    logging.info("\tstandard_min_cells = " + str(standard_min_cells))
    logging.info("\tstandard_min_counts = " + str(standard_min_counts))
    logging.info("\tstandard_n_genes = " + str(standard_n_genes))
    logging.info("\tstandard_max_counts = " + str(standard_max_counts))
    logging.info("\tstandard_percent_mito = " + str(standard_percent_mito))
