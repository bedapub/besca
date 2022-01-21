from dominate import document
from dominate.tags import *
from dominate.util import text, raw
from datetime import date
from os import remove
from os.path import join, dirname
import besca as bc
from shutil import copyfile

# from weasyprint import HTML, CSS

from ..pp._fraction_pos import top_counts_genes
from .._helper import get_raw


def write_qc(
    adata_unfiltered,
    adata_filtered,
    version,
    analysis_name,
    standard_min_genes,
    standard_min_cells,
    standard_min_counts,
    standard_percent_mito,
    standard_max_counts,
    standard_n_genes,
    filtering_output1,
    filtering_output2,
    results_folder,
    css_path,
):

    """generates the first part of the documentation"""

    adata = adata_unfiltered
    copyfile(css_path, join(results_folder, "style.css"))

    # part one before filtering
    doc = document(title="QC Report - " + analysis_name)

    with doc.head:
        link(rel="stylesheet", href="style.css")

    with doc:

        with div():
            attr(cls="body")

            h1("Quality Control Report - Single Cell Analysis Standard Pipeline")

            with p(style="font-size:12px"):
                text("single cell analysis standard workflow version " + str(version))
                br()
                text("analysis name: " + str(analysis_name))
                br()
                text("report generated on: " + str(date.today()))

            h2("Raw Dataset Properties")

            # generate a table
            with table():
                attr(cls="minimalistBlack")
                l = tr()
                l += th("total number of cells")
                l += td(adata.n_obs)

                l = tr()
                l += th("number of genes")
                l += td(adata.n_vars)

                if "donor" in adata.obs.columns.tolist():
                    l = tr()
                    l += th("number of donors")
                    l += td(len(adata.obs.donor.value_counts().index.tolist()))

                if "condition" in adata.obs.columns.tolist():
                    l = tr()
                    l += th("number of conditions")
                    l += td(len(adata.obs.condition.value_counts().index.tolist()))

                if "treatment" in adata.obs.columns.tolist():
                    l = tr()
                    l += th("treatments")
                    l += td(adata.obs.donor.value_counts().index.tolist())

            br()
            img(src="./figures/transcriptcaptureefficiency.png", width="500px")
            br()
            img(src="./figures/librarysize.png", width="500px")
            br()

            h2("Filtering")

            with div():
                attr(style="float: left")
                with table():
                    attr(cls="minimalistBlack2", style="float: left")
                    with thead():
                        l = tr()
                        l += th("used filtering parameters", colspan=2)

                    with tbody():
                        l = tr()
                        l += th("min. genes per cell: ")
                        l += td(standard_min_genes)

                        l = tr()
                        l += th("min. UMI counts per cell: ")
                        l += td(standard_min_counts)

                        l = tr()
                        l += th("max. genes per cell: ")
                        l += td(standard_n_genes)

                        l = tr()
                        l += th("max. UMI counts per cell: ")
                        l += td(standard_max_counts)

                        l = tr()
                        l += th("max. mitochondrial gene content: ")
                        l += td(standard_percent_mito)

                        l = tr()
                        l += th("min. cells expressing a gene: ")
                        l += td(standard_min_cells)

                    with tfoot():
                        l = tr()
                        l += th("", colspan=2)
                with div():
                    attr(style="float: left")

                    with p(style="font-size:10px"):
                        for x in filtering_output1.stdout.split("\n"):
                            text(x)
                            br()
                        for x in filtering_output2.stdout.split("\n"):
                            text(x)
                            br()

            with div(style="clear: both"):

                h4("visualization of filtering thresholds")
                img(src="./figures/filtering_thresholds.png", width="500px")

                br()

                with table():
                    attr(cls="minimalistBlack", style="float: none")

                    with thead():
                        l = tr()
                        l += th("")
                        l += th("number of cells before filtering")
                        l += th("number of cells after filtering")

                    with tbody():
                        l = tr()
                        l += th("total cells")
                        l += td(adata.n_obs)

                        if "donor" in adata.obs.columns.tolist():
                            table_donor = {}
                            for donor in adata.obs.donor.value_counts().index.tolist():
                                subset_cells = (
                                    adata[adata.obs.donor == donor, :].copy().n_obs
                                )
                                table_donor[donor] = tr()
                                table_donor[donor] += th(donor)
                                table_donor[donor] += td(subset_cells)

                        if "condition" in adata.obs.columns.tolist():
                            table_condition = {}
                            for (
                                condition
                            ) in adata.obs.condition.value_counts().index.tolist():
                                subset_cells = (
                                    adata[adata.obs.condition == condition, :]
                                    .copy()
                                    .n_obs
                                )
                                table_condition[condition] = tr()
                                table_condition[condition] += th(condition)
                                table_condition[condition] += td(subset_cells)

                        if "treatment" in adata.obs.columns.tolist():
                            table_treatment = {}
                            for (
                                treatment
                            ) in adata.obs.treatment.value_counts().index.tolist():
                                subset_cells = (
                                    adata[adata.obs.treatment == treatment, :]
                                    .copy()
                                    .n_obs
                                )
                                table_treatment[treatment] = tr()
                                table_treatment[treatment] += th(treatment)
                                table_treatment[treatment] += td(subset_cells)
                    with tfoot():
                        foot = tr()
                        foot += th("")
                        foot += th("")
                        foot += th("")

                    br()
    # part two after filtering
    adata = adata_filtered

    # add number of cells
    l += td(adata.n_obs)

    if "donor" in adata.obs.columns.tolist():
        for donor in adata.obs.donor.value_counts().index.tolist():
            subset_cells = adata[adata.obs.donor == donor, :].copy().n_obs
            table_donor[donor] += td(subset_cells)

    if "condition" in adata.obs.columns.tolist():
        for condition in adata.obs.condition.value_counts().index.tolist():
            subset_cells = adata[adata.obs.condition == condition, :].copy().n_obs
            table_condition[condition] += td(subset_cells)

    if "treatment" in adata.obs.columns.tolist():
        for treatment in adata.obs.treatment.value_counts().index.tolist():
            subset_cells = adata[adata.obs.treatment == treatment, :].copy().n_obs
            table_treatment[treatment] += td(subset_cells)

    # create table overview of the used filtering parameters
    with doc:

        with div(style="clear: both"):
            h2("Filtered Dataset Properties")

            img(src="./figures/violin.after_filtering.png", width="500px")

            br()

            img(src="./figures/top_genes.png", width="500px")

            topn = top_counts_genes(adata=get_raw(adata), top_n=10)

            raw(topn.to_html(classes="", header=True, index=False))

            h2("Highly Variable Gene Selection")

            img(src="./figures/filter_genes_dispersion.hvg.png", width="500px")

            h2("Principle Component Analysis")
            br()
            img(src="./figures/PCA.png", width="500px")
            br()

            h2("Clustering")
            colnames = list(adata.obs.columns)

            # get clustering algotithm to plot, if both computed, choose leiden
            clust_alg = sorted(
                [item for item in colnames if item in ["leiden", "louvain"]]
            )[0]
            img(src="./figures/umap." + str(clust_alg) + ".png", width="300px")

    with open(join(results_folder, "qc_report.html"), "w") as f:
        f.write(doc.render())

    # this commented out code is functional but was removed from besca because we did not manage to include weasybuild in the easybuild module

    # convert to a pdf
    # HTML(join(results_folder, 'qc_report.html')).write_pdf(join(results_folder, 'qc_report.pdf'), presentational_hints=True)

    # remove copied css.style
    # remove(join(results_folder, 'style.css'))
    # remove(join(results_folder, 'qc_report.html'))

    return None
