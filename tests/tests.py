import matplotlib.figure as fig
import networkx as nx
import pytest

from gedspy import (
    DSA,
    Analysis,
    Enrichment,
    Visualization,
    VisualizationDES,
)


@pytest.fixture
def gene_list():
    return [
        "CACNA1I",
        "CALD1",
        "CAMK1G",
        "CAMK2N1",
        "CAMSAP1",
        "CCL15",
        "CCL16",
        "CCNL2",
        "CCT8P1",
        "CD46",
        "CDC14A",
        "CDK18",
        "CDK19",
        "CES3",
        "CHEK2",
        "CHID1",
        "COL6A3",
        "CPVL",
        "CYP3A43",
        "CYP3A5",
        "DBNL",
        "DUSP10",
        "DUSP9",
        "ECHS1",
        "EGFR",
        "EGR2",
        "ELL2",
        "ERMP1",
        "ESR1",
        "F7",
        "FAM171A1",
        "FAM20C",
        "FGFR2",
        "FH",
        "FLAD1",
        "FUT3",
        "GAA",
        "GBA2",
        "GGCX",
        "GJB1",
        "GLRX5",
        "GNAI2",
        "GNB2",
        "GNB3",
        "GPNMB",
        "GRB10",
        "GRHPR",
        "HMGCS2",
        "HSD17B4",
        "HSP90AB1",
        "IER3IP1",
        "IGF2R",
        "IL1R1"
    ]


@pytest.fixture
def enr(gene_list):
    half = len(gene_list) // 2
    E = Enrichment()
    E.select_features(gene_list[:half])
    E.full_enrichment()
    return E


@pytest.fixture
def analysis(enr):
    results = enr.get_results
    A = Analysis(results)
    A.set_p_value(0.05)
    A.set_test("FISH")
    A.set_correction(None)
    return A


@pytest.fixture
def vis(analysis):
    analysis.full_analysis()
    results = analysis.get_full_results
    return Visualization(results)


@pytest.mark.parametrize(
    "run_method, get_method",
    [
        ("GO_overrepresentation", "get_GO_statistics"),
        ("KEGG_overrepresentation", "get_KEGG_statistics"),
        ("REACTOME_overrepresentation", "get_REACTOME_statistics"),
        ("ViMIC_overrepresentation", "get_ViMIC_statistics"),
        ("DISEASES_overrepresentation", "get_DISEASE_statistics"),
        ("features_specificity", "get_specificity_statistics"),
        ("gene_interaction", "get_features_interactions_statistics"),
    ],
)
def test_analysis_statistics_return_dict(analysis, run_method, get_method):
    getattr(analysis, run_method)()
    result = getattr(analysis, get_method)
    assert isinstance(result, dict)


@pytest.mark.parametrize(
    "run_method, get_method",
    [
        ("REACTOME_network", "get_REACTOME_network"),
        ("KEGG_network", "get_KEGG_network"),
        ("GO_network", "get_GO_network"),
    ],
)
def test_analysis_networks_return_dict(analysis, run_method, get_method):
    getattr(analysis, run_method)()
    result = getattr(analysis, get_method)
    assert isinstance(result, dict)


@pytest.mark.parametrize(
    "plot_method, kwargs",
    [
        (
            "gene_type_plot",
            {"cmap": "summer", "image_width": 6, "image_high": 6, "font_size": 15},
        ),
        (
            "GO_plot",
            {
                "p_val": 0.05,
                "test": "FISH",
                "adj": None,
                "n": 20,
                "side": "left",
                "color": "blue",
                "width": 10,
                "bar_width": 0.5,
                "stat": "p_val",
            },
        ),
        (
            "SPECIFICITY_plot",
            {
                "p_val": 0.05,
                "test": "FISH",
                "adj": None,
                "n": 25,
                "side": "right",
                "color": "bisque",
                "width": 10,
                "bar_width": 0.5,
                "stat": "p_val",
            },
        ),
        (
            "KEGG_plot",
            {
                "p_val": 0.05,
                "test": "FISH",
                "adj": None,
                "n": 25,
                "side": "right",
                "color": "orange",
                "width": 10,
                "bar_width": 0.5,
                "stat": "p_val",
            },
        ),
        (
            "REACTOME_plot",
            {
                "p_val": 0.05,
                "test": "FISH",
                "adj": None,
                "n": 25,
                "side": "right",
                "color": "silver",
                "width": 10,
                "bar_width": 0.5,
                "stat": "p_val",
            },
        ),
        (
            "DISEASES_plot",
            {
                "p_val": 0.05,
                "test": "FISH",
                "adj": None,
                "n": 25,
                "side": "right",
                "color": "thistle",
                "width": 10,
                "bar_width": 0.5,
                "stat": "p_val",
            },
        ),
        (
            "ViMIC_plot",
            {
                "p_val": 0.05,
                "test": "FISH",
                "adj": None,
                "n": 25,
                "side": "right",
                "color": "aquamarine",
                "width": 10,
                "bar_width": 0.5,
                "stat": "p_val",
            },
        ),
        ("blod_markers_plot", {}),
    ],
)
def test_visualization_plots_return_figure(vis, plot_method, kwargs):
    result = getattr(vis, plot_method)(**kwargs)
    assert isinstance(result, fig.Figure)


def test_go_network(vis):
    graph = vis.GOPa_network_create(
        data_set="GO-TERM",
        genes_inc=10,
        gene_int=True,
        genes_only=True,
        min_con=2,
        children_con=True,
        include_childrend=True,
        selected_parents=[],
        selected_genes=[],
    )
    assert isinstance(graph, nx.Graph)


def test_kegg_network(vis):
    graph = vis.GOPa_network_create(
        data_set="KEGG",
        genes_inc=10,
        gene_int=True,
        genes_only=True,
        min_con=2,
        children_con=True,
        include_childrend=True,
        selected_parents=["Nervous system", "Cell motility"],
        selected_genes=[],
    )
    assert isinstance(graph, nx.Graph)


def test_reactome_network(vis):
    graph = vis.GOPa_network_create(
        data_set="REACTOME",
        genes_inc=10,
        gene_int=True,
        genes_only=True,
        min_con=2,
        children_con=False,
        include_childrend=False,
        selected_parents=[],
        selected_genes=["F7", "LAMA5"],
    )
    assert isinstance(graph, nx.Graph)


def test_gi_network(vis):
    graph = vis.GI_network_create()
    assert isinstance(graph, nx.Graph)


def test_auto_ml(vis):
    graph = vis.AUTO_ML_network(
        genes_inc=10,
        gene_int=True,
        genes_only=True,
        min_con=1,
        children_con=True,
        include_childrend=True,
        selected_parents=["Cell motility", "Development and regeneration"],
        selected_genes=["JAK1", "KRAS"],
    )
    assert isinstance(graph, nx.Graph)


def test_gene_scatter(vis):
    result = vis.gene_scatter(
        colors="viridis",
        species="human",
        hclust="complete",
        img_width=None,
        img_high=None,
        label_size=None,
        x_lab="Genes",
        legend_lab="log(TPM + 1)",
    )
    assert isinstance(result, dict)


@pytest.fixture
def set1(gene_list):
    half = len(gene_list) // 2
    return gene_list[:half]


@pytest.fixture
def set2(gene_list):
    half = len(gene_list) // 2
    return gene_list[half:]


@pytest.fixture
def results_set1(set1):
    enr = Enrichment()
    enr.select_features(set1)
    enr.full_enrichment()
    results = enr.get_results
    return results


@pytest.fixture
def results_set2(set2):
    enr = Enrichment()
    enr.select_features(set2)
    enr.full_enrichment()
    results = enr.get_results
    return results


@pytest.fixture
def full_results_set1(results_set1):
    ans = Analysis(results_set1)
    ans.set_p_value(0.05)
    ans.set_test("FISH")
    ans.set_correction(None)
    ans.full_analysis()
    return ans.get_full_results


@pytest.fixture
def full_results_set2(results_set2):
    ans = Analysis(results_set2)
    ans.set_p_value(0.05)
    ans.set_test("FISH")
    ans.set_correction(None)
    ans.full_analysis()
    return ans.get_full_results


@pytest.fixture
def dsa(full_results_set1, full_results_set2):
    d = DSA(full_results_set1, full_results_set2)
    return d


@pytest.fixture
def vis_des(dsa):
    dsa.full_analysis()
    results = dsa.get_results
    return VisualizationDES(results)


@pytest.mark.parametrize(
    "run_method, get_method",
    [
        ("GO_diff", "get_GO_diff"),
        ("KEGG_diff", "get_KEGG_diff"),
        ("REACTOME_diff", "get_REACTOME_diff"),
        ("spec_diff", "get_specificity_diff"),
        ("gi_diff", "get_GI_diff"),
        ("network_diff", "get_networks_diff"),
        ("inter_processes", "get_inter_terms"),
        ("connections_diff", "get_set_to_set_con"),
    ],
)
def test_dsa_returns_dict(dsa, run_method, get_method):
    getattr(dsa, run_method)()
    result = getattr(dsa, get_method)
    assert isinstance(result, dict)


@pytest.mark.parametrize(
    "plot_method, kwargs",
    [
        (
            "diff_SPECIFICITY_plot",
            {
                "p_val": 0.05,
                "test": "FISH",
                "adj": "BH",
                "n": 6,
                "min_terms": 1,
                "selected_set": [],
                "width": 10,
                "bar_width": 0.5,
                "stat": "p_val",
            },
        ),
        (
            "diff_gene_type_plot",
            {
                "set1_name": "Set 1",
                "set2_name": "Set 2",
                "image_width": 12,
                "image_high": 6,
                "font_size": 15,
            },
        ),
        (
            "diff_GO_plot",
            {
                "p_val": 0.05,
                "test": "FISH",
                "adj": "BH",
                "n": 25,
                "min_terms": 5,
                "selected_parent": [],
                "width": 10,
                "bar_width": 0.5,
                "stat": "p_val",
            },
        ),
        (
            "diff_KEGG_plot",
            {
                "p_val": 0.05,
                "test": "FISH",
                "adj": "BH",
                "n": 25,
                "min_terms": 5,
                "selected_parent": [],
                "width": 10,
                "bar_width": 0.5,
                "stat": "p_val",
            },
        ),
        (
            "diff_REACTOME_plot",
            {
                "p_val": 0.05,
                "test": "FISH",
                "adj": "BH",
                "n": 25,
                "min_terms": 5,
                "selected_parent": [],
                "width": 10,
                "bar_width": 0.5,
                "stat": "n",
            },
        ),
    ],
)
def test_vis_des_plots_return_figure(vis_des, plot_method, kwargs):
    plot = getattr(vis_des, plot_method)(**kwargs)
    assert isinstance(plot, fig.Figure)


@pytest.mark.parametrize("dataset", ["GO-TERM", "KEGG", "REACTOME"])
def test_vis_des_diff_GOPa_network(vis_des, dataset):
    graph = vis_des.diff_GOPa_network_create(
        data_set=dataset,
        genes_inc=10,
        gene_int=True,
        genes_only=True,
        min_con=2,
        children_con=False,
        include_childrend=True,
        selected_parents=[],
        selected_genes=[],
    )
    assert isinstance(graph, nx.Graph)


def test_vis_des_gi_network(vis_des):
    graph = vis_des.diff_GI_network_create(min_con=2)
    assert isinstance(graph, nx.Graph)


def test_vis_des_auto_ml_network(vis_des):
    graph = vis_des.diff_AUTO_ML_network(
        genes_inc=10,
        gene_int=True,
        genes_only=True,
        min_con=2,
        children_con=False,
        include_childrend=False,
        selected_parents=[],
        selected_genes=[],
    )
    assert isinstance(graph, nx.Graph)


def test_vis_des_gene_scatter(vis_des):
    result = vis_des.diff_gene_scatter(
        set_num=1,
        colors="viridis",
        species="human",
        hclust="complete",
        img_width=None,
        img_high=None,
        label_size=None,
        x_lab="Genes",
        legend_lab="log(TPM + 1)",
        selected_list=[],
    )
    assert isinstance(result, dict)
