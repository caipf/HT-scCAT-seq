import sys
import matplotlib.pyplot as plt
if "google.colab" in sys.modules:
    print (2)
    #!pip install -q git+https://github.com/theislab/cellrank
import numpy as np
#import cellrank as cr
import scanpy as sc
import scanpy.external as sce
import scvelo as scv
import scipy.stats as st
import cellrank as cr
cr.settings.verbosity = 2
import palantir
scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")

sc.settings.set_figure_params(frameon=False, dpi=120)
#cr.settings.verbosity = 2
import warnings
#
import os
os.chdir("./04_Skin_Fibroblast")
warnings.simplefilter("ignore", category=UserWarning)

#palantir
adata = sc.read("./04_Skin_Fibroblast/root.h5ad")
sce.tl.palantir(adata, n_components=5, knn=30)
sc.write("./04_Skin_Fibroblastpalantir.h5ad",adata)

adata = sc.read("./04_Skin_Fibroblast/palantir.h5ad")

start_cell = 'E13_HU2_BC0503_N04'


pr_res = sce.tl.palantir_results(
    adata,
    early_cell=start_cell,
    ms_data='X_palantir_multiscale',
    num_waypoints=500,
    #terminal_states=terminal_states
)
adata.obs['palantir_pseudotime']=pr_res.pseudotime
sc.pl.embedding(
    adata,
    basis="X_umap",
    color=["palantir_pseudotime", "celltype_l2"],
    color_map="gnuplot2",
    show=False
)
plt.savefig("./04_Skin_Fibroblast/palantir_embedding1.pdf")

pk = cr.kernels.PseudotimeKernel(adata, time_key="palantir_pseudotime")
pk.compute_transition_matrix()
print(pk)
#compute initial and terminal states
g = cr.estimators.GPCCA(pk)
print(g)
g.compute_schur()
g.plot_spectrum(real_only=True)
# compute macrostates
g.compute_macrostates(n_states=5, cluster_key="celltype_l2")
g.plot_macrostates(which="all")
#g.plot_macrostates(which="all", discrete=True, legend_loc="right", s=100,show = False)
plt.savefig("./04_Skin_Fibroblast/macrostates_pb.pdf")
# compute terminal_states
g.predict_terminal_states()
g.plot_macrostates(which="terminal")
plt.savefig("./04_Skin_Fibroblast/terminal_pb.pdf")
#嵌入
g.set_terminal_states(states=["Fib.DC", "Fib.Deep", "Fib.Inter", 'Fib.Lower'])
g.plot_macrostates(which="terminal")
plt.savefig("./04_Skin_Fibroblast/terminal_set_pb.pdf")
# compute initial_states
g.predict_initial_states(n_states=2,allow_overlap=True)
g.plot_macrostates(which="initial",show=False)
plt.savefig("./04_Skin_Fibroblast/initial_pb.pdf")
g.set_initial_states(states=["Fib.Origin"])
print(g)

# write macrostates to AnnData
adata.obs["macrostates"] = g.macrostates
adata.uns["macrostates_colors"] = g.macrostates_memberships.colors
#adata.write("./04_Skin_Fibroblastmacrostates_pb.h5ad")

###driver gene
# compute and visualize fate probabilities
g.compute_fate_probabilities()
g.plot_fate_probabilities(same_plot=False)
plt.savefig("./04_Skin_Fibroblast/fate1_pb.pdf")
g.plot_fate_probabilities(same_plot=True)
plt.savefig("./04_Skin_Fibroblast/fate2_pb.pdf")

#plot
df = g.compute_lineage_drivers()
#print(df)
df.to_csv('driver_gene_pb.csv', index=False)
#beta_drivers = g.compute_lineage_drivers(lineages="Fib.Deep")
beta_drivers = g.compute_lineage_drivers(lineages="Fib.DC")
model = cr.models.GAMR(adata)

# plot heatmap
cr.pl.heatmap(
    adata,
    model=model,  # use the model from before
    lineages="Fib.DC",
    cluster_key="celltype_l2",
    show_fate_probabilities=True,
    genes=beta_drivers.head(40).index,
    time_key="palantir_pseudotime",
    figsize=(12, 10),
    show_all_genes=True,
    weight_threshold=(1e-3, 1e-3),
    save = "./04_Skin_Fibroblast/heatmap_DC_pb.pdf"
    )
# plot cluster_trend
cr.pl.cluster_trends(
    adata,
    model=model,  # use the model from before
    lineage="Fib.DC",
    genes=adata[:, adata.var["highly_variable"]].var_names,
    time_key="palantir_pseudotime",
    weight_threshold=(1e-3, 1e-3),
    n_jobs=8,
    random_state=0,
    clustering_kwargs={"resolution": 0.2, "random_state": 0},
    neighbors_kwargs={"random_state": 0},
    save = "./04_Skin_Fibroblast/cluster_lineage_DC_pb.pdf"
)
del beta_drivers
#Deep
#beta_drivers = g.compute_lineage_drivers(lineages="Fib.Deep")
gene = ["Fbxl7","Ptprd","Gm20631","Efna5","Nebl","Mdga2","Tenm3","Hmcn1","Gm48742","Apod","Mir100hg","Tpm1","Ctnna2","Zfhx3","Tnmd","mt-Nd5","Pde3a","Pcdh9","Pcsk5","Col8a2","Col11a1","Unc5c","Malat1","Zfhx4","Gm29478","Kif26b","Grip1","Cald1","Foxp2","Cacnalc","lgfbp2","Zfpm2","Aff3","Prkg1","Magi2","Angpt1","Bnc2","Epha7","Akap12","Alcam"]
cr.pl.heatmap(
    adata,
    model=model,  # use the model from before
    lineages="Fib.Deep",
    cluster_key="celltype_l2",
    show_fate_probabilities=True,
    genes=gene,
    time_key="palantir_pseudotime",
    figsize=(12, 10),
    show_all_genes=True,
    weight_threshold=(1e-3, 1e-3),
    save = "./04_Skin_Fibroblast/heatmap_Deep_pb.pdf"
    )
    # plot cluster_trend
cr.pl.cluster_trends(
    adata,
    model=model,  # use the model from before
    lineage="Fib.Deep",
    genes=adata[:, adata.var["highly_variable"]].var_names,
    time_key="palantir_pseudotime",
    weight_threshold=(1e-3, 1e-3),
    n_jobs=8,
    random_state=0,
    clustering_kwargs={"resolution": 0.2, "random_state": 0},
    neighbors_kwargs={"random_state": 0},
    save = "./04_Skin_Fibroblast/cluster_lineage_Deep_pb.pdf"
)
#Inter
#
beta_drivers = g.compute_lineage_drivers(lineages="Fib.Inter")
cr.pl.heatmap(
    adata,
    model=model,  # use the model from before
    lineages="Fib.Inter",
    cluster_key="celltype_l2",
    show_fate_probabilities=True,
    genes=beta_drivers.head(40).index,
    time_key="palantir_pseudotime",
    figsize=(12, 10),
    show_all_genes=True,
    weight_threshold=(1e-3, 1e-3),
    save = "./04_Skin_Fibroblast/heatmap_Inter_pb.pdf"
    )
# plot cluster_trend
cr.pl.cluster_trends(
    adata,
    model=model,  # use the model from before
    lineage="Fib.Inter",
    genes=adata[:, adata.var["highly_variable"]].var_names,
    time_key="palantir_pseudotime",
    weight_threshold=(1e-3, 1e-3),
    n_jobs=8,
    random_state=0,
    clustering_kwargs={"resolution": 0.2, "random_state": 0},
    neighbors_kwargs={"random_state": 0},
    save = "./04_Skin_Fibroblast/cluster_lineage_Inter_pb.pdf"
)
del beta_drivers
#lower
#
beta_drivers = g.compute_lineage_drivers(lineages="Fib.Lower")
cr.pl.heatmap(
    adata,
    model=model,  # use the model from before
    lineages="Fib.Lower",
    cluster_key="celltype_l2",
    show_fate_probabilities=True,
    genes=beta_drivers.head(40).index,
    time_key="palantir_pseudotime",
    figsize=(12, 10),
    show_all_genes=True,
    weight_threshold=(1e-3, 1e-3),
    save = "./04_Skin_Fibroblast/heatmap_Lower_pb.pdf"
    )
# plot cluster_trend
cr.pl.cluster_trends(
    adata,
    model=model,  # use the model from before
    lineage="Fib.Lower",
    genes=adata[:, adata.var["highly_variable"]].var_names,
    time_key="palantir_pseudotime",
    weight_threshold=(1e-3, 1e-3),
    n_jobs=8,
    random_state=0,
    clustering_kwargs={"resolution": 0.2, "random_state": 0},
    neighbors_kwargs={"random_state": 0},
    save = "./04_Skin_Fibroblast/cluster_lineage_Lower_pb.pdf"
)
sc.write("./04_Skin_Fibroblastafter_cluster_trend.h5ad",adata)
