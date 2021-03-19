from .common import *


@pytest.fixture(scope="session")
def adata_paul15_template():
    adata = sc.datasets.paul15()
    adata.obs['celltype'] = adata.obs['paul15_clusters']
    np.random.seed(42)
    adata.obs['batch'] = np.random.randint(1, 5, adata.n_obs)
    adata.obs['batch'] = adata.obs['batch'].astype(str)
    adata.obs['batch'] = adata.obs['batch'].astype("category")
    adata.layers['counts'] = adata.X
    scIB.preprocessing.reduce_data(
        adata,
        pca=False,
        n_top_genes=None,
        neighbors=False,
        umap=False
    )
    yield adata
    del adata


@pytest.fixture(scope="session")
def adata_pbmc_template():
    # adata_ref = sc.datasets.pbmc3k_processed()
    # quick fix for broken dataset paths, should be removed with scanpy>=1.6.0
    adata_ref = sc.read(
        "pbmc3k_processed.h5ad",
        backup_url="https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad"
    )
    adata = sc.datasets.pbmc68k_reduced()

    var_names = adata_ref.var_names.intersection(adata.var_names)
    adata_ref = adata_ref[:, var_names]
    adata = adata[:, var_names]

    sc.pp.pca(adata_ref)
    sc.pp.neighbors(adata_ref)
    sc.tl.umap(adata_ref)

    # merge cell type labels
    sc.tl.ingest(adata, adata_ref, obs='louvain')
    adata_concat = adata_ref.concatenate(adata, batch_categories=['ref', 'new'])
    adata_concat.obs.louvain = adata_concat.obs.louvain.astype('category')
    # fix category ordering
    adata_concat.obs.louvain.cat.reorder_categories(adata_ref.obs.louvain.cat.categories, inplace=True)
    adata_concat.obs['celltype'] = adata_concat.obs['louvain']

    del adata_concat.obs['louvain']
    del adata_concat.uns
    del adata_concat.obsm
    del adata_concat.varm

    yield adata_concat
    del adata_concat


@pytest.fixture()
def adata_paul15(adata_paul15_template):
    adata_obj = adata_paul15_template.copy()
    yield adata_obj
    del adata_obj


@pytest.fixture()
def adata(adata_pbmc_template):
    adata_obj = adata_pbmc_template.copy()
    yield adata_obj
    del adata_obj


@pytest.fixture()
def adata_pca(adata):
    adata_obj = adata
    scIB.pp.reduce_data(
        adata_obj,
        pca=True,
        n_top_genes=200,
        neighbors=False,
        umap=False
    )
    yield adata_obj


@pytest.fixture()
def adata_neighbors(adata):
    adata_obj = adata
    scIB.pp.reduce_data(
        adata_obj,
        pca=True,
        n_top_genes=200,
        neighbors=True,
        umap=False
    )
    yield adata_obj


@pytest.fixture()
def adata_clustered(adata_neighbors):
    adata_obj = adata_neighbors
    scIB.cl.opt_louvain(
        adata_obj,
        cluster_key='cluster',
        label_key='celltype',
        verbose=True
    )
    yield adata_obj
