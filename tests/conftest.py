from .common import *


@pytest.fixture(scope="session")
def adata_pbmc():
    #adata_ref = sc.datasets.pbmc3k_processed()
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


@pytest.fixture(scope="function")
def adata(adata_pbmc):
    adata = adata_pbmc.copy()
    yield adata
    del adata

@pytest.fixture(scope="session")
def adata_pca():
    def adata_pca_(adata):
        scIB.pp.reduce_data(adata, pca=True, n_top_genes=200, neighbors=False, umap=False)
        return adata
    return adata_pca_


@pytest.fixture(scope="session")
def adata_neighbors():
    def adata_neighbors_(adata):
        scIB.pp.reduce_data(adata, pca=True, n_top_genes=200, neighbors=True, umap=False)
        return adata
    return adata_neighbors_


@pytest.fixture(scope="session")
def adata_clustered():
    def adata_clustered_(adata):
        scIB.pp.reduce_data(adata, pca=True, n_top_genes=200, neighbors=True, umap=False)
        scIB.cl.opt_louvain(adata, cluster_key='cluster', label_key='celltype', verbose=True)
        return adata
    return adata_clustered_


@pytest.fixture(scope="session")
def adata_factory():
    def adata_factory_(pca=False, n_top_genes=None, neighbors=False):
        adata = sc.datasets.paul15()
        adata.obs['celltype'] = adata.obs['paul15_clusters']
        np.random.seed(42)
        adata.obs['batch'] = np.random.randint(1, 5, adata.n_obs)
        adata.obs['batch'] = adata.obs['batch'].astype(str)
        adata.obs['batch'] = adata.obs['batch'].astype("category")
        adata.layers['counts'] = adata.X
        scIB.preprocessing.reduce_data(
            adata,
            pca=pca,
            n_top_genes=n_top_genes,
            umap=False,
            neighbors=neighbors
        )
        return adata
    return adata_factory_


@pytest.fixture(scope="session")
def embed_factory():
    def adata_embed_(adata, type_):
        if type_ == 'pca':
            if 'X_pca' in adata.obsm:
                mtx = adata.obsm['X_pca']
            else:
                mtx = sc.tl.pca(adata, copy=True).obsm['X_pca']
        elif type_ == 'full':
            mtx = adata.X
        else:
            raise ValueError(f"'{type_}' not a valid embedding type")
        adata.obsm['X_emb'] = mtx
        return adata
    return adata_embed_


@pytest.fixture(scope="session")
def cluster_factory():
    def cluster_factory_(adata, label_key, cluster_key="cluster", verbose=False):
        res_max, score_max, score_all = scIB.cl.opt_louvain(
            adata,
            label_key=label_key,
            cluster_key=cluster_key,
            plot=False, inplace=True, force=True,
            verbose=verbose
        )
        return res_max, score_max, score_all
    return cluster_factory_

