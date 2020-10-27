from .common import *


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


@pytest.fixture(scope="module")
def adata_batch():
    adata_ref = sc.datasets.pbmc3k_processed()
    adata = sc.datasets.pbmc68k_reduced()

    var_names = adata_ref.var_names.intersection(adata.var_names)
    adata_ref = adata_ref[:, var_names]
    adata = adata[:, var_names]

    sc.pp.pca(adata_ref)
    sc.pp.neighbors(adata_ref)
    sc.tl.umap(adata_ref)

    sc.tl.ingest(adata, adata_ref, obs='louvain')

    adata_concat = adata_ref.concatenate(adata, batch_categories=['ref', 'new'])
    adata_concat.obs.louvain = adata_concat.obs.louvain.astype('category')
    # fix category ordering
    adata_concat.obs.louvain.cat.reorder_categories(adata_ref.obs.louvain.cat.categories, inplace=True)

    sc.tl.pca(adata_concat)
    adata_concat.obs['celltype'] = adata_concat.obs['louvain']

    return adata_concat


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
