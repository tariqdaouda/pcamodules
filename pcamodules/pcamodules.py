import numpy as np
from sklearn.linear_model import LogisticRegressionCV

class NeatObject(object):
    """docstring for NeatObject"""
    def __init__(self, **kwargs):
        super(NeatObject, self).__init__()
        for kw, val in kwargs.items():
            setattr(self, kw, val)
        
        self.keys = kwargs.keys()

    def __repr__(self):
        return "<neat object: available attributes: %s>" % self.keys

class PCAModules(object):
    """docstring for PCAModules"""
    def __init__(self, cv=2, random_state=0, max_iter=100, class_weight="balanced", **logreg_kwargs):
        super(PCAModules, self).__init__()
        self.clf = LogisticRegressionCV(cv=cv, random_state=random_state, class_weight=class_weight, max_iter=max_iter, **logreg_kwargs)
        self.train_acc = None
        self.test_acc = None
        self.coefs = None

        self.genes = {}
        self.ordered_pcs = {}
        self.weights = {}

    def fit(self, adata, obs_key="leiden", max_normalize=True, max_pcs=None):
        """
        Fit the model to find the importance of every PC for each clister and 
        weights genes accordingly
        """
        stop = int(len(adata.X) * 0.2)
        X_test = adata.obsm["X_pca"][:stop]
        y_test = adata.obs[obs_key].values[:stop].astype('int32')

        X_train = adata.obsm["X_pca"][stop:]
        y_train = adata.obs[obs_key].values[stop:].astype('int32')

        self.clf.fit(X_train, y_train)
        
        self.train_acc = self.clf.score(X_train, y_train)
        self.test_acc = self.clf.score(X_test, y_test)
        print("Train acc", self.train_acc)
        print("Test acc", self.test_acc)

        self.coefs = self.clf.coef_
    
        for aidi in np.unique(y_test):
            self.ordered_pcs[aidi] = np.argsort(self.coefs[aidi])[::-1]
            pcs = None
            if max_pcs is not None :
                pcs = self.ordered_pcs[aidi][:max_pcs]
            
            genes, weights = self._weight_genes(adata, self.coefs[aidi], max_normalize=max_normalize, pcs=pcs)
            self.genes[aidi] = genes
            self.weights[aidi] = weights

        return self
    
    def select(self, pattern="RP", exclude_genes=True, max_threshold=1, min_threshold=0.5):
        """Select only a subet of genes, returns a NeatObjet"""
        from collections import defaultdict
        res_genes = defaultdict(list)
        res_weights = defaultdict(list)
        for aidi in self.genes.keys():
            for gene, weight in zip(self.genes[aidi], self.weights[aidi]):
                if self._select_gene(gene, weight, pattern=pattern, exclude_gene=exclude_genes) \
                    and self._filter_gene(gene, weight, min_threshold=min_threshold):
                    
                    res_genes[aidi].append(gene)
                    res_weights[aidi].append(weight)

        return NeatObject(genes=res_genes, weights=res_weights)

    def _weight_genes(self, adata, coef, max_normalize=True, pcs=None):
        if pcs is not None :
            tmp_coef = np.zeros_like(coef)
            tmp_coef[pcs] = coef[pcs]
            coef = tmp_coef
        loadings = adata.varm["PCs"] * coef
        weights = np.sum(loadings, axis = 1)
    
        if max_normalize:
            weights = weights / np.max(weights)

        genes = adata.var.index.values
        idx = weights.argsort()[::-1]
        
        genes, weights = genes[idx], weights[idx]
        
        return genes, weights
    
    def _select_gene(self, gene, weight, pattern, exclude_gene=True):
        """
        exclude or keeps genes that have a certain string pattern
        """
        if pattern is None :
            return True
        
        if (exclude_gene and not (pattern in gene)) or (not exclude_gene and (pattern in gene)):
            return True
        return False

    def _filter_gene(self, gene, weight, max_threshold=1, min_threshold=0.5):
        if weight > min_threshold and weight <= max_threshold:
            return True
        return False
