import numpy as np
from sklearn.linear_model import LogisticRegressionCV
from sklearn import preprocessing
from sklearn.decomposition import PCA

class NeatObject:
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
        self.weights = {}
        self.max_norm_weights = {}
        self.global_max_norm_weights = {}
        self.ordered_pcs = {}
        self.label_encoder = preprocessing.LabelEncoder()

    def fit(self, adata, fit_on=None, labels_key="leiden", max_pcs=None):
        """
        Fit the model to find the importance of every PC for each clister and 
        weights genes accordingly
        """
        if labels_key in adata.obs:
            labels = adata.obs[labels_key].values
        elif labels_key in adata.obsm:
            labels = adata.obsm[labels_key].values
        elif labels_key in adata.uns:
            labels = adata.uns[labels_key].values
        else:
            raise KeyError("Key %s not found in obs, obsm or uns" % labels_key)

        self.label_encoder.fit(labels)

        if fit_on is None:
           the_repr = adata.obsm["X_pca"]
        else :
            print("running pca...")
            if max_pcs is None:
                raise ValueError("With fit_on, max_pcs can't be None")
            pca = PCA(n_components=max_pcs)
            pca.fit(fit_on)
            the_repr = pca.transform(fit_on)
            print("\tdone.")

        stop = int(the_repr.shape[0] * 0.2)
        X_test = the_repr[:stop]
        y_test = self.label_encoder.transform(labels[:stop])#.astype('int32')

        X_train = the_repr[stop:]
        y_train = self.label_encoder.transform(labels[stop:])#.astype('int32')

        self.clf.fit(X_train, y_train)
        
        self.train_acc = self.clf.score(X_train, y_train)
        self.test_acc = self.clf.score(X_test, y_test)
        print("Train acc", self.train_acc)
        print("Test acc", self.test_acc)

        self.coefs = self.clf.coef_
    
        global_max = None
        labels = np.unique(y_test)
        for label, aidi in zip(labels, self.label_encoder.inverse_transform(labels)):
            self.ordered_pcs[aidi] = np.argsort(self.coefs[label])[::-1]
            pcs = None
            if max_pcs is not None :
                pcs = self.ordered_pcs[aidi][:max_pcs]
            
            if fit_on is None:
                PCs = adata.varm["PCs"]
            else:
                PCs = pca.components_.T

            genes, weights, max_norm_weights = self._weight_genes(adata, PCs, self.coefs[label], pcs=pcs)
            curr_max = np.max(weights)
            if global_max is None or global_max < curr_max:
                global_max = curr_max
            
            self.genes[aidi] = genes
            self.weights[aidi] = weights
            self.max_norm_weights[aidi] = max_norm_weights

        for aidi in self.label_encoder.inverse_transform(labels):
            self.global_max_norm_weights[aidi] = self.weights[aidi] / global_max
        return self
    
    def select(self, pattern="RP", exclude_genes=True, max_threshold=1, min_threshold=0.5):
        """Select only a subet of genes, returns a NeatObjet"""
        from collections import defaultdict
        # res_genes = defaultdict(list)
        # res_max_norm_weights = defaultdict(list)
        # res_weights = defaultdict(list)
        fields = ["genes", "weights", "max_norm_weights", "global_max_norm_weights"]
        final_res = { k: defaultdict(list) for k in fields}

        for aidi in self.genes.keys():
            for gene_id, (gene, weight) in enumerate( zip(self.genes[aidi], self.max_norm_weights[aidi]) ):
                if self._select_gene(gene, weight, pattern=pattern, exclude_gene=exclude_genes) \
                    and self._filter_gene(gene, weight, min_threshold=min_threshold):
                    for field in fields :
                        final_res[field][aidi].append(
                            getattr(self, field)[aidi][gene_id]
                        )
                    # res_genes[aidi].append(gene)
                    # res_weights[aidi].append(weight)

        # return NeatObject(genes=res_genes, weights=, max_norm_weights=res_weights)
        return NeatObject(**final_res)

    def _weight_genes(self, adata, PCs, coef, pcs=None):
        if pcs is not None :
            tmp_coef = np.zeros_like(coef)
            tmp_coef[pcs] = coef[pcs]
            coef = tmp_coef
        loadings = PCs * coef
        weights = np.sum(loadings, axis = 1)
    

        genes = adata.var.index.values
        idx = weights.argsort()[::-1]
        
        genes, weights = genes[idx], weights[idx]
        max_norm_weights = weights / np.max(weights)
        
        return genes, weights, max_norm_weights
    
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
