import pickle
import os
import networkx as nx
import pandas as pd
import prody
import numpy as np
import scipy.cluster.hierarchy as sch
import copy
from tqdm import tqdm

from .visualize import plot_network_spring,heatmap_annotated
from .utils import query_goatools


class Enm():
    """This object is wrapper around prody GNM object and networkx analysis
    """

    def __init__(self, name):
        """Constructor

        :param name: name of the object for future reference, can be arbitrary
        :type name: string"""
        self.name = name
        self.figure_path = "../reports/figures/"
        self.output_path = "../data/interim/"
        self.rewired = False
        self.rewire_df = None
        self.arr = None
        self.e_list = None
        self.df = None
        self.G = None
        self.graph_gc = None
        self.nodes=None
        self.degree = None
        self.gnm = None
        self.coll=None
        self.coll_index_sorted = None
        self.prs_mat = None
        self.prs_mat_df  = None
        self.prs_mat_cl = None
        self.L = None
        self.go_groups = {}

    def print_name(self):
        """This prints name variable
        """
        print(self.name)

    def spring_pos(self, seed = None):
        """Create spring layout of the giant component network and assign positions to nodes
        """
#        rs = np.random.RandomState(2)
        try:
            pos = nx.spring_layout(self.graph_gc, k=0.6,
                                   scale=4, iterations=200, seed=seed)
            nx.set_node_attributes(self.graph_gc, pos, 'pos')
        except AttributeError:
            raise(
                'Giant component is not set yet. First call read network or giant component')

    def read_network(self, path, **kwargs):
        """Read network file and assign it to object

        :param path: Network file path
        :type path: string

        """
        sep = kwargs.pop('sep', None)
        _, ext = os.path.splitext(path)
        # return fname, ext
        if ext == '.csv' or sep == ',':
            nw_read = pd.read_csv(path)
            # Create network ##########
            G = nx.from_pandas_edgelist(
                nw_read, source='gene1', target='gene2')
        elif ext == '.gpickle':
            G = nx.read_gpickle(path)
        elif ext == '.tsv' or sep == '\t':
            nw_read = pd.read_csv(path, sep='\t')
            # Create network ##########
            G = nx.from_pandas_edgelist(
                nw_read, source='gene1', target='gene2')
        elif ext == '.gml':
            G = nx.read_gml(path)
        self.G = G
        self.giant_component()

    def giant_component(self, **kwargs):
        """From the graph variable create giant component and assing nodes and degree to Enm object
        """

        Gc = max([self.G.subgraph(c).copy()
                  for c in nx.connected_components(self.G)], key=len)
        self.graph_gc = Gc
        nodes = [_ for _ in self.graph_gc.nodes()]
        self.nodes = nodes
        degree = [deg for id, deg in list(self.graph_gc.degree)]
        self.degree = degree
        self.laplacian_matrix(**kwargs)

    def laplacian_matrix(self, normalized=False,**kwargs):
        """get Laplacian matrix of the giant component. Wrapper around networkx.laplacian_matrix
        """
        if normalized:
            self.L = nx.normalized_laplacian_matrix(self.graph_gc,weight=None).todense()
        else:
            self.L = nx.laplacian_matrix(self.graph_gc, weight=None).todense()

    def get_gnm(self, **kwargs):
        """Calculate GNM modes and collectivity values
        """
        if self.L is None:
            self.laplacian_matrix(**kwargs)
        gnm = prody.GNM()
        gnm.setKirchhoff(self.L)

        gnm.calcModes(n_modes=None, zeros=False)
        #gnm.nodes = nodes
        #gnm.degrees = degree
        self.gnm = gnm
        # showCrossCorr(gnm)
        #sqf_orig = prody.calcSqFlucts(self.gnm)
        coll = prody.calcCollectivity(self.gnm)
        self.coll = coll
        coll_index_sorted = sorted(
            range(len(self.coll)), key=lambda k: self.coll[k], reverse=True)

        self.coll_index_sorted = coll_index_sorted

    def get_prs(self,**kwargs):
        """Calculate Perturbation response matrix
        """
        try:
            no_diag = kwargs.pop('no_diag','True')
            prs_mat, _, _ = prody.calcPerturbResponse(
                self.gnm, no_diag=no_diag)
        except Exception as e:
            raise AttributeError('GNM is not calculated yet. Call get_gnm() first')
        self.prs_mat = prs_mat
        self.prs_mat_df = pd.DataFrame(prs_mat,columns=self.nodes, index=self.nodes)

    def create_df(self):
        """Create an overall dataframe to store network related and GNM related values
        """
        df = pd.DataFrame()
        df['orf_name'] = self.nodes
        df['deg'] = np.diag(self.L)
        eff_orig = np.sum(self.prs_mat, axis=1)
        sens_orig = np.sum(self.prs_mat, axis=0)
        eigvecs_df = pd.DataFrame(self.gnm.getEigvecs(
        )[:, self.coll_index_sorted[:8]], columns=[f'eig_{i}' for i in range(8)])

        df_ = df#pd.merge(df, eigvecs_df, left_index=True, right_index=True)
        df_['eff'] = eff_orig
        df_['sens'] = sens_orig
        eigenvector_centr = nx.eigenvector_centrality_numpy(self.graph_gc)
        closeness_centr = nx.closeness_centrality(self.graph_gc)
#        df_['btw'] = betweenness_nx(self.graph_gc, normalized=True)
        df_['trans'] = list((nx.clustering(self.graph_gc)).values())
        df_['eigenvec_centr'] = [eigenvector_centr[i]
                                 for i in eigenvector_centr]
        df_['closeness_centr'] = [closeness_centr[i] for i in closeness_centr]
        df_['smallest_eigenvec'] = self.gnm.getEigvecs()[:, 0]
        self.df = df_

    def gnm_analysis(self, **kwargs):
        """Wrapper to run gnm, prs and create_df
        """
        self.get_gnm(**kwargs)
        self.get_prs(**kwargs)
        self.create_df()

    def get_category(self, strain_ids_df_file):
        """ Uses costanzo strain id file to determine categories. File is created by pgsNetwork project earlier, depreceated
        """
        strain_ids = pd.read_csv(strain_ids_df_file)
        combined_df = pd.merge(self.df, strain_ids, left_on='orf_name',right_on='Allele Gene name')
        combined_df['group']=np.where(combined_df.cat.isna(),'essential','nonessential')
        combined_df = combined_df.fillna({'cat':'essential'})
        cat_change_dict = {'essential': 'Essential',
                  'na.nq.nxes':'Nonessential\nquery and array',
                  'nxes.only':'Nonessential \nquery crossed \n with Essential',
                  'nq.nxes':'Nonessential query',
                  'na.nq':'Nonessential\nquery and array',
                  'na.nxes':'Nonessential\nquery and array',
                  'nq.only':'Nonessential query',
                  'na.only':'Nonessential array'}
        combined_df['cat_']= combined_df['cat'].map(cat_change_dict)
        self.df = combined_df
        #db_connection_str = 'mysql+pymysql://oma21:dktgp2750@localhost:3306/ANNE'
        #db_connection = sql.create_engine(db_connection_str)
        #db_df = pd.read_sql('SELECT * FROM SUMMARY_2012', con=db_connection)
        #combined_df = pd.merge(combined_df, db_df, left_on='Systematic gene name', right_on='orf_name')

    def get_sensor_effector(self, use_threshold=True, quantile_threshold=0.99):
        """create sensor and effector sub dataframes using PRS matrix clusters or effectiveness sensitivity thresholds

        :param use_threshold: if true, the nodes with effectiveness/sensitivity above quantile threshold will be taken as effectors and sensors, defaults to True
        :type use_threshold: bool, optional
        :param quantile_threshold: any effectiveness/sensitivity value above this quantile will be important, defaults to 0.99
        :type quantile_threshold: float, optional
        """
        if self.df is None:
            self.gnm_analysis()
        if use_threshold:
            sensors_df = self.df.loc[self.df.sens>np.quantile(self.df.sens,quantile_threshold)]
            effectors_df = self.df.loc[self.df.eff>np.quantile(self.df.eff,quantile_threshold)]
        else:
            if self.prs_mat_cl is None:
                self.cluster_matrix(self.prs_mat)
            sensors_df = self.df.iloc[self.get_clustered_nodes('column'),:]
            effectors_df = self.df.iloc[self.get_clustered_nodes('row'),:]
        self.sensors_df = sensors_df
        self.effectors_df = effectors_df

    def analyze_components_biology(self, goea, geneid2name, sensors=True):
        """Use sets of sensors or effectors to find connected components among them and calculate GO term enrichments 

        :param goea: GOAtools object for GO enrichment analysis
        :type goea: 
        :param geneid2name: dictionary for name convention
        :type geneid2name: dict
        :param sensors: if true uses self.sensors_df, else uses self.effectors_df for analysis, defaults to True
        :type sensors: bool, optional
        """
        try:
            if sensors:
                col_name = 'sensor_cluster'
                df = self.sensors_df
            else:
                col_name = 'effector_cluster'
                df = self.effectors_df
        
        except AttributeError :
                print(
                    'Sensors or effectors dataframe might be missing. Make sure to call get_sensor_effector() first')
        components = [i for i in nx.connected_components(nx.induced_subgraph(self.graph_gc, df.orf_name))]
        dd = {}
        for i,j in dict(zip(range(len(components)),components)).items():
            id_ = i if len(j)>=3 else None
            for item in j:
                dd[item]=id_
        df.loc[:,col_name] = df['orf_name'].map(dd)
        go_terms = {}
        for i in df.loc[pd.isna(df[col_name])==False,col_name].unique():
            go_terms[i] = query_goatools(df.loc[df[col_name]==i,:], goea, geneid2name)
        df['go_group'] = None
        if sensors:
            self.go_groups['sensors_go_groups'] = go_terms
        else:
            self.go_groups['effectors_go_groups'] = go_terms
        for i,j in go_terms.items():
            if j is not None:
                df.loc[df[col_name]==i,'go_group']=j.iloc[0,3]

    def get_clustered_nodes(self, dimension='row'):
        """Uses the prs matrix clustering for row or columns separately and takes the smaller cluster of 2nd level as the effector or sensor nodes cluster

        :param dimension: which linkage should be used, defaults to 'row'
        :type dimension: str, optional
        :return: node ids
        :rtype: list
        """
        if dimension == 'row':
            root = self.root_row
        else:
            root = self.root_col
        right = root.right
        left = root.left
        if right.get_count() < left.get_count():
            nds = right.pre_order()
        else:
            nds = left.pre_order()
        return nds
         
    def get_prs_weighted_path(self, source, target=None, weight='weight', node_weight=None):
        """
            calculate shortest path using the PRS matrix values caused by perturbations on source node
            This uses an adjusted version of networkx dijkstra functions to incorporate node_weight parameter
            networkx source code could be installed using
            pip install git+git://github.com/oacar/networkx@node_weight#egg=networkx  -U 
        """
        if nx.get_edge_attributes(self.graph_gc,'weight') == {}:
            nx.set_edge_attributes(self.graph_gc,0,'weight')
        if node_weight is None:
            node_weight = (1/self.prs_mat_df.loc[source,:]).to_dict()
            node_weight[source] = 1
        if target is not None:
            distances = nx.single_source_dijkstra(self.graph_gc, source=source,target=target, weight=weight, node_weight=node_weight)
        else:
            distances = nx.single_source_dijkstra(self.graph_gc, source=source, weight=weight, node_weight=node_weight)
        return distances

    def cluster_matrix(self, mat, method='ward', distance_metric='seuclidean', quantile_threshold = 0.95, cluster_normalized=True, show_normalized=True, optimal_ordering=True):
        """create row and column linkage for given matrix `mat`

        :param mat: the input matrix to be clustered
        :type mat: numpy matrix, 2 dimensional
        :param method: clustering method. see `scipy.cluster.hierarchy.linkage` for details, defaults to 'ward'
        :type method: str, optional
        :param quantile_threshold: any values above this quantile will be equal to quantile value , defaults to 0.95
        :type quantile_threshold: float, optional
        :param cluster_normalized: whether to cluster the matrix after using quantile threshold or not, defaults to True
        :type cluster_normalized: bool, optional
        :param optimal_ordering: use optimal leaf ordering for linkage calculation. see, `scipy.cluster.hierarchy.linkage`.  defaults to True
        :type optimal_ordering: bool, optional
        """
        q99 = np.quantile(mat, quantile_threshold)
        mat_cl = copy.deepcopy(mat)
        mat_cl[mat_cl > q99] = q99
        if cluster_normalized:
            row_dist = sch.distance.pdist(mat_cl, metric=distance_metric)
            col_dist = sch.distance.pdist(mat_cl.T, metric=distance_metric)
            row_linkage = sch.linkage(row_dist, method=method,optimal_ordering=optimal_ordering)
            col_linkage = sch.linkage(col_dist, method=method,optimal_ordering=optimal_ordering)
        else:
            row_dist = sch.distance.pdist(mat, metric=distance_metric)
            col_dist = sch.distance.pdist(mat.T, metric=distance_metric)
            row_linkage = sch.linkage(row_dist, method=method,optimal_ordering=optimal_ordering)
            col_linkage = sch.linkage(col_dist, method=method,optimal_ordering=optimal_ordering)

        root_row, tree_list_row = sch.to_tree(row_linkage, True)
        root_col, tree_list_col = sch.to_tree(col_linkage, True)

        # ,optimal_ordering=True)
        self.row_linkage = row_linkage
        self.col_linkage = col_linkage
        self.row_dist = row_dist
        self.col_dist = col_dist
        self.prs_mat_cl = mat_cl
        self.root_row = root_row
        self.root_col = root_col
        self.q99 = q99

    def get_rwr_mat(self, c=0.15):
        """use pyrwr to calculate Random walk with rewiring for each nodes in the network and create rwr_mat/rwr_mat_df 

        :param c: restart probability, defaults to 0.15
        :type c: float, optional
        """
        from pyrwr.rwr import RWR
        rwr = RWR()
        rwr.A = nx.adj_matrix(self.graph_gc, weight=None)
        rwr.m , rwr.n = rwr.A.shape
        rwr_mat = np.zeros(rwr.A.shape)
        rwr.base = 0
        for i in range(rwr_mat.shape[0]):
            seed = i
            r = rwr.compute(seed, c=c)
            rwr_mat[i,:] = r
        self.rwr_mat = rwr_mat
        self.rwr_mat_df = pd.DataFrame(rwr_mat,columns=self.nodes, index=self.nodes)
        self.pagerank = np.mean(rwr_mat, axis=0)

    def plot_network_spring(self, **kwargs):
        """Plot network with spring layout
        """
        Gc = self.graph_gc
        if nx.get_node_attributes(Gc, 'pos') == {}:
            self.spring_pos()

        return plot_network_spring(Gc, self.figure_path, **kwargs)


    def simulate_rewire(self, rewired=False, rewire_df_name=None, arr_name=None, **kwargs):
        """Wrapper around enm.simulate_rewire function out of the class
        """
        output_path = kwargs.pop('output_path', self.output_path) 

        arr, rewire_df, e_list = simulate_rewire(
            self.graph_gc, rewired, rewire_df_name, arr_name, output_path=output_path, **kwargs)
        nodegseq= kwargs.pop('nodegseq',False)
        if nodegseq:
            self.arr_nodegseq = arr
            self.rewire_df_nodegseq = rewire_df
            self.e_list_nodegseq = e_list
        else:
            self.arr = arr
            self.rewire_df = rewire_df
            self.e_list = e_list

    def heatmap_annotated(self, **kwargs):
        """Plot PRS heatmap with clustering dendrograms
        """
        if self.prs_mat_cl is None:
            self.cluster_matrix(self.prs_mat, **kwargs)
        return heatmap_annotated(self.prs_mat, self.prs_mat_cl, self.figure_path, self.row_linkage, self.col_linkage, **kwargs)

    # def save(self, **kwargs):
    #     filename = kwargs.pop('filename', self.name)
    #     with open(f"../data/interim/{filename}.pickle", 'wb') as handle:
    #         pickle.dump(self, handle)


def rewire_network(Gc, **kwargs):
    """This is wrapper around networkx's rewire functions for my purposes. It can rewire by keeping degree, or if not it will generate a network with same number of nodes and arbitrary number of edges using barabasi_albert_graph or erdos_renyi_graph

    :param Gc: Network to be rewired    
    :type Gc: networkx object
    :return: Rewired network
    :rtype: networkx object
    """
    nodegseq = kwargs.pop('nodegseq', False)
    if nodegseq:
        random_network_type = kwargs.pop('random_network_type', 'ba')
        if random_network_type == 'ba':
            Gc_rewired = nx.barabasi_albert_graph(n=len(Gc.nodes), m=7)
        elif random_network_type == 'er':
            Gc_rewired = nx.erdos_renyi_graph(len(Gc.nodes), 0.004)
    else:
        Gc_rewired = Gc.copy()
        minimum_swaps = 10*len(Gc_rewired.edges)
        current_swaps = 0
        while current_swaps < minimum_swaps:
            remaining_swaps = max(minimum_swaps-current_swaps,100)
            swp_count = nx.connected_double_edge_swap(Gc_rewired, remaining_swaps)
            current_swaps+=swp_count

    return Gc_rewired


# def rewire_network():


# TODO add option to save the rewire results
def simulate_rewire(Gc, rewired=False, rewire_df_name=None, arr_name=None, **kwargs):
    """This function reads rewired network GNM data or calls rewire function

    :param Gc: network to be rewired
    :type Gc: networkx object
    :param rewired: Is there a rewired GNM data available, defaults to False
    :type rewired: bool, optional
    :param rewire_df_name: If rewired is True, related dataframe path should be given, defaults to None
    :type rewire_df_name: string, optional
    :param arr_name: If rewired is True, related numpy array txt file path should be given, defaults to None
    :type arr_name: string, optional
    :raises ValueError: If rewired is True and rewire_df_name or arr_name is not given, raises error
    :return: arr
    :rtype: numpy array
    :return: rewire_df
    :rtype: pandas dataframe    
    """
    from scipy.stats import pearsonr, spearmanr
    from tqdm import tqdm
    save = kwargs.pop('save', False)
    output_name = kwargs.pop('output_name', 'rewire_data')
    output_path = kwargs.pop('output_path', '../data/interim/') 
    if rewired:
        # rewire_df_name = kwargs.pop('rewire_df_name',None)
        # arr_name = kwargs.pop('arr_name',None)
        if rewire_df_name is None or arr_name is None:
            raise ValueError(
                'Rewired dataframe path or arr_name is not given ')
        rewire_df = pd.read_csv(rewire_df_name)
        arr = np.load(arr_name)
    else:

        sim_num = kwargs.pop('sim_num', 10)
        arr = np.zeros((len(Gc.nodes), len(Gc.nodes), int(sim_num)))
        # pearsonr(eff,degree)[0]
        rewire_df = pd.DataFrame(columns=['eff_deg_pearson', 'sens_deg_pearson', 'eff_deg_spearman', 'sens_deg_spearman',
                                          'eff_trans_pearson', 'sens_trans_pearson', 'eff_trans_spearman',
                                          'sens_trans_spearman'])
        e_list = []
#        for i in tqdm(range(int(sim_num))):
        success = 0
        i=0
        maxtry=sim_num*10
        while success<(sim_num) and i<=maxtry:
            print(i)
            i=i+1
            try:
                Gc_rew = rewire_network(Gc, **kwargs)
#                print(Gc_rew)
                enm_rew = Enm('rewired')

                
                enm_rew.G = Gc_rew
                enm_rew.giant_component()
#                print(enm_rew.graph_gc)
                enm_rew.gnm_analysis(**kwargs)
                res = enm_rew.prs_mat
#                print(res)
                #Gc_rew = enm_rew.graph_gc
                degree = enm_rew.degree
                e_list.append(enm_rew)
                success = success+1
                #i=i+1
            except Exception as e:
                print('error')
                continue

            arr[:, :, (success-1)] = res
            eff_rew = enm_rew.df.eff.values
            sens_rew = enm_rew.df.sens.values
            # eff_hist_list.append(eff_rew)
            # sens_hist_list.append(eff_rew)
            # betweenness_nx(Gc_rew, normalized=True)
            # betweenness = enm_rew.df.btw.values
            clustering_coeff = enm_rew.df.trans.values
            rewire_df_itr = pd.DataFrame([[pearsonr(eff_rew, degree)[0], pearsonr(sens_rew, degree)[0],
                                           spearmanr(eff_rew, degree)[0], spearmanr(
                                               sens_rew, degree)[0],
                                        #    pearsonr(eff_rew, betweenness)[0], pearsonr(
                                        #        sens_rew, betweenness)[0],
                                        #    spearmanr(eff_rew, betweenness)[0], spearmanr(
                                        #        sens_rew, betweenness)[0],
                                           pearsonr(eff_rew, clustering_coeff)[0], pearsonr(
                                               sens_rew, clustering_coeff)[0],
                                           spearmanr(
                                               eff_rew, clustering_coeff)[0],
                                           spearmanr(sens_rew, clustering_coeff)[0]]],
                                         columns=['eff_deg_pearson', 'sens_deg_pearson', 'eff_deg_spearman',
                                                  'sens_deg_spearman',
                                                  'eff_trans_pearson', 'sens_trans_pearson', 'eff_trans_spearman',
                                                  'sens_trans_spearman'])
            rewire_df = rewire_df.append(rewire_df_itr)
            # if i % 1 == 0:
            #     print(i)
        # rewire_df.to_csv(outpath+'/rewire_df.csv')

        if save:
            rewire_df.to_csv(f"{output_path}/{output_name}.csv")
            np.save(f"{output_path}/{output_name}.npy", arr)
            with open(f"{output_path}/{output_name}.pickle", 'wb') as handle:
                pickle.dump(e_list, handle)
    return arr, rewire_df, e_list
