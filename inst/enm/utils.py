import os 
import networkx as nx
import pandas as pd 
import prody
import copy
import numpy as np
import random

from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.anno.gaf_reader import GafReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS


#from .Enm import *

def create_goea(gaf = '../data/raw/ontology/sgd.gaf', obo_fname = '../data/raw/ontology/go-basic.obo', background='../data/raw/ontology/sgd_costanzogenes', sgd_info_tab='../data/raw/ontology/SGD_features.tab',species='yeast',goset= ['BP'] , methods =['fdr'] , **kwargs):
    
    #background='../data/raw/ontology/sgd_costanzogenes'
    obodag = GODag(obo_fname,optional_attrs={'relationship'})
    objanno = GafReader(gaf,namespaces=set(goset))
#    ns2assoc = objanno.get_ns2assc()

    ns2assoc_excl = objanno.get_ns2assc( ev_exclude = {'HGI' , 'IGI'})

    bg = pd.read_csv(background, header=None)
    bg_list = list(bg.iloc[:, 0])  # .to_list()
    geneids = bg_list
    sgd_info = pd.read_csv(sgd_info_tab, sep='\t',header=None)

    geneid2name = pd.Series(sgd_info.iloc[:,3].values,index=sgd_info.iloc[:,0]).to_dict()
    goeaobj = GOEnrichmentStudyNS(
        geneids,  # List of mouse protein-coding genes
        ns2assoc_excl,  # geneid/GO associations
        obodag,  # Ontologies
        propagate_counts=True,
        relationships=True,
        alpha=0.1,  # default significance cut-off
        methods=methods, prt=None)

    return goeaobj, geneid2name , obodag #, objanno, ns2assoc, ns2assoc_excl

def goea_to_pandas(goea_results_sig, geneid2name):
    """ Converts goea object from goatools GO enrichment test to a Pandas dataframe
    :param goea_results_sig: Significant GO term objects
    :type goea_results_sig: list of GOEnrichmentStudy
    :return: Dataframe
    :rtype: Pandas DataFrame
    """
    if len(goea_results_sig) == 0 :
        return None
    go_df_n = pd.DataFrame([[getattr(i,x) for x in i.get_prtflds_default()] for i  in goea_results_sig], columns=goea_results_sig[0].get_prtflds_default())
    orf_names = []
    for i in go_df_n.study_items:
        orf_names.append([geneid2name[_id] for _id in i])
    go_df_n.study_items = orf_names
    if 'p_fdr' in go_df_n.columns:
        go_df_n['p_fdr_fix'] = (go_df_n['p_fdr']*500+1)/501
    return go_df_n 


def query_goatools(query, goea,geneid2name):
    """get query dataframe and goa files and return enrichments

    :param query: query gene dataframe
    :type query: pandas dataframe
    :param goea: goa object that will be used to run gene ontology analysis using GOAtools
    :type goea: goatools.goea.go_enrichment_ns.GOEnrichmentStudyNS
    :param geneid2name: dictionary to map geneids to gene names. needed to convert systematic names to intended names
    :type geneid2name: dict
    :return: enrichment dataframe 
    :rtype: pandas dataframe
    """
    query_gene_ids = [key for key,value in geneid2name.items() if value in query.loc[:,'Systematic gene name'].unique()]
    goea_res_all = goea.run_study(query_gene_ids)
    goea_res_sig = [r for r in goea_res_all if r.p_fdr<0.098]
    go_df_sensor = goea_to_pandas(goea_res_sig, geneid2name)
    return go_df_sensor

def combine_data(list_of_dfs):
    """given a list of dfs return a pandas dataframe by concatenating them

    :param list_of_dfs: a list of pandas dataframes, can contain None
    :type list_of_dfs: list
    :return: pandas DataFrame
    :rtype: DataFrame
    """
    df_res_concat = pd.concat(list_of_dfs,ignore_index=False,keys=range(len(list_of_dfs)))
    df_res_concat = (df_res_concat.
                     reset_index(inplace=False).
                     rename(columns = {'level_0':'orf_name_id'}).
                    drop('level_1',axis=1))

    df_res_concat
    return df_res_concat
def sample_nodes(sample_space, size=1):
    """get a list of nodes from the sample space with given size

    :param sample_space: list of nodes
    :type sample_space: list
    :param size: number of nodes to be sampled, defaults to 1
    :type size: int, optional
    :return: a random list of nodes
    :rtype: list
    """
    return   np.random.choice(sample_space, size=size, replace=False)

def get_degree_distribution(df):
    """given a dataframe with a 'deg' column, find the counts of each degree value

    :param df: a dataframe with a 'deg' column representing degree of nodes
    :type df: pandas DataFrame
    :return: a dictionary showing counts of each degree
    :rtype: dict
    """
    return df.groupby('deg').count().to_dict()['orf_name']   

def sample_nodes_with_degree(gnm_df , nodes_df ):
    """sample nodes with given degree

    :param gnm_df: a dataframe of all nodes with 'deg' column
    :type gnm_df: pandas DataFrame
    :param nodes_df: a dataframe of nodes to calculate expected degree distribution
    :type nodes_df: pandas DataFrame
    """
    deg_dist = get_degree_distribution(nodes_df)
    sampled_nodes = []
    for deg, count in deg_dist.items():
        sample_space = gnm_df.loc[(gnm_df.deg == deg)&(gnm_df.orf_name.isin(nodes_df.orf_name.values)==False),'orf_name'].values
        nds = sample_nodes(sample_space,size=count)
        sampled_nodes.extend(nds)
#        print(nds, deg, count, len(nds))
    return gnm_df.loc[gnm_df.orf_name.isin(sampled_nodes)]

def get_in_to_out_edge_ratio(G, nodes_df):
    """calculate edge ratio of number of edges within the nodes in nodes_df and
    to number of all edges these nodes have

    :param G: full network
    :type G: networkx graph
    :param nodes_df: subset of nodes with orf_name showing node names
    :type nodes_df: pandas DataFrame
    :return: ratio of between edges to total number of edges
    :rtype: float
    """
    btw_edges = len(nx.induced_subgraph(G,nodes_df.orf_name).edges)
    flat_list_ego = np.unique([item for sublist in [list(nx.ego_graph(G,i,radius=1).nodes) for i in nodes_df.orf_name.values] for item in sublist])
    total_edges = len(nx.induced_subgraph(G,flat_list_ego).edges)- len(nx.induced_subgraph(G,np.setdiff1d(flat_list_ego,nodes_df.orf_name)).edges)
    #rat = len([i for i in flat_list_ego if i in sensors['orf_name'].values ])/len([i for i in flat_list_ego if i not in sensors['orf_name'].values ])
    rat = btw_edges/total_edges #len([i for i in flat_list_ego if i in nodes_df['orf_name'].values ])/len(flat_list_ego)
    return rat

def get_random_in_to_out_ratio(gnm_df, df , G):
    """generate null distribution of within edge ratio

    :param gnm_df: a dataframe containing all nodes in G with 'deg' column
    :type gnm_df: pandas DataFrame
    :param df: a dataframe to be used to calculate wanted degree distribution
    :type df: pandas DataFrame
    :param G: full network
    :type G: networkx graph
    :return: a list of ratios
    :rtype: list
    """
    random_ratio_list = []
    for i in range(10000):
        rand_nodes = sample_nodes_with_degree(gnm_df, df)
    #    print(len(rand_nodes))
        random_ratio_list.append(get_in_to_out_edge_ratio(G, rand_nodes))
    return random_ratio_list

def get_subnetwork(gc, subset, radius = 1):
    neighbors =[[n for n in nx.ego_graph(gc, i, radius=radius).nodes] for i in subset]
    flat_list = [item for sublist in neighbors for item in sublist]
    flat_list.extend(subset)
    sub_gc = nx.induced_subgraph(gc, flat_list).copy()
    for (n,d) in sub_gc.nodes(data=True):
        if 'pos' in d.keys():
            del d["pos"]
    return sub_gc


# def get_maximal_subsets(sets):
#     sets = sorted(map(set, sets), key=len,reverse=True)
#     maximal_subsets = []
#     for s in sets:
#         if not any(maximal_subset.issuperset(s) for maximal_subset in maximal_subsets):
#             maximal_subsets.append(s)

#     return maximal_subsets

# def jaccard_index(G, a,b):
#     a_neigh = [i for i in nx.neighbors(G, a)]
#     b_neigh = [i for i in nx.neighbors(G, b)]
#     union = np.union1d(a_neigh,b_neigh)
#     intersect = np.intersect1d(a_neigh,b_neigh)
#     jaccard = len(intersect)/len(union)
#     return jaccard

def get_max_path(enm, target_cluster, source_cluster):
    """Calculate max information carrying path between 2 clusters by finding all possible nodes and returns based on maximum distance

    :param enm: this will be used to calculate prs weighted paths
    :type enm: Enm object
    :param target_cluster: cluster which contains possible target nodes
    :type target_cluster: networkx graph
    :param source_cluster: cluster which contains possible source nodes
    :type source_cluster: networkx graph
    """
    wmin = np.inf
    lmin = 0
    for source in source_cluster.nodes:
        for target in target_cluster.nodes:
            w, l1 = enm.get_prs_weighted_path(source,target)
            #print(w)
            num_of_sources = len([i for i in l1 if i in source_cluster.nodes])
            if num_of_sources > 1 :
                continue
            if w < wmin :
                lmin=l1
    #if lmin == 0 :
     #   lmin = l1
    return wmin, lmin
def get_path_positions(enm, sensors_sub, effectors_sub):
    """
        calculate positions for effectors and sensors
        choose a random source and target form these, respectively
        get PRS path between source and target
        put the nodes on the path linearly in between clusters
        change the source and target positions based on the initial cluster positions

        :param enm: this will be used to calculate prs weighted paths
        :type enm: Enm object
    """
    sensor_pos = nx.spring_layout(sensors_sub,center=(5,0))
    effector_pos = nx.spring_layout(effectors_sub, center=(5,10),scale=2)
    wmin, l1 = get_max_path(enm, sensors_sub, effectors_sub)
#    source = random.choice([i for i in effector_pos.keys()])
#    target = random.choice([i for i in sensor_pos.keys()])
    #e_pcc.prs_mat_df=pd.DataFrame(e_pcc.prs_mat,columns=e_pcc.nodes, index=e_pcc.nodes)
#    l1 = e_pcc.get_prs_weighted_path(source,target)[1]
    #l1.extend(['fum1', 'fbp1'])
    sub = nx.induced_subgraph(enm.graph_gc, l1)
    path_init_pos=dict(zip(l1[1:-1],np.array(tuple(zip(np.repeat(5,(len(l1)-2)), np.linspace(effector_pos[l1[0]][1]-1,sensor_pos[l1[-1]][1],len(l1))[1:-1])))))
    
    for i in sub.nodes:
        if i in effector_pos.keys():
            path_init_pos[i]=effector_pos[i]
        if i in sensor_pos.keys():
            path_init_pos[i]=sensor_pos[i]
    return sensor_pos, effector_pos, path_init_pos,sub, l1, wmin

def convert_rpy2_to_pandas(df):
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
        pd_from_r_df = ro.conversion.rpy2py(df)

    return pd_from_r_df

def get_result_dfs(fname, thr_list, default_thr=None, folder_prefix= '../data/interim'):
    def read_csv(fname):
        try:
            df = pd.read_csv(fname)
        except pd.errors.EmptyDataError:
            print(f"{fname} is empty")
            df = None
        return df

    dfs = {thr: read_csv(f"{folder_prefix}_{thr}/{fname}.csv") for thr in thr_list if thr!=default_thr}
    if default_thr is not None:
        dfs[default_thr] = pd.read_csv(f"{folder_prefix}/{fname}.csv")
    if default_thr not in thr_list and default_thr is not None:
        thr_list.insert(0, default_thr)                               
    return dfs
