# Functions for replicating results in VariantInterpretationBenchmark.py
import numpy as np
# Classification performance

def get_metrics_fromscores(y_true, y_pred, threshold, scoring_dict):
    key = scoring_dict.keys()
    v_list = []
    for k in key:
        v = scoring_dict[k]
        v_list.append(v._score_func(y_true, (y_pred >= threshold).astype('int') ))
    return dict(zip(key,
                    v_list))


def get_metrics_frompred(y_true, y_pred, scoring_dict):
    key = scoring_dict.keys()
    v_list = []
    for k in key:
        v = scoring_dict[k]
        v_list.append(v._score_func(y_true, y_pred ))
    return dict(zip(key,
                    v_list))




# Prioritization performance
# Functions
# scores: not sorted values.
def dcg(y_true, rank=None):

    """Discounted cumulative gain at rank (DCG)"""
    if rank==None:
        rank=len(y_true)
    y_true = np.asarray(y_true)[:rank]
    n_relevances = len(y_true)
    if n_relevances == 0:
        return 0.

    discounts = np.log2(np.arange(n_relevances) + 2)
    return np.sum(y_true / discounts)

def ndcg(y_true, rank):
    """Normalized discounted cumulative gain (NDGC)"""
    best_dcg = dcg(sorted(y_true, reverse=True), rank)
    if best_dcg == 0:
        return 0.

    return dcg(y_true, rank) / best_dcg


def idcg(y_true, rank):
    return dcg(sorted(y_true, reverse=True), rank)

def ndcg_tie_aware(scores, relevances, k=None):
    sort_scores = list(np.flip(np.sort(scores)))
    sort_index = list(np.flip(np.argsort(scores)))
    sort_rel = list([relevances[i] for i in sort_index])
    if k is None:
        k = len(scores)
        return ndcg(sort_rel, k)

    sort_scores_k = sort_scores[:k]
    sort_rel_k = sort_rel[:k]
    if sort_scores_k[len(sort_scores_k)-1]==sort_scores[k]:
        nc = len([x for x in sort_scores[k:] if x == sort_scores[k]])
        sort_scores_k = sort_scores_k+sort_scores[k:k+nc]
        sort_rel_k = sort_rel[:k+nc]
    dcg_k = 0

    for s in np.unique(sort_scores_k):
        i_s = [i for i,x in enumerate(sort_scores) if sort_scores[i] == s]
        lindex = len(i_s)
        if i_s[len(i_s)-1] >= k:
            lindex=[i for i,x in enumerate(i_s) if x >= k][0]
        dcg_k = dcg_k + (1/len(i_s))*np.sum([sort_rel[i] for i in i_s])/np.sum(np.log2([x+2 for x in i_s[:lindex]]))

    id =idcg(sort_rel[:k], k)
    if id == 0:
        return 0
    return dcg_k/id

# Return a new index list that takes into account the fact the multiple variants could have the same score
def ranking(scores):
    desc_scores = np.flip(np.sort(scores))

    unique_scores = np.flip(np.unique(desc_scores))
    new_index = [0]*len(desc_scores)

    for i,u in enumerate(unique_scores):
        ind_u = [i for i in range(len(desc_scores)) if desc_scores[i]==u]
        if len(ind_u)==1:
            new_index[ind_u[0]] = ind_u[0]
        else:
            for j in ind_u:
                new_index[j] = np.median(ind_u)

    return new_index


def norm_discounted_cumulative_gain(y_true, scores, k=None, tie=False):
    if k is None:
        k = len(scores)

    # From highest to lowest score
    sorted_index = np.flip(np.argsort(scores))
    new_index = ranking(scores)
    # lst_i = np.flip(np.where(np.array(new_index)<=k)[0])[0]

    sorted_y_true = [y_true[i] for i in sorted_index[:k]]
    if tie:
        return(ndcg_tie_aware(scores, y_true, k))
    return(ndcg(sorted_y_true, rank=k))


# Tie-aware implementation of Precision at rank k from McSherry and Najork
# relevances -> y_pred
def precision_at_rank(relevances, scores, k=None):

    sort_scores = list(np.flip(np.sort(scores)))
    sort_index = list(np.flip(np.argsort(scores)))
    sort_rel = list([relevances[i] for i in sort_index])

    if k is None:
        k = len(scores)
        return np.sum(sort_rel) / k

    sort_scores_k = sort_scores[:k]
    if sort_scores_k[len(sort_scores_k)-1] == sort_scores[k]:
        # Tie
        tied_value = sort_scores[k]
        # Last index of not tied value across k
        li = sort_scores.index(tied_value)-1
        nc = len([x for x in sort_scores if x == tied_value])
        ntied_k = len([x for x in sort_scores_k if x==tied_value])
        prec_tie = ( np.sum( sort_rel[:li+1] ) + (ntied_k/nc)* np.sum(sort_rel[li+1:li+1+nc]) )/k
        return prec_tie

    else:

        return np.sum(sort_rel[:k])/k

# Recall Tie Aware (see https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/ecir2008.pdf)
# y_true: 0/1
# y_pred: prob predicted for the class 1
def recall_tie_aware(y_true, y_pred, k):
    last_index = k-1
    v = np.sort(y_pred)[::-1]
    rel = y_true[np.argsort(y_pred)][::-1]

    if v[last_index] != v[last_index+1]:
        # No tie
        return len(np.argwhere(rel[0:last_index+1]))/len(np.argwhere(rel))


    # Tied score value
    v_li = v[last_index]

    # First index of the tied elem across k
    tc = np.argwhere(v==v_li)[0][0]
    nc = len(np.argwhere(v==v_li))
    Rc = len(np.argwhere(rel[0:tc]==1))
    rc = len(np.argwhere(rel[np.argwhere(v==v_li)]==1))
    total_recall = len(np.argwhere(rel))

    return (Rc+(last_index-tc)*rc/nc)/(total_recall)

def pathogenic_rule(row):
    if row.nPVS > 0:
        if row.nPS >= 1:
            return True
        if row.nPM >= 2:
            return True
        if row.nPM == 1:
            if row.nPP == 1:
                return True
        if row.nPP >= 2:
            return True
        
    if row.nPS >= 2:
        return True
    
    
    
    if row.nPS == 1:
        if row.nPM >= 3:
            return True
        if row.nPM >= 2:
            if row.nPP >= 2:
                return True
        if row.nPM == 1:
            if row.nPP >= 4:
                return True
    
    
    return False


def likelypatho_rule(row):
    if row.nPS==1 and row.nPM==1:
        return True
    
    if row.nPS==1 and row.nPM <=2:
        return True
    
    if row.nPS==1 and row.nPP>=2:
        return True
    
    if row.nPM>=3:
        return True
    
    if row.nPM==2 and row.nPP>=2:
        return True
    
    if row.nPM==1 and row.nPP>=4:
        return True
    
    return False

def likelybenign_rule(row):
    
    if row.nBS==1 and row.nBP==1:
        return True
    
    if row.nBP>=2:
        return True
    
    return False

def benign_rule(row):
    
    if row.nBA==1:
        return True
    
    if row.nBS>=2:
        return True
    
    return False