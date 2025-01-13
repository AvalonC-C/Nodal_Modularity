# Nodal_Modularity
Nodal modularity (nQ) is a mathematical extension of classical modularity (Q) to individual nodes. That is, a node's contribution to overall segregation within a network. It has been implemented for both single and multi-layer networks. 

This code was created to be used with pre-calculated measures of modularity such as those from community_louvain.m of the Brain Connectivity Toolbox (https://sites.google.com/site/bctnet/list-of-measures?authuser=0) or the General Louvain method (https://github.com/GenLouvain/GenLouvain) which optimize modularity using either the single or multi-layer modularity quality functions. In these cases, sum(nQ)=Q. Nodal modularity can also be calculated with group/module assignments obtained via other methods of community detection. While potentially beneficial, the interpretability and connection to standard modularity may change significantly. 

# Citation

Campbell-Cousins, A., Guazzo, F., Bastin, M., Parra, M. & Escudero, J. Multiplex Nodal Modularity: A novel network metric for the regional analysis of amnestic mild cognitive impairment during a working memory binding task. *arXiv* [LINK] (2025)
