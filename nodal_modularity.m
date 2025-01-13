function [nQ, Q] = nodalQ(A, S, gamma, omega) 
%% 
% Calculate the per node contribution to modularity (nodal modularity or
% nQ) for a single or multi-layer network. 
%
% Version: 1.0, 27/11/2024 
%
% Avalon Campbell-Cousins - avalon.campbell-cousins@ed.ac.uk
% Javier Escudero - jescuder@ed.ac.uk
% 
%%%%%%%%%% Single-layer version %%%%%%%%%%
%
% [nQ] = nodalQ(A, S);
%
% Input:
%
% A - a weighted or binary nxn square and symmetric adjacency matrix (edges must be positive and undirected). Please ensure that A was used to
% calculate S apriori. Note: please contact me if you would like support
% for directed and/or signed networks.
% 
% S - a 1xn vector where entries correspond to the group assignment of
% nodes in A. S has been tested with the freely available BCT ver. 03-03-2019 (community_louvain.m in the Brain Connectivity Toolbox: https://sites.google.com/site/bctnet/list-of-measures?authuser=0) 
% and the General Louvain method ver. 2.2 (available here: https://github.com/GenLouvain/GenLouvain). 
%
% Output: 
%
% nQ - a 1xn vector where sum(nQ) = Q.
%
% Q - classical modularity
%
%%%%%%%%%% Multi-layer version %%%%%%%%%%
%
% [nQ] = nodalQ(A, S, omega, gamma);
%
% A - a 1xL cell array where L are the layers of the multi-layer network, and each cell contains an nxn matrix as in the form of A in the single-layer example. 
% Note that A is currently assumed to be ordered, i.e. A{1} is connected to
% A{2}, ..., A{L}. As before, please contact me if non-ordered networks are
% of interest.
%
% S - a 1x(nxL) vector where entries correspond to the group assignment.
% Note: this has only been tested using the output of the previously mentioned genLouvain
% method (using the multiord function for modularity matrix construction). 
%
% omega - inter-layer strength.
%
% gamma - intra-layer coupling parameter.
%
% Output: 
%
% nQ - a 1x(nxL) vector where sum(nQ) = Q. 
%
% Q - classical multiplex modularity.
%
%%%%%%%%%% Other %%%%%%%%%%
%
% Note on group assignment - S was intended to be calculated from modularity maximization algorithms. Be wary of using S calculated from algorithms
% which do not follow standard single/multi-layer modularity maximization.
% While potentially beneficial, nQ's interpretability and connection to standard
% modularity may change significantly. 
% 
% Note that variables in this code mirror eq. x,y. in [cite nQ paper]
% 
% If you use this code please cite as:
% [give citation]
%

%% TO DO
% include reference to nQ paper 

%% 
% initialize nQ
total_nodes = numel(S);
nQ = zeros(total_nodes,1);

% check inputs (A square and S is size of nodes A, check no self loops,
% another to make sure of this for a cell array of A, Where S now is NxL.

% check single or multi-layer 
if iscell(A)==0

    % Check if A is square
    if ~ismatrix(A) || size(A, 1) ~= size(A, 2)
        error('A must be a square matrix');
    end
    
    % Check if S and A agree
    if ~isequal(total_nodes, size(A, 1))
        error('S must be a 1xn vector, where n = #nodes');
    end
    
    % Check for self-loops
    if any(diag(A) ~= 0)
        error('A cannot contain self-loops');
    end
    
%% single-layer nQ  

	% calculate neighborhoods 
	k = sum(A,2);
	m = sum(k); % represents 2m
    
    for i=1:total_nodes

        % finds nodes of the same group
        find_group = find(S==S(i))';
        
        % k_Ni calculation
        k_Ni = k(find_group)';
        
        % k_i  calculation
        k_i = k(i);
        
        % all the edge weights between the node and its neighbors
        AiNi = A(i,find_group); 
    
        % calculate nQ
        exp_edge = (k_i*k_Ni)/(m);
        Q = (AiNi-exp_edge)/(m);
        nQ(i) = sum(Q,"all");
        
    end

else

	% # layers and nodes per layer
	L = numel(A);
	N = length(A{1});
    
    for i = 1:L
        % Check if each A is square
        if ~ismatrix(A{i}) || size(A{i}, 1) ~= size(A{i}, 2)
            error('A must contain square matrices');
        end
        
        % Check if matrix dimensions are consistent
        if size(A{i}, 1) ~= N 
            error('All matrices in A must have the same dimensions');
        end
        
        % Check for self-loops 
        if any(diag(A{i}) ~= 0)
            error('Matrices in A cannot contain self-loops)');
        end
    end
    
    % Check S dimensions for cell array case
    if numel(S) ~= L*N
        error('S must be a 1x(n*L) vector, where n is the number of nodes in A and L is the number of matrices');
    end
    
    % Check A has more than one matrix
    if L==1
        error('A is a cell with one matrix. Please input A as a double for single-layer nQ');
    end
    
%% multi-layer nQ

	% calculate u (representing 2u)
	u = 0;
	for k=1:L
		m = sum(A{k},'all');
		u = u+m; 
	end

	u = u+2*N*(L-1)*omega;

	% per layer calculation of nQ
    for j=1:L

		k = sum(A{j},2);
		m = sum(k); % represents 2m
		matrix = A{j};

		layer_index = [1+N*(j-1), N*(j)];
		S_layer = S(layer_index(1):layer_index(2));

		% upper layer
		if j<L
			S_upper = S(layer_index(2)+1:layer_index(2)+N);
		else
			S_upper = S(layer_index(1)-N:layer_index(2)-N);
		end

		% lower layer
        if j~=1 && j~=L
			S_lower = S(layer_index(1)-N:layer_index(2)-N);
        end

        for i=1:N
			% finds nodes of the same group
			find_group = find(S_layer==S_layer(i))';
		   
			% k_Ni calculation
			k_Ni = k(find_group)';
			
			% k_i calculation
			k_i = k(i);
			
			% all the edge weights between the node and its neighbors
			AiNi = matrix(i,find_group); 
			
			% inter-layer weight calculation
			if j==1 || j==L
				Cjsr=(logical(S_layer(i)==S_upper(i)*(omega)))/u;     
			else
				Cjsr=((logical(S_layer(i)==S_upper(i))+logical(S_layer(i)==S_lower(i)))*(omega))/u;
			end
		   
			% nQ calculation
			exp_edge = gamma*((k_i*k_Ni)/(m));
			Q_vector = (AiNi-exp_edge)/(u);
			
			index = i+N*(j-1);
			nQ(index) = sum(Q_vector,"all")+Cjsr;
					
        end 
    
    end 
end

% calculate modularity
Q = sum(nQ);

end