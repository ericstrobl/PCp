function [pdag,cell_p2,num_struc2,G]=PC_original(G, sep, cell_p, num_struc, k, cond_indep, varargin)

[pdag, cell_p2, num_struc2] = get_v_structures2_original(G, sep, cell_p, num_struc, k, cond_indep, varargin{:}); % get v-structures, unconstrained edge directions

[pdag, cell_p2, num_struc2] = orientation_rules_original(G, pdag, cell_p2, num_struc2); %unconstrianed edge propagation via orientation rules