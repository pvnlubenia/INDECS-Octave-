# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                             #
#    indep_decomp                                                             #
#                                                                             #
#                                                                             #
# OUTPUT: Returns a nontrivial independent decomposition of a chemical        #
#            reaction network (CRN), if it exists. If no such decomposition   #
#            exists, a message appears saying so.                             #
#         The output variables 'model', 'R', 'G', and 'P' allow the user to   #
#            view the following, respectively:                                #
#               - Complete network with all the species listed in the         #
#                    'species' field of the structure 'model'                 #
#               - Matrix of reaction vectors of the network                   #
#               - Undirected graph of R                                       #
#               - Partitions representing the decomposition of the reactions  #
# INPUT: model: a structure, representing the CRN, with the following fields  #
#           (the kinetics of the network is not needed):                      #
#           - id: name of the model                                           #
#           - species: a list of all species in the network; this is left     #
#                blank since incorporated into the function is a step which   #
#                compiles all species used in the model                       #
#           - reaction: a list of all reactions in the network, each with the #
#                following subfields:                                         #
#                   - id: a string representing the reaction                  #
#                   - reactant: has the following further subfields:          #
#                        - species: a list of strings representing the        #
#                             species in the reactant complex                 #
#                        - stoichiometry: a list of numbers representing the  #
#                             stoichiometric coefficient of each species in   #
#                             the reactant complex (listed in the same order  #
#                             of the species)                                 #
#                   - product: has the following further subfields:           #
#                        - species: a list of strings representing the        #
#                             species in the product complex                  #
#                        - stoichiometry: a list of numbers representing the  #
#                             stoichiometric coefficient of each species in   #
#                             the product complex (listed in the same order   #
#                             of the species)                                 #
#                   - reversible: has the value true or false indicating if   #
#                        the reaction is reversible or not, respectively      #
#                                                                             #
# References:                                                                 #
#   [1] Hernandez, B.S. and De la Cruz, R.J.L. (2021). Independent            #
#          decompositions of chemical reaction networks. Bulletin of          #
#          Mathematical Biology, 83(76), 1–23. doi:10.1007/s11538-021-00906-3.#
#   [2] Soranzo, N. and Altafini, C. (2009). ERNEST: a toolbox for chemical   #
#          chemical reaction network theory. Bioinformatics, 25(21),          #
#          2853–2854. doi:10.1093/bioinformatics/btp513.                      #
#                                                                             #
# Created: 19 July 2021                                                       #
# Last Modified: 22 July 2021                                                 #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



function [model, R, G, P] = indep_decomp(model)
    
    %
    % STEP 0: Add to 'model.species' all species indicated in the reactions
    %
    
    % Get all species from reactants
    for i = 1:numel(model.reaction)
        for j = 1:numel(model.reaction(i).reactant)
            model.species{end+1} = model.reaction(i).reactant(j).species;
        end
    end
    
    % Get species from products
    for i = 1:numel(model.reaction)
        for j = 1:numel(model.reaction(i).product)
            model.species{end+1} = model.reaction(i).product(j).species;
        end
    end
    
    % Get only unique species
    model.species = unique(model.species);
    
    
    
    %
    % STEP 1: Form stoichiometric matrix N (based on [2])
    %
    
    % Count the number of species
    m = numel(model.species);
    
    % Initialize the matrix of reactant complexes
    reactant_complexes = [ ];
    
    % Initialize the matrix of product complexes
    product_complexes = [ ];
    
    % Initialize the stoichiometric matrix
    N = [ ];
    
    % For each reaction in the model
    for i = 1:numel(model.reaction)
      
        % Initialize the vector for the reaction's reactant complex
        reactant_complexes(:, end+1) = zeros(m, 1);
        
        % Fill it out with the stoichiometric coefficients of the species in the reactant complex
        for j = 1:numel(model.reaction(i).reactant)
            reactant_complexes(find(strcmp(model.reaction(i).reactant(j).species, model.species), 1), end) = model.reaction(i).reactant(j).stoichiometry;
        end
        
        % Initialize the vector for the reaction's product complex
        product_complexes(:, end+1) = zeros(m, 1);
        
        % Fill it out with the stoichiometric coefficients of the species in the product complex
        for j = 1:numel(model.reaction(i).product)
            product_complexes(find(strcmp(model.reaction(i).product(j).species, model.species), 1), end) = model.reaction(i).product(j).stoichiometry;
        end
        
        % Create a vector for the stoichiometric matrix: Difference between the two previous vectors
        N(:, end+1) = product_complexes(:, end) - reactant_complexes(:, end);
        
        % If the reaction is reversible
        if model.reaction(i).reversible
          
            % Insert a new vector for the reactant complex: make it same as the product complex
            reactant_complexes(:, end+1) = product_complexes(:, end);
            
            % Insert a new vector for the product complex: make it the same as the reactant complex
            product_complexes(:, end+1) = reactant_complexes(:, end-1);
            
            % Insert a new vector in the stoichiometric matrix: make it the additive inverse of the vector formed earlier
            N(:, end+1) = -N(:, end);
        end
    end
    
    % Count the total number of reactions
    r = size(N, 2);
    
    
    
    %
    % STEP 2: Get the transpose of N: Each row now represents the reaction vector a reaction (this point onward is based on [1])
    %
    
    R = N';
    
    
    
    %
    % STEP 3: Form a basis for the rowspace of R
    %
    
    % Write R in reduced row echelon form: the transpose of R is used so 'basis_reaction_num' will give the pivot rows of R
    %    - 'A' is R in reduced row echelon form
    %    - 'basis_reaction_num' gives the row numbers of R which form a basis for the rowspace of R
    [A, basis_reaction_num] = rref(R');
    
    % Form the basis
    basis = R(basis_reaction_num, :);
    
    
    
    %
    % STEP 4: Construct the vertex set of undirected graph G
    %
    
    % Initialize an undirected graph G
    G = init_graph();
    
    % Add vertices to G: these are the reaction vectors that form a basis for the rowspace of R
    for i = 1:numel(basis_reaction_num)
        
        % Use the reaction number as label for each vertex
        G = add_vertex(G, strcat('R', num2str(basis_reaction_num(i))));
    end
    
    
    
    %
    % STEP 5: Write the nonbasis reaction vectors as a linear combination of the basis vectors
    %
    
    % Initialize matrix of linear combinations
    linear_combo = zeros(r, numel(basis_reaction_num));
    
    % Do this for the nonbasis reactions vectors
    for i = 1:r
        if ~ismember(i, basis_reaction_num)
          
          % This gives the coefficients of the linear combinations
          % The basis vectors will have a row of zeros
          linear_combo(i, :) = basis'\R(i, :)';
        end
    end
    
    % Round off values to nearest whole number to avoid round off errors
    linear_combo = round(linear_combo);
    
    
    
    %
    % STEP 6: Construct the edge set of undirected graph G
    %
    
    % Get the reactions that are linear combinations of at least 2 basis reactions
    % These are the reactions where we'll get the edges
    get_edges = find(sum(abs(linear_combo), 2) > 1);
        
    % Initialize an array for sets of vertices that will form the edges
    vertex_set = { };
     
    % Identify which vertices form edges in each reaction: get those with non-zero coefficients in the linear combinations
    for i = 1:numel(get_edges)
        vertex_set{i} = find(linear_combo(get_edges(i), :) ~= 0);
    end
    
    % Initialize the edge set
    edges = [ ];
    
    % Get all possible combinations (not permutations) of the reactions involved in the linear combinations
    for i = 1:numel(vertex_set)
        edges = [edges; nchoosek(vertex_set{i}, 2)];
    end
    
    % Get just the unique edges
    edges = unique(edges, 'rows');
    
    % Add these edges to graph G
    for i = 1:size(edges, 1)
        G = add_edge(G, strcat('R', num2str(basis_reaction_num(edges(i, 1)))), strcat('R', num2str(basis_reaction_num(edges(i, 2)))));
    end
    
    
    
    %
    % STEP 7: Check if G is connected, i.e., has only one connected component
    %
    
    % Determine to which component each vertex belongs to
    component_numbers = vertex_component(G);
    
    % Determine the number of connected components of G: this is the number of partitions R will be decomposed to
    num_components = max(component_numbers);
    
    % For the case of only one connected component
    if num_components == 1
        P = [];
        disp([model.id ' has no nontrivial independent decomposition.']);
        
        % 'return' exits the function; we don't need to continue the code
        % If we wanted to just get out of the loop, we use 'break'
        return
    end
    
    
    
    %
    % STEP 8: If G is NOT connected, form the partitions of R
    %
    
    % Initialize the list of partitions
    P = cell(1, num_components);
    
    % Basis vectors: assign them first into their respective partition based on their component number
    for i = 1:numel(component_numbers)
        P{component_numbers(i)}(end+1) = basis_reaction_num(i);
    end
    
    % Nonbasis vectors: they go to the same partition as the basis vectors that form their linear combination
    for i = 1:numel(P)
        for j = 1:numel(P{i})
            
            % Get the column number representing the basis vectors in 'linear_combo'
            col = find(basis_reaction_num == P{i}(j));
            
            % Check which reactions used a particular basis vector and assign them to their respective partition
            P{i} = [P{i} find(linear_combo(:, col) ~= 0)'];
        end
    end
    
    % Get only unique elements in each partition
    for i = 1:numel(P)
        P{i} = unique(P{i});
    end
    
    
    
    %
    % FINAL STEP: Display the independent decomposition
    %
    
    % Use 'fprintf' instead of 'disp' to interpret '\n' as 'newline'
    fprintf(['Independent decomposition of ' model.id '\n\n'])
    for i = 1:numel(P)
        fprintf(['P' num2str(i) ': ' sprintf('R%d ', P{i}) '\n\n']);
    end