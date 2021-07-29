=================================================

  INDECS: INdependent DEComposition of networkS

=================================================

GNU Octave (https://www.gnu.org/software/octave/index) was used to develop the functions used here.



=========
Functions
=========The function indep_decomp returns a nontrivial independent decomposition of a chemical reaction network (CRN), if it exists. If no such decomposition exists, a message appears saying so. Furthermore, the output variables 'model', 'R', 'G', and 'P' allows the user to  view the following, respectively:

   - Complete network with all the species listed in the 'species' field of the structure 'model'
   - Matrix of reaction vectors of the network
   - Undirected graph of R
   - Partitions representing the decomposition of the reactions

indep_decomp uses the following functions:     1. init_graph
          - OUTPUT: Creates an empty structure that represents an undirected graph. The structure has the following fields: 'vertices' and 'edges'.
          - INPUT: none     2. add_vertex
          - OUTPUT: Adds a vertex to an undirected graph. This is indicated in the 'vertices' field of the structure representing the graph.
          - INPUTS:
                    - g: a structure with fields 'vertices' and 'edges'
                    - v: a string representing the vertex     3. add_edge
          - OUTPUT: Adds an undirected edge between two vertices. The vertex connected to a vertex is indicated in the subfield 'vertex' and the label for the edge is in the subfield 'label'. Both subfields are under the field 'edges' corresponding to the vertex. The field and subfields belong to the structure representing the graph.
          - INPUTS:
                    - g: a structure with fields 'vertices' and 'edges'
                    - v1, v2: strings representing the vertices connected by an undirected edge (make sure 'add_vertex' has been used to add the vertices 'v1' and 'v2' to g)

     4. vertex_component
          - OUTPUT: Returns a vector whose entries are the component numbers where each vertex of an undirected graph belongs to. The function returns an empty value if there are no vertices in the graph.
          - INPUT: g: a structure with fields 'vertices' and 'edges'



====
Note
====

Make sure all 5 functions are in the same folder/path being used as the current working directory.



=================================
How to fill out 'model' structure
=================================

'model' is the input for the function indep_decomp. It is a structure, representing the CRN, with the following fields (the kinetics of the network is not needed):

   - id: name of the model
   - species: a list of all species in the network; this is left blank since incorporated into the function is a step which compiles all species used in the model
   - reaction: a list of all reactions in the network, each with the following subfields:
        - id: a string representing the reaction
        - reactant: has the following further subfields:
             - species: a list of strings representing the species in the reactant complex
             - stoichiometry: a list of numbers representing the stoichiometric coefficient of each species in the reactant complex (listed in the same order of the species)
        - product: has the following further subfields:
             - species: a list of strings representing the species in the product complex
             - stoichiometry: a list of numbers representing the stoichiometric coefficient of each species in the product complex (listed in the same order of the species)
        - reversible: has the value true or false indicating if the reaction is reversible or not, respectively



========
Examples
========

7 examples are included in this folder, all based on [1]:

   - Example 1: Generalized mass action model of anaerobic fermentation pathway of Saccharomyces cerevisiae

   - Example 2: A metabolic network with one positive feedforward and a negative feedback

   - Example 3: Baccam influenza virus model

   - Example 4: Baccam influenza virus model with delayed virus production

   - Example 5: Handel influenza virus model

   - Example 6: Generalized mass action model of purine metabolism in man

   - Example 7: A reaction network governed by mass action kinetics



===================
Contact Information
===================

For questions, comments, and suggestions, feel free to contact me at pvnlubenia@yahoo.co.uk.


- Patrick Lubenia



==========
References
==========

   [1] Hernandez, B.S. and De la Cruz, R.J.L. (2021). Independent decompositions of chemical reaction networks. Bulletin of Mathematical Biology, 83(76), 1Ð23. doi:10.1007/s11538-021-00906-3.

   [2] Soranzo, N. and Altafini, C. (2009). ERNEST: a toolbox for chemical chemical reaction network theory. Bioinformatics, 25(21), 2853Ð2854. doi:10.1093/bioinformatics/btp513.