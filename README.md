# epistasis_simulator

Currently working on
-----------------------------
- Creating the classes for the primogenitor and progeny
- Creating methods for the primogenitor to create a progeny (following can be done in any order contrary to how I introduce the problems :,) )
  - First thing is to transfer the functions for deep mutational analysis over from the Jupyter Notebook, a problem for mutliple screens
  - Next would be to create a system where the progeny can be identified by its sequence of mutations
    - Currently I think that the best way to do this is to just make a list of tuples, with the first index being the original acid, second being the mutation, and third
    being the location.
      - [(O, M, X)] O being original, M being mutation, X being the index
      - Maybe, just maybe make a method in the progeny that allows for recreation of the parent strand, couldn't be that hard :D
      
