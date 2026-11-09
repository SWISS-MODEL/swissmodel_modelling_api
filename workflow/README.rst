=======================
SWISS-MODEL to FoldTree
=======================

A Galaxy workflow that builds
`SWISS-MODEL <https://swissmodel.expasy.org/>`_ homology models for all
sequences of a multi-FASTA file, keeps the models with a GMQE of at least a
given threshold and builds structure-based trees from them with
`FoldTree <https://doi.org/10.1038/s41594-025-01649-8>`_.

Each sequence is modelled with SWISS-MODEL automodel through the SWISS-MODEL
API, so the runtime depends on the SWISS-MODEL server.

Inputs
------

- **Multi-FASTA file**: protein sequences to model, one record per target.
  Headers name the results and the tree leaves, so keep them short and
  unique.
- **GMQE threshold**: models with a lower GMQE are left out of the trees
  (default 0.7).
- **Project title** (optional): title of the SWISS-MODEL projects, also used
  to name the model collections.

Outputs
-------

- **FoldTree metric tree**, **TM-score tree**, **LDDT tree**: Newick trees of
  the models passing the filter. The FoldTree metric tree is the main result,
  the other two are for comparison.
- **Models passing GMQE filter**: table of sequence, model and GMQE of the
  models in the trees.
- **Model PDB files** and **Model metadata**: all SWISS-MODEL models, grouped
  by sequence, including the ones that did not pass the filter.

FoldTree needs at least two structures; if fewer models pass the filter, no
tree is built.

Requirements
------------

- A SWISS-MODEL API token, from your account on
  https://swissmodel.expasy.org. Galaxy asks for it as
  "SWISS-MODEL API Credentials" and stores it in your user preferences.
- The Galaxy server needs these Tool Shed tools installed:
  ``iuc/swissmodel_modelling_api``, ``iuc/foldtree``, ``iuc/jq``,
  ``iuc/compose_text_param``, ``rnateam/splitfasta``,
  ``nml/collapse_collections`` and ``bgruening/text_processing``.

Importing the workflow
----------------------

**For yourself:** in Galaxy, open *Workflows*, click *Import* and upload
`SWISS-MODEL_to_FoldTree.ga <SWISS-MODEL_to_FoldTree.ga>`_. The workflow then
shows up under *My workflows* and can be run from there.

**For other users of the same Galaxy server:** in *My workflows*, open the
workflow's sharing settings and make it accessible or publish it. Published
workflows are listed under *Public workflows* for everyone on that server.

**For all Galaxy servers:** submit the workflow to the
`Intergalactic Workflow Commission (IWC) <https://github.com/galaxyproject/iwc>`_.
IWC workflows are tested, listed on
`WorkflowHub <https://workflowhub.eu>`_ and
`Dockstore <https://dockstore.org>`_, and installed on the public
usegalaxy.* servers.

Attribution
-----------

When you publish results of this workflow, please cite
`SWISS-MODEL <https://doi.org/10.1093/nar/gky427>`_ and
`FoldTree <https://doi.org/10.1038/s41594-025-01649-8>`_.
