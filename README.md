Channel Capacity of T cells  
==============================

This project aims to compute the channel capacity of the T cell signaling network within the context of 
kinetic proofreading. Calculations show how proofreading increases T cell channel capacity until 1 bit of information is
extracted from the environment.

All chemical reactions are solved using a system of ordinary differential equations simulated using the PySB 
(Systems biology modeling in Python) package. Reactions can also be solved using stochastic simulation algorithms 
(Gillespie).


Data
----
All processed data is located in /data/processed/. Plots shown in the main text can be reproduced using the jupyter notebooks
in the respective folders.

Data and plots for histograms shown in Figure 1C is located in:  
data/processed/LATP_Extension/latpp_k_on_2/1_step  
data/processed/LATP_Extension/latpp_k_on_2/4_step  
data/processed/LATP_Extension/latpp_k_on_2/7_step  

Data and plotting for Figure 2 is located in:    
data/processed/LATP_Extension/latpp_k_on_2/   
Slow_Early_Step_Bound_LATPP.ipynb   

Data and plots for Figure 5 is located in:    
data/processed/New_Negative_Fb_Param_Search/  
Neg_Fb_Parameter_Search.ipynb  
data/processed/Late_Negative_Fb/  
data/processed/Positive_Fb_Loop_New_Model/

Questions? Email me at rganti(at)mit.edu

Project Organization
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── data           <- Simulation scripts to generate data
    │   │   └── __init__.py
    │   │   ├── negative_feedback_search.py     <- Scripts for doing parameter search with negative feedback loops
    │   │   ├── pysb_t_cell_network.py     <- Scripts for running tcr kinetic proofreading network with pysb
        
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │   └── build_features.py
    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │   │                 predictions
    │   │   ├── predict_model.py
    │   │   └── train_model.py
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │       └── visualize.py
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io


--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
