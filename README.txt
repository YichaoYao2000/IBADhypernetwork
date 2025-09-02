# Higher-Order Dynamic Network Modeling with IBAD

This repository contains the core simulation codes used in the manuscript:

**"The Inheritance of Connectivity: Modeling Higher-Order Dynamic Networks with Birth-Death Process"**

The codes implement:
1. A **higher-order dynamic network model** governed by an inheritance-based deletion (IBAD) mechanism, including continuous-time birth–death processes, preferential attachment, and inheritance rules.  
2. **Empirical validation** by modeling real-world networks such as food webs, plant–pollinator systems, and email communication networks, with comparisons of degree distributions and topological measures.

---

## File Descriptions

### `higher_order_network_modeling.py`
- Simulates the dynamic evolution of a higher-order network (hypergraph) with IBAD.  
- Implements node/hyperedge birth and deletion with inheritance-based redistribution.  
- Computes structural measures: degree distributions, clustering coefficient, average path length, PLCC, and triadic interactions.  
- Produces simulation results that verify theoretical stationarity and scale-free properties.

### `modeling_real_network.py`
- Loads empirical network data (food webs, plant–pollinator networks, email communication networks).  
- Builds hypergraphs from data and compares their properties with IBAD simulations.  
- Provides degree distribution fitting (power-law) and correlation analysis (Spearman’s ρ).

---

## Requirements
- Python 3.x  
- Packages: `networkx`, `numpy`, `scipy`, `matplotlib`, `pandas`  

---

## Usage

### 1. Simulated Networks
```bash
python higher_order_network_modeling.py
```
Generates simulations of IBAD-based higher-order networks and outputs degree distributions, clustering, PLCC, etc.

### 2. Empirical Networks
Place network files (e.g., `FW_008.csv`, `M_PL_015.csv`, `email-Eu-full-nverts.txt`, `email-Eu-full-simplices.txt`) in the working directory and run:
```bash
python modeling_real_network.py
```
Builds real hypergraphs, fits degree distributions, and computes similarity metrics.

---

## Notes
- The IBAD mechanism uses degree-history correlations to reassign connections after deletions.  
- Parameters (birth/death rates, attachment numbers) can be adjusted inside the scripts.  
- Codes are provided to reproduce the results in the manuscript and can be extended to other datasets.

---

## Contact
For questions, please contact:  
- **Yichao Yao** – *yichaoyao2000@gmail.com*  
