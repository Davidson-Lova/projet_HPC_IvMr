# projet_HPC_IvMr

Projet pour HPC calcul parallèle sur la parallélisation d'un code pour la résolution d'un problème de convection/diffusion de la chaleur en 2 dimensions, pour $(x,y) \in [0,10] \times [0,10]$ : 

$$\frac{\partial T}{\partial t} + U \frac{\partial T}{\partial x} + V \frac{\partial T}{\partial y} = k \left(\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} \right)$$

Le sujet complet est [ici](https://github.com/Davidson-Lova/projet_HPC_IvMr/blob/master/projet_cours_parallele.pdf).

Ce repertoire comporte les quatres fichiers demander :
- [Advection_diffusion_seq.cpp](https://github.com/Davidson-Lova/projet_HPC_IvMr/blob/master/Advection_diffusion_seq.cpp)
- [Advection_diffusion_omp.cpp](https://github.com/Davidson-Lova/projet_HPC_IvMr/blob/master/Advection_diffusion_omp.cpp)
- [Advection_diffusion_mpiv1.cpp](https://github.com/Davidson-Lova/projet_HPC_IvMr/blob/master/Advection_diffusion_mpiv1.cpp)

qui fonctionnent comme leur nom l'indique.

Le fichier suivant :
- [Advection_diffusion_mpiv2.cpp](https://github.com/Davidson-Lova/projet_HPC_IvMr/blob/master/Advection_diffusion_mpiv2.cpp)
doit être utiliser avec [res_norm.py](https://github.com/Davidson-Lova/projet_HPC_IvMr/blob/master/res_norm.py) pour renormaliser le format du res.txt au bon format pour gnuplot.

Ce [Notebook](https://github.com/Davidson-Lova/projet_HPC_IvMr/blob/master/note.ipynb) est utiliser pour tester l'integralité des fonctions demander.

Le Notebook [sandbox.ipynb](https://github.com/Davidson-Lova/projet_HPC_IvMr/blob/master/sandbox.ipynb) est un brouillon.