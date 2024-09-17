# misclassifyr

`misclassifyr` is an R package for estimating misclassification models, as discussed in Mattheis (2024). 

Given a discrete regressor $X$ and two noisy measures $Y_1, Y_2$  of a discrete outcome $Y^*$, this package presents a menu of options for estimation and inference on:

- the joint distribution of $X$ and $Y^*$, denoted in the matrix $\Pi$. 
- the distribution of measurement error $(Y_1,Y_2) \hspace{0.25em} | \hspace{0.25em} Y^*$ . denoted by the matrix $\Delta$.
- a linear regression coefficient $\beta$  corresponding to the projection of $Y^*$ on $X$, assuming some mapping from $X$ and $Y^*$ to scalars. 

  

---

To install, run: `devtools::install_github("ramattheis/misclassifyr")`

*note:* This package is still at an early stage of development. No promises are made about backwards compatibility at this point. 

