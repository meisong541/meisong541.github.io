---
layout: post
title:  "Statistical physics and random matrix theory: Replica method I"
date:   2019-07-20 00:09:06 -0700
categories: jekyll update
---



$$\def\bsigma{\boldsymbol \sigma}$$
$$\def\bv{\boldsymbol v}$$
$$\def\bu{\boldsymbol u}$$
$$\def\sT{\mathsf T}$$
$$\def\bW{\boldsymbol W}$$
$$\def\bA{\boldsymbol A}$$
$$\def\R{\mathbb R}$$
$$\def\S{\mathbb S}$$
$$\def\GOE{\text{GOE}}$$
$$\def\|{\Vert}$$
$$\def\bx{\boldsymbol x}$$
$$\def\cN{\mathcal N}$$
$$\def\E{\mathbb E}$$
$$\def\de{\text{d}}$$
$$\def\vphi{\varphi}$$
$$\def\bQ{\boldsymbol Q}$$
$$\def\diag{\text{diag}}$$
$$\def\bzero{\boldsymbol 0}$$
$$\def\id{\mathbf I}$$
$$\def\ones{\mathbf 1}$$
$$\def\ext{\text{ext}}$$


## 1. Introduction

In my research, I encountered many interesting problems that is related to random matrices. There are many rigorous and powerful tools dealing with random matrices, like concentration inequalities, Stieltjes transform, leave-one-out method, etc. Besides these methods, there is one tool that I find it to be extremely powerful. It is called the replica method, which is a non-rigorous method originated from statistical physics. I used several times this tool in my research, and it helped me a lot to "guess" the answer quickly. I gave two talks about this method, respectively in Andrea's and Tengyu's group meetings. Now I think it is the time to summarize my talks and present it to the large audiences.


Replica method was used to calculate the free energy associated to some random Hamiltonians in statistical physics. Its first well-known application is calculating the free energy of Sherrington-Kirkpatrick model. Beyond the physicists' community, replica method is a very powerful tool studying random matrices. Its applications includes coding theory, compressed sensing, and high dimensional statistical inference problems. In this blog, I will discuss about how to use replica method to predict the spectral norm of spiked GOE matrix.

For the references for replica method, I suggest people reading: .....



## 2. A motivating example: spiked GOE matrix

### 2.1. The spiked GOE model

We consider a symmetric random matrix $$\bA \in \R^{n \times n}$$ with
\\[
\bA= \lambda \bu \bu^\sT + \bW,
\\]
where $$\lambda \in \R_+$$ is the signal to noise ratio, $$\bu \in \S^{n-1} = \{ \bx \in \R^n: \| \bx \|_2 = 1 \}$$ is a spike vector, and $$\bW \sim \GOE(n)$$: that means, $$\bW$$ is a symmetric matrix, with $$W_{ij} \sim \cN(0, 1/n)$$ for $$1 \le i < j \le n$$ and $$W_{ii} \sim \cN(0, 2/n)$$ for $$1 \le i \le n$$. 


### 2.2. The BBP phase transitions

The following properties of the spiked GOE matrix is called BBP phase
transition.

-   The largest eigenvalue of $$\bA$$: <br/>
    For $$0 \le \lambda \le 1$$, we have
    \\[
    \lim_{n \to \infty} \E[\lambda_{\max} ({\bA})] = 2;
    \\]
    For $$\lambda > 1$$, we have
    \\[
    \lim_{n \to \infty} \E[\lambda_{\max} ({\bA})] = \lambda + 1/\lambda.
    \\]

-   The correlation of $$\bv_{\max}(\bA)$$ with $$\bu$$: <br/>
    For $$0 \le \lambda \le 1$$, we have
    \\[
    \lim_{n \to \infty} \E[\langle \bv_{\max}(\bA), \bu \rangle^2] = 0;
    \\]
    For $$\lambda > 1$$, we have
    \\[
    \lim_{n \to \infty} \mathbb{E}[\langle {\boldsymbol v}_{\max}({\boldsymbol A}), {\boldsymbol u}\rangle^2] = 1 - 1/\lambda^2.
    \\]

In the following section, we will discuss how to derive this result using replica method. 
    

## 3. Tricks in statistical physics


### 3.1. The free energy trick


Suppose $$H: \Omega \to \R$$ and $$f: \Omega \to \R$$ are two functions. We would like to analytically calculate the following quantities
\\[
\begin{aligned}
&\max_{\bsigma \in \Omega} H(\bsigma),\\
&f(\arg \max_{\bsigma} H(\bsigma)).
\end{aligned}
\\]
The free energy trick is the following. Define a Gibbs measure
\\[
\mu_{\beta, h}(\de \bsigma) = \frac{1}{Z(\beta, h)} \exp\\{ \beta H(\bsigma) + h f(\bsigma)\\} \nu_0(\de \bsigma),
\\]
where $$Z(\beta, h)$$ is the normalizing constant
\\[
Z(\beta, h) = \int_{\Omega} \exp\\{ \beta H(\bsigma) + h f(\bsigma)\\} \nu_0(\de \bsigma). 
\\]
In physics, $$\beta$$ stands for inverse temperature, and $$h$$ stands for the strength of external field. Define the free energy function $$\Psi$$ to be 
\\[
\Psi(\beta, h) = \log Z(\beta, h). 
\\]
The following lemma shows that, as long as we can derive a formula for the free energy function $\Psi$, we can derive many interesting quantities including $$\max_{\bsigma \in \Omega} H(\bsigma)$$ and $$f(\arg \max_{\bsigma} H(\bsigma))$$. This lemma follows by basic calculus, which is the same as the calculus for exponential family. 

<div class="lemma">
Denote the ensemble average operator <script type="math/tex">\langle \cdot \rangle</script> to be the following:
\[
\langle g \rangle_{\beta, h} = \int_\Omega g(\bsigma) \mu_{\beta, h}(\de \bsigma).
\]
Then we have
\begin{align}\tag{1}\label{eq:1}
\partial_\beta \Psi(\beta, h) =& \langle H \rangle_{\beta, h},\\
\partial_h \Psi(\beta, h) =& \langle f \rangle_{\beta, h},\\
\end{align}
and
\begin{align}\tag{2}\label{eq:2}
\lim_{\beta \to \infty} \partial_\beta \Psi(\beta, h) =& \max_{\bsigma \in \Omega} H(\bsigma),\\
\lim_{\beta \to \infty} \partial_h \Psi(\beta, h) =& f(\arg\max_{\bsigma \in \Omega} H(\bsigma) ).
\end{align}
</div>
Equations \eqref{eq:1} follow from basic calculus of exponential family, and Equations \eqref{eq:2} follow from the fact that, for $$h = 0$$ and large $$\beta$$, the measure $$\mu_{\beta, 0}$$ will concentrate at the maxima of $$H$$. 



### 3.2. Free energy trick for spiked GOE model

In the spiked GOE model, the Hamiltonian associated with the random matrix $$\bA$$ is defined as
\\[
H_{n, \lambda}(\bsigma) = \langle \bsigma, \bA \bsigma \rangle = \lambda \langle \bu, \bsigma\rangle^2 + \langle \bsigma, \bW \bsigma\rangle. 
\\]
Denote the (random) partition function associated with the Hamiltonian $$H_{n, \lambda}$$:
\\[
Z(\beta, \lambda, n) = \int_{\S^{n-1}} \exp\\{ \beta n H_{n, \lambda}( \bsigma ) \\} \nu_0(\de \bsigma).
\\]
Here $$\nu_0$$ is the uniform probability measure on $$\S^{n-1}$$, and the $$n$$ factor on the exponent serves for proper normalization. 

Define 
\\[
\begin{align}\tag{3}\label{eq:3}
\vphi_n(\lambda) =& \lim_{\beta \to \infty} \frac{1}{\beta n}\E[\log Z(\beta, \lambda, n)],\\\
\vphi(\lambda) =& \lim_{n\to \infty} \vphi_n(\lambda) = \lim_{n \to \infty} \lim_{\beta \to \infty} \frac{1}{\beta n}\E[\log Z(\beta, \lambda, n)].
\end{align}
\\]
The following lemma tells us that, in order to get the asymptotics of $$\E[\lambda_{\max}(\bA)]$$ and $$\E[\langle \bv_{\max}(\bA), \bu\rangle^2]$$, we just need to calculate $$\vphi(\lambda)$$. 

<div class="lemma">
For <script type="math/tex">\bA = \lambda \bu \bu^\sT + \bW \in \R^{n \times n}</script>, we have
\[
\begin{aligned}
\vphi_n(\lambda) =& \E[\sup_{\bsigma} H_{n, \lambda}(\bsigma)]= \E[\lambda_{\max} (\bA)],\\
\vphi_n'(\lambda) =& \E[\langle \bv_{\max}(\bA), \bu \rangle^2]. 
\end{aligned}
\]
</div>
We leave the proof of this Lemma to the readers. 




### 3.3. The replica trick

To calculate the function $$\vphi$$, we need to calculate $$\E[\log Z]$$. However, it is hard to calculate this quantity, since we don't know how to deal with the expectation of $$\log$$. Instead, there is an easy way to calculate the moments of $$Z$$. The replica trick is to convert the problem of calculating the expectation of $$\log
$$ into the moments calculations. 

The following formula is called replica trick, 
\\[\tag{4} \label{eq:4}
\E[\log Z] = \lim_{k \to 0} \frac{1}{k}\log \E[Z^k ], 
\\]
with justification 
\\[
\begin{aligned}
\E[\log Z] =& \E[(\log Z^k) / k] = \lim_{k \to 0} \E [ \log(1 + Z^k - 1) / k] = \lim_{k \to 0} \E [ ( Z^k - 1 ) / k ] \\\
=&\lim_{k \to 0} (\E [ Z^k] - 1)/k = \lim_{k \to 0} \frac{1}{k} \log (1 + \E[Z^k] - 1) = \lim_{k \to 0} \frac{1}{k}\log \E[Z^k ]. 
\end{aligned}
\\]

### 3.4. The replica trick for spiked GOE model

Using the definition of $$\vphi$$ in Eq. \eqref{eq:3} and the replica trick Eq. \eqref{eq:4}, we can rewrite $$\vphi$$ to be
\\[
\begin{aligned}
\vphi(\lambda) \equiv& \lim_{n\to \infty} \lim_{\beta \to \infty} \frac{1}{\beta n}\E[\log Z(\beta, \lambda, n)] = \lim_{n \to \infty}   \lim_{\beta \to \infty} \lim_{k \to 0}\frac{1}{\beta k n}\log \E[Z(\beta, \lambda, n)^k]. \\\
\end{aligned}
\\]
With some heuristic change of limit (which we will not justify), we get
\\[
\vphi(\lambda) \stackrel{\cdot}{=} \lim_{\beta \to \infty}  \lim_{k \to 0} \lim_{n \to \infty} \frac{1}{\beta k n}\log \E[Z(\beta, \lambda, n)^k]. 
\\]
There are three limits in the right hand side. Define
\\[
S(k, \beta, \lambda) = \lim_{n \to \infty} \frac{1}{n}\log \E[Z(\beta, \lambda, n)^k],
\\]
and
\\[
\Psi(\beta, \lambda) = \lim_{k \to 0} \frac{1}{k} S(k, \beta, \lambda).
\\]
This gives
\\[
\vphi(\lambda) = \lim_{\beta \to \infty} \frac{1}{\beta}\Psi(\beta, \lambda). 
\\]
In the following section, we calculate these $$S$$, $$\Psi$$, and $$\vphi$$ functions sequentially. 




## 4. The replica calculations

### 4.1. The $$n \to \infty$$ limit: moments calculation

#### 4.1.1. The result for $$S$$ function

We claim the following expression for function $$S$$
\\[
S(k, \beta, \lambda) \equiv \lim_{n \to \infty} \frac{1}{n}\log\E[Z(\beta, \lambda, n)^k]= \sup_{\bQ \in \R^{(k+1) \times (k+1)}, \diag(\bQ) = 1, \bQ \succeq \bzero} \{U(\bQ) \}, 
\\]
where
\\[
U(\bQ) = \beta \lambda \sum_{i=1}^k q_{0i}^2 + \beta^2 \sum_{ij = 1}^k q_{ij}^2 + \frac{1}{2} \log(\det(\bQ))
\\]
and $$\bQ \in \R^{(k + 1) \times (k + 1)}$$ is symmetric, with (the index of $$\bQ$$ start from $$0$$ to $$k$$)
\\[
\bQ = \begin{bmatrix}
1& q_{01} & \ldots & q_{0k}\\\
q_{01} & 1 & \ldots & q_{1k}\\\
\ldots & \ldots & \ldots & \ldots\\\
q_{0k} & q_{k1} & \ldots & 1\\\
\end{bmatrix}.
\\]
This result is exact and rigorous, but the computation is involved. To derive this result, we calculate the moments of $$Z$$ directly. After a long calculations, we can get
\\[
\E[Z(\beta, \lambda, n)^k] = \int \exp\\{ n U(\bQ) + o(n)\\} \de \bQ, 
\\]
and applying Laplace method. We will discuss about how to do this calculation at the end of the talk. First let us look at what is Physicists' philosophy after getting to this point. 





### 4.2. The $$k \to 0$$ limit: replica symmetric ansatz

#### 4.2.1. Replica symmetric ansatz

We would like to calculate $$\Psi(\beta, \lambda) = \lim_{k \to 0} S(\beta, \lambda, k)/k$$. However, our expression for $$S$$ is only valid for integer $$k$$. Hence, this limit taken here is completely non-rigorous. 

Note $$k$$ serves as the dimension of variable $$\bQ$$ in the expression of $$S$$, so we cannot take the limit. The trick here is to plug in an ansatz for matrix $$\bQ$$ to simplify the expression for $$S$$. Note the maximization problem in the expression of $$S$$ has some symmetric properties, this motivate us to take the "replica symmetric ansatz": 
\\[
\bQ = \begin{bmatrix}
1& \mu & \ldots & \mu\\\
\mu & 1 & \ldots & q\\\
\ldots & \ldots & \ldots & \ldots\\\
\mu & q & \ldots & 1\\\
\end{bmatrix}.
\\]
Under this ansatz, using matrix determinant lemma, we have
\\[
\begin{aligned}
&\log(\det(\bQ)) = \log \det( (1-q) \id_k + (q - \mu^2) \ones \ones^\sT) \\\
=&  \log \det( \id_k + [(q - \mu^2)/(1-q) ] \ones \ones^\sT) + k \log(1-q)\\\
=& \log(1 + k [(q - \mu^2)/(1-q) ]) + k \log(1-q)\\\
=& \log\Big( 1 - \frac{\mu^2 k}{1 + (k - 1)q} \Big) + \log \Big(1 + \frac{kq}{1-q} \Big) + k \log(1 - q),
\end{aligned}
\\] 
and this gives
\\[
S(k, \beta, \lambda) \stackrel{\cdot}{=} \sup_{\mu, q} U(\mu, q, k).
\\]
where
\\[
U(\mu, q, k) = \beta \lambda k \mu^2 + \beta^2 k + \beta^2 k (k-1) q^2 + \frac{1}{2} \Big[ \log\Big( 1 - \frac{\mu^2 k}{1 + (k - 1)q} \Big) + \log \Big(1 + \frac{kq}{1-q} \Big) + k \log(1 - q) \Big].
\\]


#### 4.2.2. The $$k \to 0$$ limit
Now we would like to take the limit $$k \to 0$$ using the expression above. Define 
\\[
u(\mu, q) = \lim_{k \to 0} \frac{1}{k} U(\mu, q, k) = \beta \lambda \mu^2  + \beta^2 (1 - q^2) - \frac{\mu^2}{2(1 - q)} + \frac{q}{2 (1 - q)} + \frac{1}{2}\log(1 - q). 
\\]
Here the convention is that, the $$\lim_{k \to 0} \sup_{\mu, q}$$ becomes $$\exp_{\mu, q} \lim_{k \to 0}$$: 
\\[
\Psi(\beta, \lambda) =  \lim_{k \to 0} S(k, \beta, \lambda) / k \stackrel{\cdot}{=} \ext_{\mu, q}u(\mu, q),
\\]
where $$\ext_{\mu, q}$$ means setting stationery (extrema) with respect to $$\mu$$ and $$q$$. 

We have 
\\[
\begin{aligned}
\partial_\mu u(\mu, q) =& 2 \{\beta \lambda - 1/[2 (1-q)] \} \mu, \\\
\partial_q u(\mu, q) =& -2 \beta^2 q -  \frac{\mu^2}{2(1 - q)^2} + \frac{1}{2 (1 - q)^2}.
\end{aligned}
\\]
When $$\beta > 1$$, there are two branches of solutions of the extrema of $$u$$: $$\mu = 0$$ and $$\mu \neq 0$$. 



#### 4.2.3. The first branch $$\mu = 0$$

When $$\mu = 0$$, 
\\[
\Psi_1(\beta, \lambda) = \ext_{q}\Big\\{ \beta^2 (1 - q^2) + \frac{q}{2 (1 - q)} + \frac{1}{2}\log(1 - q)\Big\\}. 
\\]
the extrema is taken at $$q = 1 - 1/(2 \beta)$$, hence we have 
\\[
\Psi_1(\beta, \lambda) =2 \beta - 3/4 - \frac{1}{2} \log(2 \beta).
\\]

#### 4.2.4. The second branch $$\mu \neq 0$$

When $$\mu \neq 0$$, we must have 
\\[
\begin{aligned}
\mu =& \Big( (1 - 1 / \lambda^2) (1 - 1/(2 \beta \lambda))  \Big)^{1/2},\\\
q =& 1 - 1/(2 \lambda \beta),
\end{aligned}
\\]
and hence
\\[
\Psi_2(\beta, \lambda) =  \beta( \lambda + 1/ \lambda) - [1/(4 \lambda^2)+ 1/2] - \frac{1}{2}\log(2 \lambda \beta). 
\\]




#### 4.3. The $$\beta$$ limit

In the first branch, we have 
\\[
\vphi_1(\lambda) = \lim_{\beta \to \infty} \frac{1}{\beta} \Psi_1(\beta, \lambda) = 2. 
\\]
In the second branch, we have 
\\[
\vphi_2(\lambda) = \lim_{\beta \to \infty} \frac{1}{\beta} \Psi_2(\beta, \lambda) = \lambda + 1/\lambda. 
\\]

#### 4.4. Select a branch

We would like to select a branch that makes sense: $$\vphi$$ should be non-decreasing and $$\lim_{\lambda \to \infty} \vphi(\lambda) = \infty$$. 
- When $$\lambda\le 1$$, the branch $$\vphi_2(\lambda) = \lambda + 1/\lambda$$ is decreasing, and hence the branch $$\vphi_1(\lambda) = 2$$ makes sense. 
- When $$\lambda \ge 1$$, the branch $$\vphi_1(\lambda) = 2$$ stays constant, hence the branch $$\vphi_2(\lambda) = \lambda + 1/\lambda$$ makes sense. 


Combining the above discussions, we predict that 
\\[
\vphi(\lambda) = \begin{cases}
2, ~~~& \lambda \le 1, \\\
\lambda + 1/\lambda, ~~~& \lambda > 1.
\end{cases}
\\]
This turned out to be the correct answer! 



## 5. Summary

In summary, to calculate the extreme eigenvalue of the spiked GOE matrix, we used free energy trick and replica trick to transform the problem into sequentially calculating three limits: the $$n \to \infty$$ limit, the $$k \to 0$$ limit, and the $$\beta \to \infty$$ limit. The $$n \to \infty$$ limit can be dealt with rigorously, and there are some principled way to carry out the calculations. However, the longest calculations happen here. The $$k \to 0$$ limit is non-rigorous and no methods till now can make it rigorous. Some mysterious things happen in this part, where after exchanged with the limit, the $$\sup$$ operation becomes an $$\ext$$ operation. The final $$\beta \to \infty$$ limit is easy to deal with. After sequentially taking these limit, one need to use some intuitions to select a solution that is physical. 









