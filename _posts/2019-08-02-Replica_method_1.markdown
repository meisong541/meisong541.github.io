---
layout: post
title:  "Replica method and random matrices (I)"
date:   2019-08-04 00:09:06 -0700
categories: jekyll update
---

$$\def\blambda{\boldsymbol \lambda}$$
$$\def\bX{\boldsymbol X}$$
$$\def\N{\mathbb N}$$
$$\def\P{\mathbb P}$$
$$\def\bi{\boldsymbol i}$$
$$\def\bG{\boldsymbol G}$$
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
$$\def\|{\Vert}$$
$$\def\bLambda{\boldsymbol \Lambda}$$
$$\def\const{\text{const}}$$
$$\def\Unif{\text{Unif}}$$
$$\def\bSigma{\boldsymbol \Sigma}$$



## 1. Introduction


In my first several posts, I will discuss about the replica method. The replica method is a non-rigorous but powerful approach that originated within the statistical physics literature. It has been widely applied in other fields, including coding theory and high dimensional statistics. I used this method several times in my research. It helped me quickly figure out the answer of numerous problems. 


Most physicists know the replica method better than me. Here I would like to introduce the replica method to theorists outside the physics community. I will not assume knowledges of statistical physics from the audience, but I will occasionally use terminology from physics. For a more sophisticated introduction to the replica method, I recommend Chapter 8 of this <a href = "https://web.stanford.edu/~montanar/RESEARCH/book.html">fantastic book</a> [[1]](#ref1). 

 
In this post, I will discuss how to use the replica method to calculate the spectral norm of certain random matrices. In the next few posts, I will discuss how to use it to calculate the spectral density (Stieltjes transforms) of certain random matrices. 



## 2. A motivating example: spiked GOE matrix

We consider a symmetric random matrix $$\bA \in \R^{n \times n}$$, 
\\[
\bA= \lambda \bu \bu^\sT + \bW,
\\]
where $$\lambda \ge 0$$ is the signal to noise ratio, $$\bu \in \S^{n-1} \equiv \{ \bx \in \R^n: \| \bx \|_2 = 1 \}$$ is a spike vector, and $$\bW \sim \GOE(n)$$: meaning that $$\bW \in \R^{n \times n}$$ is a symmetric matrix, $$W_{ij} \sim \cN(0, 1/n)$$ for $$1 \le i < j \le n$$, and $$W_{ii} \sim \cN(0, 2/n)$$ for $$1 \le i \le n$$. 

We are interested in the largest eigenvalue and the corresponding eigenvector of $$\bA$$, which are respectively denoted by $$\lambda_{\max}(\bA)$$ and $$\bv_{\max}(\bA)$$. The result below is called the BBP phase transition phenomenon (see [[2]](#ref2) for spiked Wishart matrix and [[3]](#ref3) for spiked Wigner matrix). 

-   The largest eigenvalue $$\lambda_{\max}(\bA)$$: <br/>
    For $$0 \le \lambda \le 1$$, we have
    \\[
    \lim_{n \to \infty} \E[\lambda_{\max} ({\bA})] = 2;
    \\]
    For $$\lambda > 1$$, we have
    \\[
    \lim_{n \to \infty} \E[\lambda_{\max} ({\bA})] = \lambda + 1/\lambda.
    \\]

-   The correlation of top eigenvector $$\bv_{\max}(\bA)$$ with $$\bu$$: <br/>
    For $$0 \le \lambda \le 1$$, we have
    \\[
    \lim_{n \to \infty} \E[\langle \bv_{\max}(\bA), \bu \rangle^2] = 0;
    \\]
    For $$\lambda > 1$$, we have
    \\[
    \lim_{n \to \infty} \mathbb{E}[\langle {\boldsymbol v}_{\max}({\boldsymbol A}), {\boldsymbol u}\rangle^2] = 1 - 1/\lambda^2.
    \\]

In the rest of this blog post, we will derive this result using the replica method. 
    

## 3. Tricks in statistical physics


### 3.1. The free energy trick


Let $$H: \Omega \to \R$$ and $$f: \Omega \to \R$$ be two functions. We would like to calculate the following quantities
\\[
\begin{aligned}
&\max_{\bsigma \in \Omega} H(\bsigma),\\
&f(\arg \max_{\bsigma} H(\bsigma)).
\end{aligned}
\\]

The free energy trick works as follows. Define a Gibbs measure
\\[
\mu_{\beta, h}(\de \bsigma) = \frac{1}{Z(\beta, h)} \exp\\{ \beta H(\bsigma) + h f(\bsigma)\\} \nu_0(\de \bsigma),
\\]
where $$Z(\beta, h)$$ denotes the normalizing constant
\\[
\begin{align}\tag{Def Z} \label{eqn:def_Z}
Z(\beta, h) = \int_{\Omega} \exp\\{ \beta H(\bsigma) + h f(\bsigma)\\} \nu_0(\de \bsigma),
\end{align}
\\]
and $$\nu_0$$ denotes a reference measure on $$\Omega$$. In physics, $$\beta$$ stands for the inverse temperature, and $$h$$ stands for the strength of external field. Define the free energy function by 
\\[
\Psi(\beta, h) = \log Z(\beta, h). 
\\]
The following lemma shows that, the free energy function $$\Psi$$ can generate many interesting quantities, including $$\max_{\bsigma \in \Omega} H(\bsigma)$$ and $$f(\arg \max_{\bsigma} H(\bsigma))$$. 

<div class="lemma">
Denote the ensemble average operator <script type="math/tex">\langle \cdot \rangle</script> by
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
Equation \eqref{eq:1} follows from basic calculus of exponential family. The intuition for Equation \eqref{eq:2} gives the following: for $$h = 0$$ and large $$\beta$$, the measure $$\mu_{\beta, 0}$$ concentrates at the maxima of $$H$$. 



<!-- ### 3.2. Free energy trick for spiked GOE model -->

In the spiked GOE model, the Hamiltonian associated with the random matrix $$\bA$$ is defined as
\\[
H_{n, \lambda}(\bsigma) = \langle \bsigma, \bA \bsigma \rangle = \lambda \langle \bu, \bsigma\rangle^2 + \langle \bsigma, \bW \bsigma\rangle. 
\\]
Denote the (random) partition function associated with the Hamiltonian $$H_{n, \lambda}$$ by
\\[
Z(\beta, \lambda, n) = \int_{\S^{n-1}} \exp\\{ \beta n H_{n, \lambda}( \bsigma ) \\} \nu_0(\de \bsigma).
\\]
Here $$\nu_0$$ is the uniform probability measure on $$\S^{n-1}$$, and the $$n$$ factor in the exponent serves for proper normalization. 

Define 
\\[
\begin{align}\tag{3}\label{eq:3}
\vphi_n(\lambda) =& \lim_{\beta \to \infty} \frac{1}{\beta n}\E[\log Z(\beta, \lambda, n)],\\\
\vphi(\lambda) =& \lim_{n\to \infty} \vphi_n(\lambda) = \lim_{n \to \infty} \lim_{\beta \to \infty} \frac{1}{\beta n}\E[\log Z(\beta, \lambda, n)]. 
\end{align}
\\]
The function $$\vphi(\lambda)$$ describes the asymptotics of $$\E[\lambda_{\max}(\bA)]$$ and $$\E[\langle \bv_{\max}(\bA), \bu\rangle^2]$$, as shown in the lemma below. 

<div class="lemma">
For <script type="math/tex">\bA = \lambda \bu \bu^\sT + \bW \in \R^{n \times n}</script>, we have
\[
\begin{aligned}
\vphi_n(\lambda) =& \E[\sup_{\bsigma} H_{n, \lambda}(\bsigma)]= \E[\lambda_{\max} (\bA)],\\
\vphi_n'(\lambda) =& \E[\langle \bv_{\max}(\bA), \bu \rangle^2]. 
\end{aligned}
\]
</div>
We leave the proof of this lemma to readers. 


### 3.2. The replica trick

The difficulty of calculating $$\vphi(\lambda)$$ originated from the expectation of $$\log Z$$ as in Eq. \eqref{eq:3}. The replica trick can convert the problem of calculating the expectation of $$\log Z$$ into calculating the moments of $$Z$$. Calculating the moments of $$Z$$ is an easier problem. 

In the language of mathematics, the replica trick gives the following formula
\\[\tag{4} \label{eq:4}
\E[\log Z] = \lim_{k \to 0} \frac{1}{k}\log \E[Z^k ]. 
\\]
Here is a simple proof of this equality:
\\[
\begin{aligned}
\E[\log Z] =& \E[(\log Z^k) / k] = \lim_{k \to 0} \E [ \log(1 + Z^k - 1) / k] = \lim_{k \to 0} \E [ ( Z^k - 1 ) / k ] \\\
=&\lim_{k \to 0} (\E [ Z^k] - 1)/k = \lim_{k \to 0} \frac{1}{k} \log (1 + \E[Z^k] - 1) = \lim_{k \to 0} \frac{1}{k}\log \E[Z^k ]. 
\end{aligned}
\\]


<!-- ### 3.4. The replica trick for spiked GOE model -->

Using the definition of $$\vphi$$ in Eq. \eqref{eq:3} and the replica trick in Eq. \eqref{eq:4}, we can rewrite $$\vphi$$ as
\\[
\begin{aligned}
\vphi(\lambda) \equiv& \lim_{n\to \infty} \lim_{\beta \to \infty} \frac{1}{\beta n}\E[\log Z(\beta, \lambda, n)] = \lim_{n \to \infty}   \lim_{\beta \to \infty} \lim_{k \to 0}\frac{1}{\beta k n}\log \E[Z(\beta, \lambda, n)^k]. \\\
\end{aligned}
\\]
With some heuristic change of limit (which we will not justify), we get
\\[
\vphi(\lambda) \stackrel{\cdot}{=} \lim_{\beta \to \infty}  \lim_{k \to 0} \lim_{n \to \infty} \frac{1}{\beta k n}\log \E[Z(\beta, \lambda, n)^k]. 
\\]
There are three limits in the right hand side of the above equation. Define
\\[
\begin{align}\tag{Def S} \label{eqn:defS}
S(k, \beta, \lambda) = \lim_{n \to \infty} \frac{1}{n}\log \E[Z(\beta, \lambda, n)^k],
\end{align}
\\]
and
\\[
\begin{align}\tag{Def $\Psi$} \label{eqn:defPsi}
\Psi(\beta, \lambda) = \lim_{k \to 0} \frac{1}{k} S(k, \beta, \lambda).
\end{align}
\\]
This gives
\\[
\begin{align}\tag{Def $\vphi$} \label{eqn:defvphi}
\vphi(\lambda) = \lim_{\beta \to \infty} \frac{1}{\beta}\Psi(\beta, \lambda). 
\end{align}
\\]
In the following section, we calculate $$S$$, $$\Psi$$, and $$\vphi$$ functions sequentially. 



## 4. The replica calculations

### 4.1. The large $$n$$ limit: moments calculation

#### 4.1.1. A variational formula for $$S$$ function

Recall the definition of $$S$$ function given in Eq. \eqref{eqn:defS}. We claim the following equality for $$k \in \N_+$$
\begin{align}\tag{5}\label{eqn:S_function}
S(k, \beta, \lambda) \equiv \lim_{n \to \infty} \frac{1}{n}\log\E[Z(\beta, \lambda, n)^k]= \sup_{\bQ \in \R^{(k+1) \times (k+1)}, \diag(\bQ) = 1, \bQ \succeq \bzero} \{U(\bQ) \}, 
\end{align}
where
\begin{align}\tag{6}\label{eqn:U_function}
U(\bQ) = \beta \lambda \sum_{i=1}^k q_{0i}^2 + \beta^2 \sum_{ij = 1}^k q_{ij}^2 + \frac{1}{2} \log(\det(\bQ)), 
\end{align}
and $$\bQ \in \R^{(k + 1) \times (k + 1)}$$ is symmetric, with (the index of $$\bQ$$ starts from $$0$$ and ends at $$k$$)
\\[
\begin{align}\tag{Def $\bQ$} \label{eqn:def_Q}
\bQ = \begin{bmatrix}
1& q_{01} & q_{02} & \ldots & q_{0k}\\\
q_{01} & 1 & q_{12} & \ldots & q_{1k}\\\
q_{02} & q_{12} & 1 & \ldots & q_{2k} \\\
\ldots & \ldots & \ldots & \ldots & \ldots\\\
q_{0k} & q_{1k} & q_{2k}& \ldots & 1\\\
\end{bmatrix}.
\end{align}
\\]
Equation \eqref{eqn:S_function} is exact and rigorous, but the computation is involved. By calculating the moments of $$Z$$ directly (the calculation is given in Section [4.1.2](#sec:derivation)), we get
\begin{align}\tag{7} \label{eqn:moments_result}
\E[Z(\beta, \lambda, n)^k] = \int \exp\\{ n U(\bQ) + o(n)\\} \de \bQ. 
\end{align}
Eq. \eqref{eqn:S_function} follows from Eq. \eqref{eqn:moments_result} and Laplace method. We suggest the readers skipping the derivation of Eq. \eqref{eqn:moments_result} at the first time reading, and continue to read Section [4.2](#sec:rsa). 


#### 4.1.2. The replica calculations 
{: #sec:derivation}

In this section, we derive Eq. \eqref{eqn:moments_result}. At the first time reading, we suggest the readers skipping this section and continue to read Section [4.2](#sec:rsa). 

Recall the definition of $$Z$$ given in Eq. \eqref{eqn:def_Z}. The $$k$$'th moment of $$Z(\beta, \lambda, n)$$ gives 
\\[
\begin{aligned}
\E[Z(\beta, \lambda, n)^k] = & \E\Big[\Big(\int_{\S^{n-1}} \exp\\{ \beta n H_{n, \lambda}( \bsigma ) \\} \nu_0(\de \bsigma)\Big)^k \Big].\\
\end{aligned}
\\]
Here, the expectation $$\E$$ is taken with respect to the randomness of matrix $$\bA$$ in $$H_{n, \lambda}(\bsigma) = \langle \bsigma, \bA \bsigma\rangle$$. The first trick is to rewrite the $$k$$'th power of integration into integration with respect to $$k$$ replicas of $$\bsigma$$. In this way, we can exchange the expectation with integrals. 
\\[
\begin{aligned}
\E[Z(\beta, \lambda, n)^k] =&\E\Big[\Big(\int_{\S^{n-1}} \exp\\{ \beta n H_{n, \lambda}( \bsigma ) \\} \nu_0(\de \bsigma)\Big)^k \Big] = \E\Big[\int_{(\S^{n-1})^k} \exp\Big\\{ \beta n \sum_{i = 1}^k H_{n, \lambda}( \bsigma_i ) \Big\\} \prod_{i \in [k]} \nu_0(\de \bsigma_i) \Big]\\\
=&\int_{(\S^{n-1})^k} \E\Big[ \exp\Big\\{ \beta n \sum_{i = 1}^k  [\lambda \langle \bu, \bsigma_i\rangle^2 + \langle \bsigma_i, \bW \bsigma_i\rangle ] \Big\\}  \Big] \prod_{i \in [k]} \nu_0( \de \bsigma_i).
\end{aligned}
\\]
Let $$\bG = (G_{ij})_{i, j \in [n]} \in \R^{n \times n}$$ with $$G_{ij} \sim_{iid} \cN(0, 1)$$. Noting that the GOE matrix $$\bW$$ shares the same distribution with matrix $$(\bG + \bG^\sT) / \sqrt{2 n}$$, we have
\\[
\begin{aligned}
\E[Z(\beta, \lambda, n)^k] =& \int_{(\S^{n-1})^k} \E_{G_{ij} \sim \cN(0,1)}\Big[ \exp\Big\\{ \beta n \lambda \sum_{i = 1}^k  \langle \bu, \bsigma_i\rangle^2 + \beta n \Big\langle  (\bG + \bG^\sT) / \sqrt{2n},  \sum_{i=1}^k \bsigma_i\bsigma_i^\sT \Big\rangle \Big\\} \Big] \prod_{i \in [k]} \nu_0(\de \bsigma_i) \\\
=& \int_{(\S^{n-1})^k} \exp\Big\\{ \beta n \lambda \sum_{i = 1}^k  \langle \bu, \bsigma_i\rangle^2 \Big\\} \times \E_{G_{ij} \sim \cN(0,1)}\Big[ \exp\Big\\{ \beta \sqrt{2 n} \Big\langle \bG,  \sum_{i=1}^k \bsigma_i\bsigma_i^\sT \Big\rangle \Big\\} \Big] \prod_{i \in [k]} \nu_0(\de \bsigma_i). 
\end{aligned}
\\]
Using the formula of Gaussian moment generating function $$\E_{G \sim \cN(0, 1)}[e^{tG}] = e^{t^2/2}$$, we get
\\[
\begin{aligned}
\E[Z(\beta, \lambda, n)^k] =& \int_{(\S^{n-1})^k} \exp\Big\\{ \beta n \lambda \sum_{i = 1}^k  \langle \bu, \bsigma_i \rangle^2  + n \beta^2 \Big \Vert \sum_{i=1}^k \bsigma_i \bsigma_i^\sT \Big\Vert_F^2 \Big\\}  \prod_{i \in [k]} \nu_0(\de \bsigma_i) \\\
=& \int_{(\S^{n-1})^k} \exp\Big\\{ \beta n \lambda \sum_{i = 1}^k  \langle \bu, \bsigma_i \rangle^2 + n \beta^2 \sum_{i, j=1}^k\Big\langle  \bsigma_i\bsigma_i^\sT, \bsigma_j \bsigma_j^\sT\Big \rangle \Big\\} \prod_{i \in [k]} \nu_0(\de \bsigma_i) \\\
=& \int_{(\S^{n-1})^k} \exp\Big\\{ \beta n \lambda \sum_{i = 1}^k  \langle \bu, \bsigma_i \rangle^2 + n \beta^2 \sum_{i, j=1}^k \langle \bsigma_i, \bsigma_j \rangle^2 \Big\\} \prod_{i \in [k]} \nu_0(\de \bsigma_i). 
\end{aligned}
\\]




In the following, we make use a heuristic argument (using Dirac delta function). Note we have
\\[
1 = \int \prod_{i=1}^k \delta \Big( \langle \bu,\bsigma_i\rangle - n q_{0i} \Big) \prod_{1 \le i < j\le k} \delta\Big( \langle\bsigma_i, \bsigma_j\rangle - n q_{ij} \Big) \de \bQ
\\]
where $$\de \bQ \equiv \prod_{i \in [k]} \de q_{0i} \prod_{1 \le i < j \le k} \de q_{ij}$$. Using this equality, we get
\\[
\begin{aligned}
\E[Z(\beta, \lambda, n)^k] =&  \int \de \bQ \exp\Big\\{\beta n \lambda \sum_{i = 1}^k  q_{0i}^2 + n \beta^2  \sum_{i, j=1}^k q_{ij}^2 \Big\\} \cdot \int_{(\R^n)^k} \prod_{i=1}^k \delta \Big( \langle \bu,\bsigma_i\rangle - n q_{0i} \Big) \prod_{1 \le i < j\le k} \delta\Big( \langle\bsigma_i, \bsigma_j\rangle - n q_{ij} \Big) \prod_{i \in [k]} \nu_0(\de \bsigma_i) \\\
=& \int \de \bQ  \exp\Big\\{\beta n \lambda \sum_{i = 1}^k  q_{0i}^2 + n \beta^2  \sum_{i, j=1}^k q_{ij}^2 + n H_n(\bQ) \Big\\}, \\
\end{aligned}
\\]
where $$\bQ \in \R^{(k + 1) \times (k + 1)}$$ is given by Eq. \eqref{eqn:def_Q} and
\\\[
\begin{aligned}
H_n(\bQ) =& \frac{1}{n} \log \int_{(\R^n)^k}  \prod_{i=1}^k \delta \Big( \langle \bu,\bsigma_i\rangle - n q_{0i} \Big) \prod_{1 \le i < j\le k} \delta\Big( \langle\bsigma_i, \bsigma_j\rangle - n q_{ij} \Big) \prod_{i \in [k]} \nu_0(\de \bsigma_i) \\\
=& \frac{1}{n} \log \int_{(\R^n)^{k+1}} \prod_{0 \le i < j\le k} \delta\Big( \langle\bsigma_i, \bsigma_j\rangle - n q_{ij} \Big) \prod_{0 \le i \le k} \nu_0( \bsigma_i), 
\end{aligned}
\\\]
which serves as the rate function of the random matrix $$\bSigma = (\langle \bsigma_i, \bsigma_j \rangle / n)_{0 \le i, j \le k} \in \R^{(k+1) \times (k+1)}$$ when $$(\bsigma_i)_{0 \le i \le k} \sim_{iid} \Unif(\S^{n-1})$$. This rate function can be found in standard textbooks
\\[
\lim_{n \to \infty} H_n(\bQ) = \frac{1}{2}\log \det(\bQ). 
\\]
Here we give a heuristic argument to calculate the entropy term $$H_n(\bQ)$$, (here the argument is handwaving, but can be made rigorous)
<p>
\[
\begin{aligned}
H_n(\bQ) \stackrel{\cdot}{=}& \frac{1}{n} \P_{(\bsigma_i)_{i \in [k]} \sim \Unif(\S^{n-1})} \Big(  \big(\langle \bsigma_i, \bsigma_j \rangle / n \big)_{i, j \in [k]}\approx \bQ  \Big)\\\
\stackrel{\cdot}{=}& \inf_{\lambda_{ij}} \frac{1}{n} \log  \int_{(\R^n)^{k+1}} \prod_{0 \le i \le j\le k} \exp\Big\{ - \lambda_{ij}\langle\bsigma_i, \bsigma_j\rangle / 2 + n q_{ij} \lambda_{ij} / 2 \Big\} \prod_{0 \le i \le k} \de \bsigma_i + \const\\\
=& \inf_{\lambda_{ij}} \log  \int_{\R^{k+1}} \prod_{0 \le i \le j\le k} \exp\Big\{ -  \lambda_{ij} \sigma_i\sigma_j / 2 + q_{ij} \lambda_{ij}/2 \Big\} \prod_{0 \le i \le k} \de \sigma_i + \const\\\
=& \inf_{\bLambda = (\lambda_{ij})_{0 \le i \le j \le k}} \Big[ \langle \bQ, \bLambda\rangle / 2 - \log(\det(\bLambda))/2  \Big] + \const\\\
\stackrel{\cdot}{=}& \frac{1}{2} \log \det(\bQ),
\end{aligned}
\]
</p>
where we use the fact that, for $$\bX_i \sim_{iid} \bX$$, we have
\\[
\lim_{n \to \infty} \frac{1}{n} \log \P_{\bX_i}\Big( \frac{1}{n} \sum_{i=1}^n \bX_i \approx \bx \Big) = \inf_{\blambda} \log \E_\bX \Big[\exp\\{ \langle\blambda, \bX - \bx\rangle \\} \Big]. 
\\]
This gives Eq. \eqref{eqn:moments_result}. 




### 4.2. The small $$k$$ limit: replica symmetric ansatz 
{: #sec:rsa}

#### 4.2.1. Replica symmetric ansatz

Next we calculate (recall Eq. \eqref{eqn:defPsi})
\\[
\Psi(\beta, \lambda) = \lim_{k \to 0} S(\beta, \lambda, k)/k. 
\\]
We gave the expression of $$S$$ when $$k$$ is integral (c.f. Eq. \eqref{eqn:S_function}), where $$k$$ serves as the dimension of variable $$\bQ$$ in the variational formula. However, to take the limit $$k \to 0$$, we need the expression for $$S$$ when $$k$$ takes real numbers. This is the difficulty to calculate the small $$k$$ limit using Eq. \eqref{eqn:S_function}. 


The physicists' trick is to plug in an ansatz for the matrix $$\bQ$$ to simplify the expression of $$S$$. Note that $$U(\bQ)$$ defined in Eq. \eqref{eqn:U_function} has some symmetry: if we permute the $$1$$'st to $$k$$'th rows of $$\bQ$$, and perform the same permutation to its columns, the function $$U$$ stays the same. This motivates us to consider the "replica symmetric ansatz": 
\\[
\bQ = \begin{bmatrix}
1& \mu & \mu & \ldots & \mu\\\
\mu & 1 & q& \ldots & q\\\
\mu & q & 1 & \ldots & q \\\
\ldots & \ldots &\ldots & \ldots & \ldots\\\
\mu & q & q & \ldots & 1\\\
\end{bmatrix},
\\]
or equivalently 
\\[
q_{0i} = \mu,  ~~ \forall i \in [k], ~~~~~~~~~~~~~~ q_{ij} = q, ~~ \forall 1 \le i < j \le k.
\\]
Plugging this ansats into Eq. \eqref{eqn:U_function}, and simplifying the log determinant term using the matrix determinant lemma 
\\[
\begin{aligned}
&\log(\det(\bQ)) = \log \det( (1-q) \id_k + (q - \mu^2) \ones \ones^\sT) \\\
=&  \log \det( \id_k + [(q - \mu^2)/(1-q) ] \ones \ones^\sT) + k \log(1-q)\\\
=& \log(1 + k [(q - \mu^2)/(1-q) ]) + k \log(1-q)\\\
=& \log\Big( 1 - \frac{\mu^2 k}{1 + (k - 1)q} \Big) + \log \Big(1 + \frac{kq}{1-q} \Big) + k \log(1 - q),
\end{aligned}
\\]
we get
\\[
U(\bQ) = \beta \lambda k \mu^2 + \beta^2 k + \beta^2 k (k-1) q^2 + \frac{1}{2} \Big[ \log\Big( 1 - \frac{\mu^2 k}{1 + (k - 1)q} \Big) + \log \Big(1 + \frac{kq}{1-q} \Big) + k \log(1 - q) \Big] \equiv U(\mu, q, k). 
\\]

Assuming that the maximum of $$U$$ is taken at a "replica symmetric" form of $$\bQ$$, we get
\\[
S(\beta, \lambda, k) \stackrel{\cdot}{=} \sup_{\mu, q} U(\mu, q, k).
\\]



#### 4.2.2. The small $$k$$ limit
Using the replica symmetric ansatz, we gave an expression for $$S$$ for any $$k > 0$$. Next we would like to calculate the $$k \to 0$$ limit (recall the definition of $$\Psi$$ given in Eq. \eqref{eqn:defPsi})
\\[
\Psi(\beta, \lambda) = \lim_{k \to 0} S(\beta, \lambda, k) / k = \lim_{k \to 0} \sup_{\mu, q} \frac{1}{k} U(\mu, q, k). 
\\]
Define 
\\[
u(\mu, q) \equiv \lim_{k \to 0} \frac{1}{k} U(\mu, q, k) = \beta \lambda \mu^2  + \beta^2 (1 - q^2) - \frac{\mu^2}{2(1 - q)} + \frac{q}{2 (1 - q)} + \frac{1}{2}\log(1 - q). 
\\]
We need to exchange the operations $$\lim_{k \to 0}$$ and $$\sup_{\mu, q}$$ in some way. Another heuristic physicists' convention comes here: the $$\lim_{k \to 0} \sup_{\mu, q}$$ operation becomes the $$\ext_{\mu, q} \lim_{k \to 0}$$ operation, where $$\ext_{\bx} f(\bx)$$ is a set defined as
\\[
\ext_{\bx} f(\bx) = \Big\\{ f(\bx_\star): \nabla f(\bx_\star) = \bzero\Big\\}. 
\\]
Using this convention, we get
\\[
\Psi(\beta, \lambda) =  \lim_{k \to 0} S(\beta, \lambda, k) / k \in \ext_{\mu, q}u(\mu, q),
\\]

Next, we calculate all the stationary points of $$u$$, and collect the value of $$u$$ evaluated at these stationary points. 

By basic calculus, we have 
\\[
\begin{aligned}
\partial_\mu u(\mu, q) =& 2 \Big( \beta \lambda - \frac{1}{2 (1-q)} \Big) \mu, \\\
\partial_q u(\mu, q) =& -2 \beta^2 q -  \frac{\mu^2}{2(1 - q)^2} + \frac{1}{2 (1 - q)^2}.
\end{aligned}
\\]
When $$\beta > 1$$, by basic algebra, we find two solutions of $$\nabla u = \bzero$$. 

- One branch of solution of $$\nabla u = \bzero$$ gives 
	\\[
	\begin{aligned}
	\mu_1 =& 0, \\\
	q_1 =& 1 - \frac{1}{2 \beta}.
	\end{aligned}
	\\]
	At this solution, we have 
	\\[
	\Psi_1(\beta, \lambda) = u(\mu_1, q_1) = 2 \beta - \frac{3}{4} - \frac{1}{2} \log(2 \beta).
	\\]
- Another branch of solution of $$\nabla u = \bzero$$ gives 
	\\[
	\begin{aligned}
	\mu_2 =& \Big( \Big(1 - \frac{1}{\lambda^2}\Big) \Big(1 - \frac{1}{2 \beta \lambda}\Big)  \Big)^{1/2},\\\
	q_2 =& 1 - \frac{1}{2 \lambda \beta},
	\end{aligned}
	\\]
	At this solution, we have
	\\[
	\Psi_2(\beta, \lambda) = u(\mu_2, q_2) = \beta\Big( \lambda + \frac{1}{\lambda}\Big) - \Big( \frac{1}{4 \lambda^2} + \frac{1}{2} \Big) - \frac{1}{2}\log(2 \lambda \beta). 
	\\]

Hence we get
\\[
\Psi \in \ext_{\mu, q}u(\mu, q) = \\{\Psi_1, \Psi_2 \\}.
\\]



### 4.3. The large $$\beta$$ limit

Recall the definition of $$\vphi$$ in Eq. \eqref{eqn:defvphi}. For the first branch $$\Psi_1$$, the corresponding $$\vphi$$ function gives
\\[
\vphi_1(\lambda) = \lim_{\beta \to \infty} \frac{1}{\beta} \Psi_1(\beta, \lambda) = 2. 
\\]
For the second branch $$\Psi_2$$, the corresponding $$\vphi$$ function gives 
\\[
\vphi_2(\lambda) = \lim_{\beta \to \infty} \frac{1}{\beta} \Psi_2(\beta, \lambda) = \lambda + \frac{1}{\lambda}. 
\\]

### 4.4. Select a branch

There are two branches of solutions. We need to select the branch that makes sense: we hope that $$\vphi(\lambda)$$ gives an expression for $$\lim_{n\to \infty}\E[\lambda_\max(\bA)] = \lim_{n\to \infty}\E[\lambda_\max(\lambda \bu \bu^\sT + \bW)]$$, hence it should be non-decreasing in $$\lambda$$, and it should satisfy $$\lim_{\lambda \to \infty} \vphi(\lambda) = \infty$$. 
- When $$\lambda\le 1$$, the branch $$\vphi_2(\lambda) = \lambda + 1/\lambda$$ is decreasing. Hence, we should select the branch $$\vphi_1(\lambda) = 2$$. 
- When $$\lambda \ge 1$$, the branch $$\vphi_1(\lambda) = 2$$ stays constant. Hence, we should select the branch $$\vphi_2(\lambda) = \lambda + 1/\lambda$$. 


The discussion above gives the following prediction 
\\[
\lim_{n \to \infty} \E[\lambda_{\max}(\bA)]= \vphi(\lambda) = \begin{cases}
2, ~~~& \lambda \le 1, \\\
\lambda + \frac{1}{\lambda}, ~~~& \lambda > 1.
\end{cases}
\\]
This turns out to be the correct answer! 



## 5. Summary

In summary, we calculated the largest eigenvalue of the spiked GOE matrix $$\bA = \lambda \bu \bu^\sT + \bW$$. Using the free energy trick and the replica trick, the problem becomes calculating three limits sequentially: the large $$n$$ limit, the small $$k$$ limit, and the large $$\beta$$ limit. The large $$n$$ limit calculation is rigorous, but is also complicated. The small $$k$$ limit calculation is non-rigorous: we pluged in the replica symmetric ansatz, and changed the $$\lim_{k \to 0} \inf_{\mu, q}$$ operation to the $$\ext_{\mu, q} \lim_{k \to 0}$$ operation. The large $$\beta$$ limit calculation is straightforward. Finally, we get two branches of solutions. We used some simple properties of the largest eigenvalue of $$\bA$$ to select the branch when $$\lambda$$ is in different regime. The prediction formula turns out to be the correct answer. 



## References

[1] Mezard, Marc, and Andrea Montanari. Information, physics, and computation. Oxford University Press, 2009. 
{: #ref1}

[2] Baik, Jinho, Gérard Ben Arous, and Sandrine Péché. "Phase transition of the largest eigenvalue for nonnull complex sample covariance matrices." The Annals of Probability 33.5 (2005): 1643-1697. 
{: #ref2}

[3] Péché, Sandrine. "The largest eigenvalue of small rank perturbations of Hermitian random matrices." Probability Theory and Related Fields 134.1 (2006): 127-173.
{: #ref3}
