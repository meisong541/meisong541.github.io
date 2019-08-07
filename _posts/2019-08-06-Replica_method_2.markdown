---
layout: post
title:  "Replica method and random matrices (II)"
date:   2020-08-10 00:09:06 -0700
categories: jekyll update
---


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
$$\def\C{\mathbb C}$$
$$\def\tr{\text{tr}}$$


## 1. Introduction


## 2. Spectral density and Stieltjes transforms

Let $$\bA_n \in \R^{n \times n}$$ be a symmetric matrix. Let $$\xi \in \C_+ = \{ z \in \C: \Im z > 0 \}$$ ($$\Im z$$ is the imaginary part of $$z$$), we denote the resolvent of $$\bA$$ by 
<p>
\[
R_{\bA_n}(\xi) = (\bA_n - \xi \id_n)^{-1}. 
\]
</p>
We denote the Stieltjes transform of $$\bA_n$$ by
<p>
\[
s_{\bA_n}(\xi) = \frac{1}{n}\tr(R_{\bA_n}(\xi)) = \frac{1}{n}\tr\Big((\bA_n - \xi \id_n)^{-1}\Big) = \frac{1}{n} \sum_{i=1}^n \frac{1}{\lambda_i(\bA_n) - \xi}. 
\]
</p>

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