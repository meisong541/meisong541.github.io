---
layout: post
title:  "Replica method and random matrices (II)"
date:   2019-09-12 00:09:06 -0700
categories: jekyll update
---





$$\def\cP{\mathcal P}$$
$$\def\de{\text{d}}$$
$$\def\N{\mathbb N}$$
$$\def\P{\mathbb P}$$
$$\def\bi{\boldsymbol i}$$
$$\def\bG{\boldsymbol G}$$
$$\def\bsigma{\boldsymbol \sigma}$$
$$\def\bv{\boldsymbol v}$$
$$\def\bu{\boldsymbol u}$$
$$\def\bt{\boldsymbol t}$$
$$\def\bw{\boldsymbol w}$$
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

In my previous post, we saw how to use the replica method to calculate the spectral norm of the spiked GOE matrix. Besides the spectral norm, another interesting quantity of random matrices is the spectral density: the limiting empirical distribution of eigenvalues. Similarly using the replica method, we show in this post how to calculate the spectral density of random matrices. In particular, we will illustrate the method through the example of GOE matrix. 

The mathematics of this post builds upon Stieltjes transforms, which we will not dive into for the sake of brevity. For a detailed introduction to Stieltjes transforms and its applications to the GOE matrix, I recommend reading <a href = "https://terrytao.wordpress.com/2010/02/02/254a-notes-4-the-semi-circular-law/#more-3426">Terrence Tao's blog</a>. It is worth noting that the shape of the spectral density of the GOE matrix is a semicircle, and hence this famous result is called the <a href="https://en.wikipedia.org/wiki/Wigner_semicircle_distribution">semicircle law</a>. 



## 2. Spectral density and Stieltjes transforms

Let $$\bA_n \in \R^{n \times n}$$ be a symmetric matrix. Let $$z \in \C_+ = \{ \xi \in \C: \Im \xi > 0 \}$$ ($$\Im \xi$$ is the imaginary part of $$\xi$$). We denote the resolvent of $$\bA_n$$ by 
<p>
\[
R_{\bA_n}(z) = (\bA_n - z \id_n)^{-1}.
\]
</p>
Denoting the $$i$$'th largest eigenvalue of $$\bA_n$$ by $$\lambda_i(\bA_n)$$, and the empirical eigenvalue distribution (spectral density) of $$\bA_n$$ by 
\\[
\mu_{\bA_n} = \frac{1}{n} \sum_{i=1}^n \delta_{\lambda_i(\bA_n)}. 
\\]


For a probability measure $$\mu \in \cP(\R)$$ on the real line, we denote its Stieltjes transform by
\\[
s_{\mu}(z) = \int \frac{1}{\lambda - z} \mu(\de \lambda). 
\\]
The next two lemmas characterize the properties of Stieltjes transforms. 

<div class="lemma"><ri>[Inverse Stieltjes transforms]</ri>
The scaled imaginary part of the Stieltjes transform <script type="math/tex">(1/ \pi) \Im s_{\mu}(\cdot + \bi b)</script> for <script type="math/tex">b>0</script> is a probability measure on the real line <script type="math/tex">\R</script>. The sequence of probability measures <script type="math/tex">(1/ \pi) \Im s_{\mu}(\cdot + \bi b)</script> indexed by <script type="math/tex">b</script> converges weakly to <script type="math/tex">\mu</script> as <script type="math/tex">b \to 0+</script>. 
</div>

<div class="lemma"><ri>[Stieltjes continuity]</ri>
Let <script type="math/tex">\mu_n</script> be a sequence of random probability measures on the real line, and let <script type="math/tex">\mu</script> be a deterministic probability measure. Then <script type="math/tex">\mu_n</script> converges in expectation to <script type="math/tex">\mu</script> in the vague topology if and only if <script type="math/tex">\E s_{\mu_n}(z)</script> converges to <script type="math/tex">s_\mu(z)</script> for every <script type="math/tex">z</script> in the upper half-plane.
</div>



Building upon the previous two lemmas, the next theorem gives the limiting Stieltjes transform of the GOE matrix, and hence recover the semicircle law. 


<div class="theorem"><ri>[Stieltjes transform of the GOE matrix]</ri>
Let <script type="math/tex">\bW_n \sim \GOE(n)</script>. That is, <script type="math/tex">\bW_n \in \R^{n \times n}</script> is a symmetric matrix with <script type="math/tex">(W_{n, ii})_{1 \le i \le n} \sim_{i.i.d.} \cN(0, 2/n)</script> and <script type="math/tex">(W_{n, ij})_{1 \le i < j \le n} \sim_{i.i.d.} \cN(0, 1/n)</script>. Denote by <script type="math/tex">s_{\bW_n} = s_{\mu_{\bW_n}}</script> the Stieltjes transform of empirical eigenvalue distribution of <script type = "math/tex">\bW_n</script>. Then for <script type="math/tex">z \in \C_+</script>, we have 
\[
\lim_{n \to \infty}\E[s_{\bW_n}(z)] = \frac{-z + \sqrt{z^2 - 4}}{2}. 
\]
</div>



Combining the above theorem with the Stieltjes continuity lemma, for $$a < b$$, we have
\\[
\lim_{n \to \infty} \E[\\# \\{ \lambda_i(\bW_n): \lambda_i(\bW_n) \in [a, b] \\} / n] =  \sigma_{\rm sc}([a, b]),
\\]
where
\\[
\sigma_{\rm sc}(\de x) = \frac{1}{2\pi} \sqrt{4 - x^2} \cdot \ones\\{ x \in [-2, 2] \\} \de x. 
\\]

We have so far established the background of this example. It remains to calculate formally the limiting Stieltjes transform $$\lim_{n \to \infty}\E[s_{\bW_n}(z)]$$ using the replica method, which will be the focus of the following sections. 


## 3. The determinant trick 

### 3.1. The connection of determinant and Stieltjes transform

Let $$\log$$ be the complex log function defined on $$\C$$, with branch cut along the negative real axis. Define function $$D_{\bW_n}$$ as 
\\[
D_{\bW_n}(\xi) = \frac{1}{n} \sum_{i=1}^n \log(\lambda_i(\bW_n) - \xi).
\\]
Then the derivative of $$D_{\bW_n}$$ gives the negative Stieltjes transform
\\[
\frac{\de}{\de \xi} D_{\bW_n}(\xi) = - \frac{1}{n} \sum_{i=1}^n \frac{1}{\lambda_i(\bW_n) - \xi} = - s_{\bW_n}(\xi).
\\]
The function $$D_{\bW_n}$$ is almost the normalized log-determinant of $$\bW_n$$, up to a phase shift
\\[
D_{\bW_n}(\xi) = \frac{1}{n}\log \det(\bW_n - \xi \id_n) + \frac{2 \pi \bi k(\bW_n, \xi)}{n},
\\]
where $$k(\bW_n, \xi)$$ is an integer. Moreover, $$k(\bW_n, \xi)$$ will remain the same under an infinitesimal change of $$\xi$$. Hence we have 
\begin{align}
&~\lim_{n \to \infty} \E[s_{\bW_n}(\xi)] = - \lim_{n \to \infty} \frac{\de}{\de \xi} \E[D_{\bW_n}(\xi)] =  - \lim_{n \to \infty} \frac{\de}{\de \xi} \frac{1}{n} \E[\log \det(\bW_n - \xi \id_n)]\nonumber \\\
 \stackrel{\cdot}{=}&~ - \frac{\de}{\de \xi} \lim_{n \to \infty}  \frac{1}{n} \E[\log \det(\bW_n - \xi \id_n)]\tag{1}\label{eq:1}
\end{align}
In the last step, we heuristically exchanged the limit operator and the differential operator. 

The point of this trick is that, determinants are easier to work with than Stieltjes transforms, since the power of determinant can be expressed in terms of integration of exponentials
\\[
\det(\bA)^{-k/2} = \int_{\R^n} \frac{1}{(2 \pi)^{nk/2}} \exp\Big\\{ - \frac{1}{2} \sum_{j = 1}^k \langle \bx_j, \bA \bx_j \rangle  \Big\\} \prod_{j \in [k]} \de \bx_j. 
\\]
This identity holds whenever $$\bA \in \R^{n \times n}$$ is positive semi-definite. Though, we will formally use this identity for a complex matrix $$\bA$$ which may not always hold rigorously. 


### 3.2. The replica approach

In order to calculate $$\E[s_{\bW_n}(\xi)]$$, by Eq. \eqref{eq:1}, we need to calculate $$\lim_{n\to \infty}(1/n)\E[\log \det(\bW_n - \xi \id_n)]$$, and then differentiate with respect to $$\xi$$. However, there is no straightforward way to calculate the expectation of $$\log$$. One possible way is to use the replica formula introduced in my last post
\\[
\E[\log Z] = \lim_{k \to 0} \frac{1}{k} \log \E[Z^k].
\\]
The replica formula reduces the problem to calculating the moments $$\E[\det(\bW_n - \xi \id_n)^k]$$. 

### 3.3. An easier approach

Instead of calculating $$\E[\det(\bW_n - \xi \id_n)^k]$$ for a sequence of $$k$$, we just calculate $$\E[\det(\bW_n - \xi \id_n)^{-1/2}]$$. Note
\\[
\frac{1}{n} \log \det(\bW_n  - \xi \id_n ) = -\frac{2}{n} \log [\det(\bW_n  - \xi \id_n)^{-1/2}].
\\]
We expect that $$n^{-1} \log \det(\bW_n  - \xi \id_n )$$ concentrates tightly around its mean, so that 
\begin{align}
&~\lim_{n \to \infty} \frac{1}{n} \E[ \log \det(\bW_n - \xi \id_n)] \stackrel{\cdot}{=} \lim_{n \to \infty} \frac{1}{n}\log \det(\bW_n - \xi \id_n)] \nonumber\\\
 \stackrel{\cdot}{=}&~ \lim_{n \to \infty} - \frac{2}{n}\log \det(\bW_n - \xi \id_n)^{-1/2}] \stackrel{\cdot}{=} \lim_{n \to \infty} -\frac{2}{n} \log \E[\det(\bW_n  - \xi \id_n)^{-1/2}]. \tag{2}\label{eq:2}
\end{align}
This is completely heuristic, but it finally gives the correct answer. 


## 4. The replica calculations

In this section, we use the replica method to calculate the last expression of Eq. \eqref{eq:2}, and then derive the expression of the limiting Stieltjes transform. We will mostly use formal calculations that are not rigorous. 
 
### 4.1. Step 1: Get rid of the expectation operator
 
The first step of the calculation is to get rid of the expectation operator. Suppose we need to calculate $$\E[f(\bw)]$$, where $$\bw \sim \cN(0, \id_d)$$. First we need to find a way to express $$f$$ as the following integral form
\\[
f(\bw) = \int_{\R^n} \exp\\{ \langle \bw, \bt(\bx)\rangle \\} h(\bx) \de \bx,
\\]
such that using the formula of the moment generating function of Gaussian random variables $$\E_{G \sim \cN(0, 1)}[e^{tG}] = e^{t^2/2}$$, we get
\\[
\E[f(\bw)] = \int_{\R^n} \exp\\{ \Vert \bt(\bx) \Vert_2^2/2 \\} h(\bx) \de \bx. 
\\]

For the current example, note a formal identity gives (it is formal because $$\bW_n - \xi \id_n$$ is a complex matrix)
<p>
\[
\det(\bW_n - \xi \id_n)^{-1/2} = \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \exp\Big\{ - \bx^\sT (\bW_n - \xi \id_n) \bx / 2 \Big\} \de \bx.
\]
</p>
Let $$\bG = (G_{ij})_{ij \in [n]} \in \R^{n \times n}$$ with $$G_{ij} \sim_{iid} \cN(0, 1)$$. Then the distributions of $$\bW_n$$ and $$(\bG + \bG^\sT)/\sqrt{2n}$$ are identital. As a consequence, we get
\\[
\begin{aligned}
\E[\det(\bW_n - \xi \id_n)^{-1/2}] =&~ \E\Big[ \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \exp\Big\\{ - \Big\langle (\bG + \bG^\sT) / \sqrt {2n} - \xi \id_n, \bx \bx^\sT \Big\rangle / 2 \Big\\} \de \bx \Big]\\\
=&~ \E\Big[ \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \exp\Big\\{ - \Big\langle \bG , \bx \bx^\sT \Big\rangle / \sqrt {2n} + \xi \Vert \bx \Vert_2^2 /2 \Big\\} \de \bx \Big]\\\
=&~ \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \exp\Big\\{  \xi \Vert \bx \Vert_2^2 /2 \Big\\} \E\Big[\exp\Big\\{ - \Big\langle \bG , \bx \bx^\sT \Big\rangle / \sqrt {2n} \Big\\}\Big]  \de \bx. 
\end{aligned}
\\]
Using the formula of the moment generating function of Gaussian random variables, we get
\begin{align}
\E[\det(\bW_n - \xi \id_n)^{-1/2}] =&~ \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \exp\Big\\{  \xi \Vert \bx \Vert_2^2 /2 + \Vert \bx \bx^\sT \Vert_F^2 / (4n) \Big\\}  \de \bx \nonumber \\\
=&~ \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \exp\Big\\{  \xi \Vert \bx \Vert_2^2 /2 + \Vert \bx \Vert_2^4 / (4n) \Big\\}  \de \bx.  \tag{3}\label{eq:3} \\\
\end{align}
Now, we have got rid of the expectation operator, and expressed the quantity of interest in terms of an integral. 
 
### 4.2. Step 2: Calculate the integral

Eq. \eqref{eq:3} is an intractable high dimensional integration. In statistical physics, there is a systematic approach to heuristically deal with such an integration. 

#### 4.2.1. A systematic approach to deal with the integration

Suppose $$I_n$$ is an integration of form 
\\[
I_n = \int_{\R^n} \exp\\{ n \cdot f(h_1(\bx), \ldots, h_k(\bx)) \\} \de \bx, 
\\]
where $$h_j$$'s are in the form (with a slight abuse of notation)
\\[
h_j(x_1, \ldots, x_n) = \frac{1}{n}\sum_{i = 1}^n h_j(x_i). 
\\]
First by the identity that 
\begin{align}
\int_\R \delta(n h_j(\bx) - n s_j) \de s_j = 1,
\end{align}
we have 
\\[
I_n =  \int_{\R^k} \prod_{j \in [k]} \de s_j \cdot \exp\\{ n \cdot f(s_1, \ldots, s_k) \\} \int_{\R^n} \prod_{j \in [k]} \delta( n h_j(\bx) - n s_j) \de \bx. 
\\]
Then, by the delta identity formula 
\\[
\delta(x) = \int_{\bi \R} \exp\\{ \lambda x \\} \de [\lambda/ (2 \pi)], 
\\]
we get 
\\[
\begin{aligned}
I_n =&~ \int_{\R^k} \prod_{j \in [k]} \de s_j \cdot \exp\\{ n \cdot f(s_1, \ldots, s_k) \\}  \int_{(\bi \R)^k} \prod_{j \in [k]}  \de [\lambda_j / (2\pi)]  \int_{\R^n} \exp\Big\\{  \sum_{j = 1}^k \sum_{i = 1}^n \Big[ \lambda_j h_j(x_i) - \lambda_j s_j \Big] \Big\\} \de \bx  \\\
=&~ \int_{\R^k} \prod_{j \in [k]} \de s_j \int_{(\bi \R)^k} \prod_{j = 1}^k \de [\lambda_j/(2\pi)]  \exp \Big\\{ n \cdot \Big[ f(s_1, \ldots, s_k) - \sum_{j=1}^k \lambda_j s_j \Big] \Big\\} \cdot \Big( \int \exp\Big\\{ \sum_{j = 1}^n \lambda_j h_j(x) \Big\\} \de x \Big)^n \\\
=&~ \int_{\R^k} \prod_{j \in [k]} \de s_j \int_{(\bi \R)^k} \prod_{j = 1}^k \de [\lambda_j/(2\pi)]  \exp \Big\\{ n \cdot \Big[ f(s_1, \ldots, s_k) - \sum_{j=1}^k \lambda_j s_j  + \log J(\lambda_1, \ldots, \lambda_k) \Big] \Big\\},\\\
\end{aligned}
\\]
where
\\[
J(\lambda_1, \ldots, \lambda_k) = \int \exp\Big\\{ \sum_{j = 1}^n \lambda_j h_j(x) \Big\\} \de x. 
\\]
Using the method of steepest descent (my understanding is that, it is like the Laplace method but it deals with integration of complex functions), we have
\\[
\lim_{n \to \infty} n^{-1} \log I_n \stackrel{\cdot}{=} \ext_{s_j, \lambda_j} \Big[ f(s_1, \ldots, s_k) - \sum_{j=1}^k \lambda_j s_j  + \log J(\lambda_1, \ldots, \lambda_k) \Big],
\\]
where $$\ext$$ denotes the extremum operator. 


#### 4.2.2. Application to this example

<!-- 
Note we have 
\begin{align}\tag{4}\label{eq:4}
\int_\R \delta(\Vert \bx \Vert_2^2 - n s) \de s = 1. 
\end{align}
 -->
In this example, we need to deal with an integration of form Eq. \eqref{eq:3}. We take $$f(s) = \xi s + s^2/4$$ and $$h(\bx) = \sum_{i = 1}^n x_i^2 / n$$. Using the method introduced above, we have
\\[
\begin{aligned}
\E[\det(\bW_n - \xi \id_n)^{-1/2}] =&~ \int_{\R} \de s  \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \exp\Big\\{  \xi \Vert \bx \Vert_2^2 /2 + \Vert \bx \Vert_2^4 / (4n) \Big\\} \cdot \delta(\Vert \bx \Vert_2^2 - n s)  \de \bx \\\
=&~ \int_{\R} \de s  \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \exp \\{  n \xi s /2 + n s^2 / 4 \\} \cdot \delta(\Vert \bx \Vert_2^2 - n s)  \de \bx \\\
=&~ \int_{\R} \de s \int_{\bi R} [\de \lambda / (2 \pi)] \exp \Big\\{  n \Big[\xi s /2 +  s^2 / 4 - \lambda s \Big] \Big\\} \cdot \Big( \int_R \frac{1}{(2 \pi)^{1/2}} \exp\\{ \lambda x^2 \\}  \de x \Big)^n \\\
=&~ \int_{\R} \de s \int_{\bi R} [\de \lambda / (2 \pi)] \exp \Big\\{  n \Big[\xi s /2 +  s^2 / 4 - \lambda s + J(\lambda) \Big] \Big\\} \\\
\end{aligned}
\\]
where 
\\[
J(\lambda) = \log \int_R \frac{1}{(2 \pi)^{1/2}} \exp\\{ \lambda x^2\\}  \de x = - \frac{1}{2} \log ( - 2 \lambda). 
\\]
Therefore, we have
\\[
\lim_{n \to \infty} \frac{1}{n} \log \E[\det(\bW_n - \xi \id_n)^{-1/2}] = \ext_{s, \lambda} \Big[  \frac{ \xi s}{2} + \frac{s^2}{4} - \lambda s - \frac{1}{2} \log( - 2\lambda) \Big]. 
\\]
Define 
\\[
P(s, \lambda) = \frac{ \xi s}{2} + \frac{s^2}{4} - \lambda s - \frac{1}{2} \log( - 2 \lambda). 
\\]
Letting $$\partial_\lambda P(s, \lambda) = 0$$, we have $$\lambda = - 1 / (2s)$$, which gives
\\[
P(s) \equiv \ext_{\lambda} P(s, \lambda) = \frac{ \xi s}{2} + \frac{s^2}{4} + \frac{1}{2} \log(s) + \frac{1}{2}. 
\\]
Then we differentiate $$P(s)$$, its extremum $$s_\star$$ satisfies
\begin{align}\tag{5}\label{eq:5}
\xi + s_\star + \frac{1}{s_\star} = 0. 
\end{align}
As a result, we get
\\[
\lim_{n \to \infty} \frac{1}{n} \log \E[\det(\bW_n - \xi \id_n)^{-1/2}] \stackrel{\cdot}{=} \frac{ \xi s_\star}{2} + \frac{s_\star^2}{4} + \frac{1}{2} \log(s_\star) + \frac{1}{2},
\\]
and by Eq. \eqref{eq:2} we get
\\[
\lim_{n \to \infty} \frac{1}{n}\E[\log \det(\bW_n - \xi \id_n)] \stackrel{\cdot}{=} - \xi s_\star - \frac{s_\star^2}{2} - \log(s_\star) - 1.
\\]
where $$s_\star$$ gives the solution of Eq. \eqref{eq:5}. Finally, by Eq. \eqref{eq:1}, we have 
\\[
\begin{aligned}
\lim_{n \to \infty} \E[s_{\bW_n}(\xi)] \stackrel{\cdot}{=}&~ - \frac{\de}{\de \xi} \lim_{n \to \infty} \frac{1}{n} \E[ \log \det(\bW_n - \xi \id_n)] \\\
=&~ - \frac{\de }{\de \xi} \Big[- \xi s_\star(\xi) - \frac{s_\star^2(\xi)}{2} - \log s_\star(\xi) - 1\Big] = - \frac{\partial }{\partial \xi} \Big[- \xi s - \frac{s^2}{2} - \log s - 1\Big] \Big\vert_{s = s_\star(\xi)} \\\
=&~ s_\star(\xi). 
\end{aligned}
\\]
Let $$\xi \in \C_+$$, and finding the solution $$s_\star$$ of Eq. \eqref{eq:5} with positive imaginary part, we get 
\\[
s_\star(\xi) = \frac{- z + \sqrt{z^2 - 4}}{2}. 
\\]
This is exactly the limiting Stieltjes transform of the GOE matrix. Of course, to rigorously prove this formula, we need to adopt other techniques. 

## 5. Summary 

We showed a systematic approach for formally calculating the spectral density of random matrices. Here is an exercise for the readers: calculate the limiting Stieltjes transform of the empirical covariance matrix $$\hat \Sigma_n$$ of isotropic Gaussian data, where $$\hat \Sigma_n = n^{-1} \sum_{i = 1}^n \bx_i \bx_i^\sT$$ with $$\bx_i \sim_{iid} \cN(0, \id_d)$$, in the asymptotic regime $$n / d \to \gamma$$ as $$d \to \infty$$. The spectral density of the empirical covariance matrix, which can be derived from the exercise above, gives the famous <a href="https://en.wikipedia.org/wiki/Marchenko%E2%80%93Pastur_distribution">Marchenko-Pastur law</a>. 

In my next post, I will introduce some applications of the replica method in machine learning. 

