---
layout: post
title:  "Replica method and random matrices (II)"
date:   2020-10-01 00:09:06 -0700
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

In my last blog, the replica method is used to calculate the spectral norm of the spiked GOE matrix. In this blog, I will show how to use the replica method to calculate the spectral density of a random matrix. I will use the GOE matrix as an example. 

Knowledge of Stieltjes transforms is a preliminary requirement for readers. For a basic introduction to Stieltjes transforms and its applications to the GOE matrix, I recommend <a href = "https://terrytao.wordpress.com/2010/02/02/254a-notes-4-the-semi-circular-law/#more-3426">Terrence Tao's blog</a>. 



## 2. Spectral density and Stieltjes transforms

Let $$\bA_n \in \R^{n \times n}$$ be a symmetric matrix. Let $$z \in \C_+ = \{ \xi \in \C: \Im \xi > 0 \}$$ ($$\Im \xi$$ is the imaginary part of $$\xi$$). We denote the resolvent of $$\bA_n$$ by 
<p>
\[
R_{\bA_n}(z) = (\bA_n - z \id_n)^{-1},
\]
</p>
and denote the empirical eigenvalue distribution (spectral density) of $$\bA_n$$ by 
\\[
\mu_{\bA_n} = \frac{1}{n} \sum_{i=1}^n \delta_{\lambda_i(\bA_n)}. 
\\]


For a probability measure $$\mu \in \cP(\R)$$ on the real line, we denote its Stieltjes transform by
\\[
s_{\mu}(z) = \int \frac{1}{\lambda - z} \mu(\de \lambda). 
\\]
The next two lemmas characterize some properties of Stieltjes transforms. 

<div class="lemma"><ri>[Inverse Stieltjes transforms]</ri>
The scaled imaginary part of Stieltjes transform <script type="math/tex">(1/ \pi) \Im s_{\mu}(\cdot + \bi b)</script> is a probability measure on the real line <script type="math/tex">\R</script>. The sequence of probability measure <script type="math/tex">(1/ \pi) \Im s_{\mu}(\cdot + \bi b)</script> indexed by <script type="math/tex">b</script> converges weakly to <script type="math/tex">\mu</script> as <script type="math/tex">b \to 0+</script>. 
</div>

<div class="lemma"><ri>[Stieltjes continuity]</ri>
Let <script type="math/tex">\mu_n</script> be a sequence of random probability measures on the real line, and let <script type="math/tex">\mu</script> be a deterministic probability measure. Then <script type="math/tex">\mu_n</script> converges in expectation to <script type="math/tex">\mu</script> in the vague topology if and only if <script type="math/tex">\E s_{\mu_n}(z)</script> converges to <script type="math/tex">s_\mu(z)</script> for every <script type="math/tex">z</script> in the upper half-plane.
</div>



The next theorem is our focus in this blog. It gives the limiting Stieltjes transform of the GOE matrix. Combining with the two lemmas above, we can recover the semi-circle law of the GOE matrix. 


<div class="theorem"><ri>[Stieltjes transform of the GOE matrix]</ri>
Let <script type="math/tex">\bW_n \sim \GOE(n)</script>. That is, <script type="math/tex">\bW_n \in \R^{n \times n}</script> is a symmetric matrix with <script type="math/tex">(W_{n, ii})_{1 \le i \le n} \sim_{i.i.d.} \cN(0, 2/n)</script> and <script type="math/tex">(W_{n, ij})_{1 \le i \neq j \le n} \sim_{i.i.d.} \cN(0, 1/n)</script>. Denote by <script type="math/tex">s_{\bW_n} = s_{\mu_{\bW_n}}</script> the Stieltjes transform of empirical eigenvalue distribution of <script type = "math/tex">\bW_n</script>. Then we have 
\[
\lim_{n \to \infty}\E[s_{\bW_n}(z)] = \frac{-z + \sqrt{z^2 - 4}}{2}. 
\]
</div>

In the following sections, I will derive this result using the replica method. 


## 3. The determinant trick 

### 3.1. The connection of determinant and Stieltjes transform

Let $$\log$$ be the log function defined on $$\C$$ with cut on the negative real line. Define 
\\[
D_{\bW_n}(\xi) = \frac{1}{n} \sum_{i=1}^n \log(\lambda_i(\bW_n) - \xi),
\\]
By simple calculus we have 
\\[
\frac{\de}{\de \xi} D_{\bW_n}(\xi) = - \frac{1}{n} \sum_{i=1}^n \frac{1}{\lambda_i(\bW_n) - \xi} = - s_{\bW_n}(\xi), 
\\]
Note 
\\[
D_{\bW_n}(\xi) = \frac{1}{n}\log \det(\bW_n - \xi \id_n) + \frac{2 \pi k(\bW_n, \xi)}{n},
\\]
where $$k(\bW_n, \xi)$$ is an integer. Moreover, $$k(\bW_n, \xi)$$ will remain unchanged for an infinitesimal change of $$\xi$$. Hence we have 
\begin{align}
&~\lim_{n \to \infty} \E[s_{\bW_n}(\xi)] = - \lim_{n \to \infty} \frac{\de}{\de \xi} \E[D_{\bW_n}(\xi)] =  - \lim_{n \to \infty} \frac{\de}{\de \xi} \frac{1}{n} \E[\log \det(\bW_n - \xi \id_n)]\nonumber \\\
 \stackrel{\cdot}{=}&~ - \frac{\de}{\de \xi} \lim_{n \to \infty}  \frac{1}{n} \E[\log \det(\bW_n - \xi \id_n)]\tag{1}\label{eq:1}
\end{align}
In the last step, we heuristically exchanged the limit and the derivative operator. 

Why the determinant is easier to deal with than the Stieltjes transform? Because we can make use of the identity calculating the determinant
\\[
\det(\bA_n)^{-1/2} = \int_{\R^n} \frac{1}{(2 \pi)^{n/2}} \exp\Big\\{ - \frac{1}{2} \langle \bx, \bA_n \bx \rangle  \Big\\} \de \bx. 
\\]
This identity is correct when $$\bA_n$$ is positive semi-definite. However, we will use this identity formally when $$\bA_n$$ is a complex matrix. 


### 3.2. The replica approach

To calculate $$\E[s_{\bW_n}(\xi)]$$, we just need to calculate $$\lim_{n\to \infty}(1/n)\E[\log \det(\bW_n - \xi \id_n)]$$ and then take the derivative with respect to $$\xi$$. However, there is not a direct way to calculate the expectation of $$\log$$. One way to handle this difficulty is to use the replica formula introduced in my last blog
\\[
\E[\log Z] = \lim_{k \to 0} \frac{1}{k} \log \E[Z^k].
\\]
The replica formula reduce the problem to calculating the moments $$\E[\det(\bW_n - \xi \id_n)^k]$$, and then calculate the $$k \to 0$$ limit. 

### 3.3. An even easier approach

Instead of calculating $$\log \E[\det(\bW_n - \xi \id_n)^k]$$ for a sequence of $$k$$, we will just calculate $$\log \E[\det(\bW_n - \xi \id_n)^{-1/2}]$$. Note
\\[
\frac{1}{n} \log \det(\bW_n  - \xi \id_n ) = -\frac{2}{n} \log [\det(\bW_n  - \xi \id_n)^{-1/2}], 
\\]
and we expect that $$\log \det(\bW_n  - \xi \id_n )$$ tightly concentrate around its mean so that 
\begin{align}
&~\lim_{n \to \infty} \frac{1}{n} \E[ \log \det(\bW_n - \xi \id_n)] \stackrel{\cdot}{=} \lim_{n \to \infty} \frac{1}{n}\log \det(\bW_n - \xi \id_n)] \nonumber\\\
 \stackrel{\cdot}{=}&~ \lim_{n \to \infty} - \frac{2}{n}\log \det(\bW_n - \xi \id_n)^{-1/2}] \stackrel{\cdot}{=} \lim_{n \to \infty} -\frac{2}{n} \log \E[\det(\bW_n  - \xi \id_n)^{-1/2}]. \tag{2}\label{eq:2}
\end{align}
This is completely heuristic. 


## 4. The replica calculations

From now on, all the calculations are non-rigorous. We will perform formal computations following specific rules. When you want to use replica method yourself, follow the same rules and you can get some useful answers. 

We have formal identity (It is formal because $$\bW_n - \xi \id_n$$ is a complex matrix. Let us for now pretend that $$\xi$$ is a negative real number that is very large in absolute value, then the equality holds rigorous. )
<p>
\[
\E[\det(\bW_n - \xi \id_n)^{-1/2}] = \E\Big[ \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \exp\Big\{ - \bx^\sT (\bW_n - \xi \id_n) \bx / 2 \Big\} \de \bx \Big].
\]
</p>
Define $$\bG = (G_{ij})_{ij \in [n]} \in \R^{n \times n}$$ with $$G_{ij} \sim_{iid} \cN(0, 1)$$. Then $$\bW_n$$ and $$(\bG + \bG^\sT)/\sqrt{2n}$$ shares the same distribution. Hence, we get
\\[
\begin{aligned}
\E[\det(\bW_n - \xi \id_n)^{-1/2}] =&~ \E\Big[ \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \exp\Big\\{ - \Big\langle (\bG + \bG^\sT) / \sqrt {2n} - \xi \id_n, \bx \bx^\sT \Big\rangle / 2 \Big\\} \de \bx \Big]\\\
=&~ \E\Big[ \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \exp\Big\\{ - \Big\langle \bG , \bx \bx^\sT \Big\rangle / \sqrt {2n} + \xi \Vert \bx \Vert_2^2 /2 \Big\\} \de \bx \Big]\\\
=&~ \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \exp\Big\\{  \xi \Vert \bx \Vert_2^2 /2 \Big\\} \E\Big[\exp\Big\\{ - \Big\langle \bG , \bx \bx^\sT \Big\rangle / \sqrt {2n} \Big\\}\Big]  \de \bx. 
\end{aligned}
\\]
Using the formula of moment generating function of Gaussian variable $$\E_{G \sim \cN(0, 1)}[e^{tG}] = e^{t^2/2}$$, we get
\begin{align}
\E[\det(\bW_n - \xi \id_n)^{-1/2}] =&~ \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \exp\Big\\{  \xi \Vert \bx \Vert_2^2 /2 + \Vert \bx \bx^\sT \Vert_F^2 / (4n) \Big\\}  \de \bx \nonumber \\\
=&~ \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \exp\Big\\{  \xi \Vert \bx \Vert_2^2 /2 + \Vert \bx \Vert_2^4 / (4n) \Big\\}  \de \bx.  \tag{3}\label{eq:3} \\\
\end{align}
By now, the readers can completely forget about the rigor of our calculation, because the last integral in the equation above diverges for any real $$\xi$$, and is an oscillatory integral for any complex $$\xi$$. But we still continue our formal computations, hoping to extract something useful from this formal computations. 


Note we have 
\begin{align}\tag{4}\label{eq:4}
\int_\R \delta(\Vert \bx \Vert_2^2 - n s) \de s = 1. 
\end{align}
Plugging Eq. \eqref{eq:4} into Eq. \eqref{eq:3}, we get
\\[
\begin{aligned}
\E[\det(\bW_n - \xi \id_n)^{-1/2}] =&~ \int_{\R} \de s \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \exp\Big\\{  \xi \Vert \bx \Vert_2^2 /2 + \Vert \bx \Vert_2^4 / (4n) \Big\\} \cdot \delta(\Vert \bx \Vert_2^2 - n s)  \de \bx \\\
=&~ \int_{\R} \de s  \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \exp \\{  n \xi s /2 + n s^2 / 4 \\} \cdot \delta(\Vert \bx \Vert_2^2 - n s)  \de \bx \\\
=&~ \int_{\R}  \exp \\{ n \xi s /2 + n s^2 / 4 + n H_n(s) \\}  \de s,\\\
\stackrel{\cdot}{=}& \exp\Big\\{ n \cdot \ext_{s} \Big[  \frac{ \xi s}{2} + \frac{s^2}{4} + \lim_{n \to \infty} H_n(s) \Big] \Big\\},
\end{aligned}
\\]
where 
\\[
H_n(s) = \frac{1}{n} \log \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \cdot \delta(\Vert \bx \Vert_2^2 - n s)  \de \bx. 
\\]
We use the delta identity formula to carry on the calculations (let's pretend $$s$$ is a positive real number in this calculations; eventually we will take $$s$$ to be a complex number)
\\[
\begin{aligned}
H_n(s) =&~ \frac{1}{n}\log \int_{\R^n} \frac{1}{(2 \pi)^{n/2}} \cdot \delta(\Vert \bx \Vert_2^2 - n s)  \de \bx \\\
=&~ \inf_{\lambda} \frac{1}{n} \log \int_{\R^{n}} \frac{1}{(2 \pi)^{n/2}} \exp\\{ - \lambda \Vert \bx \Vert_2^2 / 2 + n s \lambda /2 \\} \de \bx\\\
=&~ \inf_{\lambda} \Big[ \frac{\lambda s}{2} - \frac{1}{2}\log \lambda \Big] \\\
=&~ \frac{1}{2} + \frac{1}{2} \log s. 
\end{aligned}
\\]
Therefore, we have
\\[
\E[\det(\bW_n - \xi \id_n)^{-1/2}] = \exp\Big\\{ n \cdot \ext_{s} \Big[  \frac{ \xi s}{2} + \frac{s^2}{4} + \frac{1}{2} \log s + \frac{1}{2} \Big] \Big\\}
\\]
and by Eq. \eqref{eq:2}, we have
\begin{align}
&~\lim_{n \to \infty} \frac{1}{n} \E[ \log \det(\bW_n - \xi \id_n)] \\\
=&~ - \lim_{n \to \infty} \frac{2}{n} \log \E[ \det(\bW_n - \xi \id_n)^{- 1/2}] \\\
=&~ -2 \cdot \ext_{s} \Big[  \frac{\xi s}{2} + \frac{s^2}{4} + \frac{1}{2} \log s + \frac{1}{2} \Big] \\\
=&~ \ext_{s} \Big[ - \xi s - \frac{s^2}{2} - \log s - 1 \Big] \\\
=&~ - \xi s_\star(\xi) - \frac{s_\star^2(\xi)}{2} - \log s_\star(\xi) - 1. 
\end{align}
where $$s_\star = s_\star(\xi)$$ satisfies 
\begin{align}\tag{5}\label{eq:5}
\xi + s_\star + \frac{1}{s_\star} = 0. 
\end{align}
By Eq. \eqref{eq:1}, we get 
\\[
\begin{aligned}
\E[s_{\bW_n}(\xi)] =&~ - \frac{\de}{\de \xi} \lim_{n \to \infty} \frac{1}{n} \E[ \log \det(\bW_n - \xi \id_n)] \\\
=&~ - \frac{\de }{\de \xi} \Big[\xi s_\star(\xi) - \frac{s_\star^2(\xi)}{2} - \log s_\star(\xi) - 1\Big] = - \frac{\partial }{\partial \xi} \Big[\xi s - \frac{s^2}{2} - \log s - 1\Big] \Big\vert_{s = s_\star(\xi)} \\\
=&~ s_\star(\xi). 
\end{aligned}
\\]
Let $$\xi \in \C_+$$, and finding the solution $$s_\star$$ of Eq. \eqref{eq:5} with positive imaginary part, we get 
\\[
s_\star(\xi) = \frac{- z + \sqrt{z^2 - 4}}{2}. 
\\]
This is exactly the limiting Stieltjes transform of GOE matrices. Of course, after we get this answer, we will find another way to derive it rigorously and sweep this heuristic calculations under our carpet. 

## 5. Summary 

In this blog, I used a lot of formal computations. I know there are some rules I can use in the computations, but I don't even know how these equations makes sense. But I know if I follow the formal computations, I will eventually get something useful. This is the mystery of replica method. I will summarize the rules to follow during the replica calculations in the next few blogs. Perhaps by then I will understand better why we can do these formal computations. 



