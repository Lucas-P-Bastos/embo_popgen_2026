### Transition probabilities under the Wright-Fisher model

**Q1)**. Let’s warm up by refreshing the Wright-Fisher model.

**a)** What is the average number of offspring per individual?

<details>
<summary>
Solutions
</summary>
Since the population size is constant, each individuals has on average 1
offspring.
</details>

<br/>

**b)** What is the probability that an individual has no offspring at
all? Does that depend on the population size?

<details>
<summary>
Solutions
</summary>
This is happening when each of the 2N individuals of the next generation
picks another individual as its parent. Hence,
$P(\mathrm{no\\offspring}) = \left(1-\frac{1}{2N}\right)^{2N}$, which
does depend on the population size but becomes approximately 0.37 if
2*N* is large.
</details>

<br/>

**Q2)** Consider a population of size 2*N* = 100 and a current allele
frequency of *f* = 10/100 = 0.1.

**a)** What is the probability that the allele frequency remains
unchanged after one generation of random mating?

<details>
<summary>
Solutions
</summary>

Under the Wright-Fisher Model, transition probabilities are binomial.
Hence, $P(n\|2N, f)=\choose{2N}{n}f^n(1-f)^{2N})$. In R, we can get that
probability using the function `dbinom()`:

``` r
dbinom(x=10, size = 100, prob = 0.1)
```

\[1\] 0.1318653
</details>

<br/>

**b)** What is the probability that the allele is lost after one
generation of random mating?

<details>
<summary>
Solutions
</summary>

Again using `dbinom()`:

``` r
dbinom(x=0, size = 100, prob = 0.1)
```

\[1\] 2.65614e-05
</details>

<br/>

**c)** What is the probability that the allele is fixed after one
generation of random mating?

<details>
<summary>
Solutions
</summary>

Again using `dbinom()`:

``` r
dbinom(x=100, size = 100, prob = 0.1)
```

    ## [1] 1e-100

</details>

<br/>

**d)** What is the probability that the allele frequency decreases after
one generation of random mating?

<details>
<summary>
Solutions
</summary>

To calculate the probability of a decrease, we need to integrate the
transition probability from 0 to n-1. In R, we can use the function
`pbinom()`:

``` r
pbinom(9, size = 100, prob = 0.1)
```

    ## [1] 0.4512902

</details>

<br/>

**d)** What is the probability that the allele frequency increases after
one generation of random mating?

<details>
<summary>
Solutions
</summary>

We can calculate the probability of increase analogously as the integral
from 201 to 1000, or as 1 minus the integral from 0 to 200:

``` r
1 - pbinom(10, size = 100, prob = 0.1)
```

    ## [1] 0.4168445

Note that the probabilities for an increase and decrease are not equal,
i.e. this distribution is not symmetric.
</details>

<br/>

**Q3)** Plot the probability distribution on the allele frequency for a
current allele frequency *f* = 0.1 and 2*N* = 10, 2*N* = 100,
2*N* = 1000 and 2*N* = 10000. What is the effect of population size?

<details>
<summary>
Solutions
</summary>

``` r
par(mfrow=c(2,2))
for(twoN in c(10, 100, 1000, 10000)){
  p <- 0:twoN
  d <- dbinom(p, size = twoN, prob = 0.1)
  plot(p, d, type = 'b')
}
```

![](exercises_popgen_files/figure-markdown_github/unnamed-chunk-6-1.png)
</details>

<br/>

### Simulating neutral allele frequency trajectories

**Q4)** Write a function `simulateWF()` to simulate allele trajectories
under the Wright-Fisher model. Your function should take as input i) the
population size 2*N*, ii) the initial allele frequency *f* and iii) the
number *G* of generations to simulate. It should then return the allele
frequency (between 0 and 1) for each generation as a vector.

<details>
<summary>
Solutions
</summary>

We can use the function `rbinom()` to simulate trajectories:

``` r
simulateWF <- function(twoN, f, G){
  p <- numeric(G + 1)
  p[1] <- f
  for(i in 1:G){
    p[i+1] <- rbinom(1, size = twoN, prob = p[i]) / twoN
  }
  return(p)
}
```

</details>

<br/>

**Q5)** Use your function `simulateWF()` to simulate 1000 trajectories
with 2*N* = 100 and *f* = 0.1 for *G* = 1000 and plot them in one plot.
In how many cases was the allele lost? Does this match the expectation?
Repeat for different population sizes and initial allele frequencies.
What is the effect of the population size?

<details>
<summary>
Solutions
</summary>

``` r
trajectories <- replicate(1000, simulateWF(twoN = 100, f = 0.1, G = 1000))
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))
```

![](exercises_popgen_files/figure-markdown_github/unnamed-chunk-8-1.png)

``` r
print(paste("Allele was lost in", sum(trajectories[1000,] == 0), "/", ncol(trajectories), "cases."))
```

    ## [1] "Allele was lost in 897 / 1000 cases."

The allele is expected to go to fixation with the probability of its
initial allele frequency *f*, and hence to be lost with probability
1 − *f*.

If you use larger population sizes, you should observe less fluctuation
and eventually many sites will remain polymorphic even after *G* = 1000
generations. However, if run for enough generations, the same fraction
of alleles will be lost.
</details>

<br/>

**Q6)** Use your function `simulateWF()` to study fixation probability
of a new mutation under different population sizes (ensure *G* is large
enough). Does the fixation probability depend on the population size?
How does this affect the substitution rate?

<details>
<summary>
Solutions
</summary>

To study the fixation probability of a new mutation, set
$f=\frac{1}{2N}$.

``` r
trajectories <- replicate(1000, simulateWF(twoN = 100, f = 1/twoN, G = 1000))
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))
```

![](exercises_popgen_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
print(paste("Allele was fixed in", sum(trajectories[1000,] == 1), "/", ncol(trajectories), "cases."))
```

    ## [1] "Allele was fixed in 0 / 1000 cases."

The fixation probability of any mutation is given by its current
frequency *f* and is hence markedly different between small and large
populations. The substitution rate, hower, is not since there are
2*N**μ* mutations occuring per generation, and hence many more in large
than small populations.
</details>

<br/>

### Simulating selection and genetic drift

**Q7)** Write a function `simulateWFWithSelection()` that simulates
allele trajectories under both genetic drift and viability selection.
Similar to your function `simulateWF()` you wrote above, it should take
as input i) the population size 2*N*, ii) the initial allele frequency
*f*, iii) the number *G* of generations to simulate and iv) also a
vector *v* of viabilities for the genotypes AA, Aa and aa. It should
then return the allele frequency (between 0 and 1) for each generation
as a vector. In each generation, your function should apply selection to
alter the allele frequency, and then use binomial sampling to simulate
genetic drift in that modified allele frequency.

<details>
<summary>
Solutions
</summary>

We can use the function `rbinom()` to simulate trajectories:

``` r
simulateWFWithSelection <- function(twoN, f, G, v){
  p <- numeric(G + 1)
  p[1] <- f
  for(i in 1:G){
    # selection
    fA <- p[i]
    fa <- 1-fA
    fPrime <- (v[1]*fA*fA + v[2]*fA*fa)/(v[1]*fA*fA + v[2]*2*fA*fa + v[3]*fa*fa);
    # drift
    p[i+1] <- rbinom(1, size = twoN, prob = fPrime) / twoN
  }
  return(p)
}
```

</details>

<br/>

**Q8)** Use your function `simulateWFWithSelection()` to simulate 1000
trajectories with 2*N* = 100, and *f* = 0.1 and genic selection with
viabilities *v* = (1, 1 − *s*, (1 − *s*)<sup>2</sup>) for *G* = 1000 and
plot them in one plot. Start with *s* = 0.01. How strong does the
selection coefficient *s* have to be to see differences to the neutral
case simulated above (or simulated with `simulateWFWithSelection()` when
setting *v* = (1, 1, 1))? How is this affected by the population size
2*N*?

<details>
<summary>
Solutions
</summary>

``` r
s <- 0.01
trajectories <- replicate(100, simulateWFWithSelection(twoN = 1000, f = 0.1, G = 1000, v=c(1,1-s,(1-s)^2)))
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))
```

![](exercises_popgen_files/figure-markdown_github/unnamed-chunk-11-1.png)

``` r
print(paste("Allele was fixed in", sum(trajectories[1000,] == 1), "/", ncol(trajectories), "cases."))
```

    ## [1] "Allele was fixed in 86 / 100 cases."

A weak selection coefficient of *s* = 0.01 leads to only very marginal
difference to the neutral case (allele is lost in about 80% rather than
90% of the cases). With *s* = 0.1, the allele goes to fixation in the
majority of cases. When using 2*N* = 1000, the allele goes to fixation
in about 90% of the cases even with *s* = 0.01.

</details>

<br/>

**Q9)** Set 2*N* = 10<sup>6</sup> and dominant *v* = (1, 1, 1 − *s*)
with *s* = 0.05. How often does the allele go to fixation within
*G* = 1000? And in case of a smaller population size? How does it look
like for the recessive case?

<details>
<summary>
Solutions
</summary>

``` r
s <- 0.05
trajectories <- replicate(100, simulateWFWithSelection(twoN = 10^6, f = 0.1, G = 1000, v=c(1,1,1-s)))
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))
```

![](exercises_popgen_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
print(paste("Allele was fixed in", sum(trajectories[1000,] == 1), "/", ncol(trajectories), "cases."))
```

    ## [1] "Allele was fixed in 0 / 100 cases."

Despite a large population 2*N* = 10<sup>6</sup> and strong selection
*s* = 0.05, the allele is essentially never fixed. This is because there
is no selection benefiting the heterozygous Aa over the homozygous AA
genotype and drift is too weak to push the allele to fixation. When
using a smaller population size such as 2*N* = 10<sup>3</sup> the allele
is essentially always fixed because selection is strong enough to push
it quickly to high frequencies, and drift strong enough to fix it by
random fluctuations.

In the recessive case *v* = *c*(1, 1 − *s*, 1 − *s*) with *s* = 0.05 and
2*N* = 10<sup>6</sup>, the allele is always fixed. When using a smaller
population size (e.g. 2*N* = 1<sup>3</sup>), the allele is occasionally
lost early on because of drift (random fluctuations) before selection
could push the allele to a high enough frequency, also because there is
no selection benefiting the heterozygous Aa over the homozygous aa
genotype.
</details>

<br/>
