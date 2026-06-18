### Some math to warm up

**Q1)** Given a constant population of size 2*N* = 100, what is the
probability that two individuals pick the same parent one generation
before? What that they pick the same parent two generations before?
<details>
<summary>
Solutions
</summary>
The probability that the pick the same parent the previous generation is
$p(t_2=1) =\frac{1}{2N} = 0.01$. The probability that the pick the same
parent two generations ago implies that they did not do so one
generation ago. Hence
$P(t_2=2)=\left(1-\frac{1}{2N}\right)\frac{1}{2N} = 0.0099$.
</details>

<br/>

**Q2)** A researcher sequences a DNA fragment of length 10kb of a
diploid individual and detects 21 heterozygous sites. Assume that
*μ* = 10<sup>−9</sup> per site for this segment and that the assumptions
of the standard coalescent and the infinites sites model hold.

**a)** What is the effective size of this population?

<details>
<summary>
Solutions
</summary>
Under the infinite sites model, E \[*H*\] = 4*N*<sub>*e*</sub>*μ*, and
hence we can get an estimate $\hat{N_e}=\frac{H}{4\mu}$. Here,
$H=\frac{21}{10^4}=2.1$, and hence
$\hat{N_e}=\frac{2.1\cdot 10^{-3}}{4\cdot 10^{-9}}=\frac{2.1}{4}10^6 = 5.25 \cdot 10^5$.}
</details>

<br/>

**b)** What is the expected time to the most recent common ancestor of
the two sequences in this individual?

<details>
<summary>
Solutions
</summary>
The
E \[*T*<sub>*M**R**C**A*</sub>\] = 2*N*<sub>*e*</sub> = 1.05 ⋅ 10<sup>6</sup>
generations.}
</details>

<br/>

**Q3)** Consider a sample of two diploid individuals from a population
for which the assumptions of the standard coalescence model holds true.
What is the probability that the two sequences of the first individual
are more closely related to each other than they are to any of the
sequences of the other individual?

<details>
<summary>
Solutions
</summary>
There are two cases (topologies) for which this is true:<br/> Case 1: If
the two lineages from the first individual are the first to coalesce,
which happens with probability
$\frac{1}{\binom{4}{2}} = \frac{1}{6}$.<br/> Case 2: If the lineages of
the other individual coalesce first, and then the two lineages of the
first individual coalesce. The probability that the lineages of the
other individual coalesce first is also $\frac{1}{6}$. The probability
that the second coalescent event is then between the lineages of the
first individual is $\frac{1}{\binom{3}{2}}=\frac{1}{3}$. The total
probabilty of this case is thus
$\frac{1}{6} \cdot \frac{1}{3} = \frac{1}{18}$.<br/> Hence, the total
probability that the lineages of the first individual are more closely
related than any of them is to a sequence of the other individual is
$\frac{1}{6} + \frac{1}{18} =   \frac{3}{18} + \frac{1}{18} = \frac{2}{9}$.
</details>

<br/>

**Q4)** The following DNA sequences were obtained from a population:

<pre><code>A<span style="color:red;">T</span>CGTGCAC<span style="color:red;">A</span> AC<span style="color:red;">T</span>TGCAACA</code>
<code>A<span style="color:red;">T</span>CGTG<span style="color:red;">G</span>ACC AC<span style="color:red;">T</span>TGCAAC<span style="color:red;">T</span></code>
<code>AGCGTG<span style="color:red;">G</span>ACC AC<span style="color:red;">T</span>TGCAAC<span style="color:red;">T</span></code></pre>

**a)** Calculate the Tajima estimator *θ*<sub>*T*</sub> per site from
this data.

<details>
<summary>
Solutions
</summary>
Pairwise differences *d*<sub>*i**j*</sub> between individuals *i* and
*j*: *d*<sub>12</sub> = 3, *d*<sub>13</sub> = 4, *d*<sub>14</sub> = 3,
*d*<sub>23</sub> = 3, *d*<sub>24</sub> = 4 and
*d*<sub>34</sub> = 1.<br/> The average number of pair-wise differences
is
$\pi = \frac{d\_{12}+d\_{13}+d\_{14}+d\_{23}+d\_{24}+d\_{34}}{6} = \frac{18}{6} = 3$<br/>
Hence, the Tajima estimator is
$\theta_T = \frac{\pi}{20}= \frac{3}{20} = 0.15$ per base pair.
</details>

<br/>

**b)** Calculate the Watterson estimator *θ*<sub>*w*</sub> per site from
this data.

<details>
<summary>
Solutions
</summary>
There are a total of *S* = 5 segregating sites in this data.<br/>
$\theta_W = S \frac{1}{\displaystyle \sum\_{k=1}^{n-1}\frac{1}{k}} = \frac{5}{1 + \frac{1}{2} + \frac{1}{3}} = 5\frac{6}{11} = \frac{30}{11} = 2.7$,
or *θ*<sub>*W*</sub> = 0.14 per base pair.
</details>

<br/>

**c)** Why are these two estimators not identical?

<details>
<summary>
Solutions
</summary>
Since they are both estimators, they are often slightly different even
if the population was constant and evolving under the assumptions of
Wright-Fisher model. Only if there was an infinite amount of data, the
two must be identical. Importantly, however, these estimates will likely
be different if the assumption of a standard neutral coalescent is
violated.
</details>

<br/>

**d)** Plot a histogram of the site frequency specturm (SFS) for this
data. Do the data contain more, fewer or the same number of singletons
as expected under the standard neutral coalescent model?

<details>
<summary>
Solutions
</summary>

``` r
barplot(c(1,3,1), xlab = "Derived allele frequency", ylab = "Counts")
```

![](exercises_day4_afternoon_coalescent_files/figure-markdown_github/unnamed-chunk-1-1.png)
The expected number of mutations with frequency *k* is given by
$\operatorname{E}\[f_k\]=\frac{\theta}{k}$. Hence, we expect
$\E\[f_1\]=\theta$ singletons per base pair. Among the 20 base pairs, we
thus expect between 20*θ*<sub>*W*</sub> = 2.8 and
20*θ*<sub>*T*</sub> = 3 singletons. Thus, it seems we have too few
singletons in the data. In contrast, we see too many doubltons: We
expect $\operatorname{E}\[f_2\]=\frac{\theta}{2}$ doubletons per base
pair, which woudl correspond to 10*θ*<sub>*W*</sub> = 1.4 and
10*θ*<sub>*T*</sub> = 1.5 doubletons, not 3. Or to put differently: we
dould expect twice as many singletons that doubletons. However, there is
so little data that we expect differences between expectations and
observations. The only way to tell if there are too few singletons is
thus a statistical test.
</details>

<br/>

### Installing `fastsimcoal2`

In the following exercises you will simulate genealogies using the
coalescent simulator `fastsimcoal2`. You can download the latest (Linux)
version of `fastsimcoal2` as follows (make sure to be in the directory
in which you would like to have it):

``` bash
wget https://cmpg.unibe.ch/software/fastsimcoal2/downloads/fsc28_linux64.zip
unzip fsc28_linux64.zip
```

The `fastsimcoal2` software comes with a lot of extra files (manual,
example files, …). Let’s move the executable to the current directory
and make sure it is executable. When doing so, we will also rename the
binary so the exercises below work regardless of the version you have.

``` bash
mv fsc2[0-9]*_linux64/fsc2[0-9]* fastsimcoal2
chmod +x fastsimcoal2
```

You can easily test if that worked by trying to launch `fastsimcoal2`

``` bash
./fastsimcoal2
```

In the following, you may want to work in different folders to keep your
work organized. Remember to always specify the absolute or relative path
to your executable to make sure the solutions work. Alternatively, you
can also add the directory with your executable to your `PATH` variable
so you can use it from any location (needs to be redone if you reconnect
to the system):

``` bash
PATH="`pwd`:${PATH}"
```

### Generating your first simulations

Let’s begin by simulating data under a constant size population and
thereby learn how to use `fastsimcoal2`. The demographic model you want
to simulate needs to be specified in an in an input file. Here is an
example of such a file for a constant size population:

    //Parameters for the coalescence simulation program : simcoal.exe
    1 samples to simulate
    //Deme sizes (haploid number of genes)
    10000
    //Sample sizes
    20
    //Growth rates
    0
    //Number of migration matrices : If 0 : No migration between demes
    0
    //Historical event: time, source, sink, migrants, new deme size, new growth rate, new migration matrix
    0 historical events
    //Number of independent chromosome
    1 0
    //Number of contiguous linkage blocks
    1
    //Per Block: Data type, No. of loci, Recombination rate to the right-side locus, plus optional parameters
    DNA 10000 0 1.0e-07 0.5

Have a look at that file. It first specified the number of populations
(samples), then the initial (backward in time) haploid population size
*N* (NOT 2*N*!), the number of haploid lineages (Sample size) to sample
from that population and the initial growth rate of the population.
After that, it lists the number of migration matrices (we do not model
migration, hence 0), any historical events (we do not habe any events)
and finally the genetic markers to simulate: 1 chromosome with one fully
linked block of 10kb with a mutation rate of 10<sup>−7</sup>.

You may also consult the manual of fastsimcoal2, which is available
[here](http://cmpg.unibe.ch/software/fastsimcoal2/man/fastsimcoal28.pdf).

Save the above content in a file `constsize.par`. You can then use this
file to generate simulations as follows:

``` bash
fastsimcoal2 -i constsize.par -n 20 -T
```

Here, `-n 20` tells `fastsimcoal2` to generate 20 replicate simulations
and `-T` to also write the simulated genealogies as we want to look at
them.

Have a look at the output, which you find in a folder with the same name
as the input file, so `constasize`. Inside the folder, you find files
with the extension `.arp` that contain the generated genetic data (have
a look at one file) as well as two files containing the trees:

1.  `constsize_1_true_trees.trees` with the true trees, i.e. with branch
    length in generations.

2.  `constsize_1_mut_trees.trees` with the mutation trees, i.e. with
    branch length in the number of mutations that fell onto them.

Let’s visualize these genealogies. You can do so in R using the package
`phytools`. Here is a plotting function that arranges the trees nicely
so they can be compared:

``` r
plotTrees <- function(trees, popCol=c("black", "orange2", "purple")){
  maxHeight <- max(unlist(lapply(trees, function(x){ max(branching.times(x)) })))
  nCols <- ceiling(sqrt(length(trees)))
  nRows <- ceiling(length(trees) / nCols)
  par(mfrow = c(nRows, nCols), oma=c(0,4,0,0), las=1, xpd=NA)
  
  for(tr in 1:length(trees)){
    plotTree(trees[[tr]], direction="downwards", ylim=c(0, maxHeight),
             lwd=0.6, ftype="off", mar=c(0.1,0.7,0.1,0.7))
    
    if(tr %% nCols == 1){ axis(side = 2)}
    
    #add tips and color by population
    pop <- as.numeric(unlist(lapply(strsplit(trees[[tr]]$tip.label, "[.]"), '[', 2) ))
    
    nTips <- length(trees[[tr]]$tip.label);
    symbols(1:nTips, rep(0, nTips), circles=rep(0.4, nTips), add=TRUE, inches=0.03, fg=NA, bg=popCol[pop])
  }
}
```

You can use this function on the trees you just simulates as follows
(make sure the directory is correct):

``` r
library(phytools)
trees <- read.nexus("constsize/constsize_1_true_trees.trees")
plotTrees(trees)
```

**Q5)** Use the above setup to simulate and compare trees under constant
size populations of different sizes. Always compute the average tree
height using

``` r
mean(unlist(lapply(trees, function(x){ max(branching.times(x)) })))
```

Does it match the expectation?

<details>
<summary>
Solutions
</summary>
The topologies are independent of the population size and should not be
visually different between simulations conducted with a small or large
population size. The expected average tree height, however, should be
affected by population size. It is 4*N* where *N* is the number of
diploid individuals. The `Deme size` specified in `fastsimocal2`,
however, is in the number of haploid individuals and hence twice that of
the number of diploid individuals. The expected tree height should
therefore match twice the number of haploid individuals.
</details>

<br/>

### Simulating an SFS from many loci

Let us next use `fastsimcoal2` to simulate site frequency spectra. We do
so by asking `fastsimcoal2` to simulate a larger number of independent
loci (DNA sequences) and to then tabulate the frequencies of all
polymorphic sites. Specifically, copy the above input file to a new file
`constsize_sfs.par` and change the following line:

    //Number of independent chromosome
    10000 0

This tells `fastsimcoal2` to simulate 10<sup>4</sup> independent loci
(each as specified in the lined that follow, i.e. of 10kb each). To then
generate the SFS call `fastsimcoal2` using the options `-s0` (to tell it
to use all polymorphic sites) and `-d` to tell it to write the derived
SFS:

``` bash
./fastsimcoal2 -i constsize_sfs.par -n1 -s0 -d
```

**Q6)** Have a look at the simulated SFS (found in the file
`constsize_sfs_DAFpop0.obs`) and understand its format. What is the
fraction of singletons to doubletons?

<details>
<summary>
Solutions
</summary>
The singletons and doubletons are shown in the second and third column
with header `d0_1` and `d0_2`. Under a constant size population, there
should be twice as many singletons as doubletons.
</details>

<br/>

### Simulating exponential population growth

Simulate a model of exponential population growth. You can do so my
creating a new input file `expgrowth.par` that differs from your
constant size file on these lines:

    //Growth rates
    -0.01

Note that we specify the growth rate as negative because `fastsimcoal2`
is running backward in time - a growing population is shrinking backward
in time.

**Q7)** Do the genealogies look different compared to the constant size
case? And how does the SFS change?

<details>
<summary>
Solutions
</summary>
You should observe that the genealogies under population growth have
much longer terminal branches. That also results in an excess of rare
compared to common variants. There are, for instance, way more than
twice as many singletons than doubletons.
</details>

<br/>

### Simulating population splits

**Q8)** Let us next simulate a model of a simple population split that
occurred 100 generations ago. For this, again copy the original
constsize input file to `split.par` and edit it as follows:

1.  Specify that you want to simulate 2 samples.
2.  Provide an additional deme size (1,000), sample size (use 10 for
    both) and growth rate (0).
3.  Add a historical event (set the count to 1) and specify the event as
    follows:

<!-- -->

    //Historical event: time, source, sink, migrants, new deme size, new growth rate, new migration matrix
    1 historical events
    10000 1 0 1 1 0 0

Here, the first number indicates the time (10,000 generations) of the
event. The next numbers indicate the source and sink population and the
fraction of migrants to move. Here we specify that we want to move all
(1.0) of all individuals (migrants) from population 1 to population 0
(note that `fastsimcoal2` starts counting populations at 0). Thus, to
simulate a split, we actually simulate two populations and at the time
of the split we move the samples from one population to the other (as we
are going backward in time). The last two number specify the new deme
size relative to the previous one (1, i.e. no change) and the new
migration matrix (0, i.e. still no migration).

Note that `fastsimcoal2` is super picky regarding the settings file and
often fails to give proper errors if something is off. So if you get
strange errors or the model seems off, check your input file carefully!

<details>
<summary>
Solutions
</summary>

    //Parameters for the coalescence simulation program : simcoal.exe
    2 samples to simulate
    //Deme sizes (haploid number of genes)
    10000
    1000
    //Sample sizes
    10
    10
    //Growth rates
    0
    0
    //Number of migration matrices : If 0 : No migration between demes
    0
    //Historical event: time, source, sink, migrants, new deme size, new growth rate, new migration matrix
    1 historical events
    10000 1 0 1 1 0 0
    //Number of independent chromosome
    1 0
    //Number of contiguous linkage blocks
    1
    //Per Block: Data type, No. of loci, Recombination rate to the right-side locus, plus optional parameters
    DNA 10000 0 1.0e-07 0.5

</details>

<br/>

**Q9)** Run this model and plot some genealogies (the function
`plotTrees()` will use different colors for the samples from different
populations). Are the populations well separated? And can you recognize
which population was smaller in size? How do the trees differ if the
split time was much younger, say 100 generations?

<details>
<summary>
Solutions
</summary>
With a split time of 10,000 generations, the samples from population 1
will mostly have coalesced prior to the split (the expected time to the
MRCA is 2,000 generations). So the samples from that population should
mostly form a monophyletic cluster and reach their MRCA earlier than the
samples of population 0. But do appreciate the variation in
genealogies!<br/> If you use a much younger split time, then the samples
within each population will not have time to coalesce prior to the split
and appear mingled in the genealogies.
</details>

<br/>

**Q10)** Now run it with a split time of 50,000 generations. What is the
average time (across replicate genealogies) to the MRCA of all samples?
Does it match the expectation?

<details>
<summary>
Solutions
</summary>
The expected time to the MRCA should be the split time plus the time to
the MRCA within the ancestral population. With a split time of 50,000
generations, there is likely only one lineage per population surviving
until the split. The time to MRCA in the ancestral population (which is
if size 10,000) will thus be 10,000 generations and the expected total
time 50,000 + 10,000 = 60,000 generations.
</details>

<br/>

### Simulating a bottleneck

**Q11)** Let us finally simulate a strong bottleneck. For this, start
again with your initial `constsize.par` file and modify it by adding two
historical events (remember that `fastsimcoal2` works backward in time):

1.  At 90 generation in the past, the population size is reduced to 10
    individuals.
2.  At 100 generations ago, the population size recovers to the original
    size of 10,000 individuals.

<details>
<summary>
Solutions
</summary>

    //Parameters for the coalescence simulation program : simcoal.exe
    1 samples to simulate
    //Deme sizes (haploid number of genes)
    10000
    //Sample sizes
    20
    //Growth rates
    0
    //Number of migration matrices : If 0 : No migration between demes
    0
    //Historical event: time, source, sink, migrants, new deme size, new growth rate, new migration matrix
    2 historical events
    90 0 0 1 0.001 0 0
    100 0 0 1 1000 0 0
    //Number of independent chromosome
    1 0
    //Number of contiguous linkage blocks
    1
    //Per Block: Data type, No. of loci, Recombination rate to the right-side locus, plus optional parameters
    DNA 10000 0 1.0e-07 0.5

</details>

<br/>

**Q12)** Run this model an look at the genealogies and compare them to
those obtained for a constant size population. What do you notice?

<details>
<summary>
Solutions
</summary>
During the bottleneck, many coalescent events occur. This results in
very many terminal branches that are very short.
</details>

<br/>

**Q13)** How do you think this will affect the SFS? Confirm your
intuition with simulations!

<details>
<summary>
Solutions
</summary>
Check above to see how to simulate the expected SFS. Since the terminal
branches are shorter than under a constant size model, there are too few
rare mutations compared to the common mutations. Under a constant size
model, we expect 1/5 as many 5-tons than singletons, for instance. But
under this severe bottleneck model, we obtain almost half as many as
singletons.
</details>

<br/>
