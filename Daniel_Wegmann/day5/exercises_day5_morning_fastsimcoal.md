## Inferring Demography with `fastsimcoal2`

### Generating some observed data

The goal of this first set of exercises is to learn how to do
demographic inference with `fastsimcoal2`. We will start with a very
question: can we date a population size change?

**Q1)** Write an input file `bottleneck_obs.par` for `fastsimcoal2` to
simulate an SFS (10,000 independent loci) under a single population of
size 10,000 in the past but that was reduced to a size of only 500
individuals 100 generations ago. We use the tag `obs` so we do not
forget that this file was used to generate the observed data. Do you
remember how to write such an input file?

<details>
<summary>
Solutions
</summary>

Your input file should look like this:

    //Parameters for the coalescence simulation program : simcoal.exe
    1 samples to simulate
    //Deme sizes (haploid number of genes)
    500
    //Sample sizes
    20
    //Growth rates
    0
    //Number of migration matrices : If 0 : No migration between demes
    0
    //Historical event: time, source, sink, migrants, new deme size, new growth rate, new migration matrix
    1 historical events
    100 0 0 1 20 0 0
    //Number of independent chromosome
    10000 0
    //Number of contiguous linkage blocks
    1
    //Per Block: Data type, No. of loci, Recombination rate to the right-side locus, plus optional parameters
    DNA 10000 0 1.0e-07 0.5

</details>

<br/>

**Q2)** Use your input file to simulate an SFS under this model. Do you
remember the command line?

<details>
<summary>
Solutions
</summary>

``` bash
./fastsimcoal2 -i bottleneck_obs.par -n1 -s0 -d
```

</details>

<br/>

### Your first inference with `fastsimcoal2`

Let us now pretend that we do not know the time of the bottleneck this
population experienced, but that we wish to infer it from the SFS. To do
so, we need to prepare three files:

**File 1**: The file `bottleneck.tpl` which specified the demographic
model we wish to infer. It differs from your `bottleneck_obs.par` file
in two ways:

1.  It uses tags rather than values for all parameters you wish to
    infer. Since we aim to infer the timing of the bottleneck, replace
    that value (100) with the tag `T_BOTTLENECK` (you may use any tag,
    this is just a suggestion).

2.  It uses a special data type `FREQ` so the expected SFS is outputted.
    So change the part of the input file regarding the genetic loci to

<!-- -->

    //Number of independent chromosome
    1 0
    //Number of contiguous linkage blocks
    1
    //Per Block: Data type, No. of loci, Recombination rate to the right-side locus, plus optional parameters
    FREQ 1 0 1.0e-07

<details>
<summary>
Solutions
</summary>

Your `bottleneck.tpl` file should look like this:

    //Parameters for the coalescence simulation program : simcoal.exe
    1 samples to simulate
    //Deme sizes (haploid number of genes)
    500
    //Sample sizes
    20
    //Growth rates
    0
    //Number of migration matrices : If 0 : No migration between demes
    0
    //Historical event: time, source, sink, migrants, new deme size, new growth rate, new migration matrix
    1 historical events
    T_BOTTLENECK 0 0 1 20 0 0
    //Number of independent chromosome
    1 0
    //Number of contiguous linkage blocks
    1
    //Per Block: Data type, No. of loci, Recombination rate to the right-side locus, plus optional parameters
    FREQ 1 0 1.0e-07

</details>

<br/>

**File 2**: The file `bottleneck.est` which specifies the parameters to
be learned. It should look like this (make sure there is a new line
after the “output” tag):

    // Parameter file for the bottleneck model
    // ***************************************
    [PARAMETERS]
    //#isInt? #name #dist. #min #max
    1 T_BOTTLENECK unif 1 10000 output

This file defines all parameters that need to be learned on one line
under `[PARAMETERS]`. Here, we only want to infer `T_BOTTLENECK`, so it
only contains that one line. The first number (1) specifies that this
parameter is an integer (a 0 would mean that it is a floating point
number). This is followed by the name of the parameter, the “prior”
distribution (does not matter for the MLE inference we perform), the
minimum and maximum value, and the tag `output` to specify that this
parameter should be reported.

**File 3**: The file with the obsered SFS. `fastsimcoal2` expects this
to have a very specific name, namely `bottleneck_DAFpop0.obs`, which
consists of the name of the `tpl` file followed by `_DAFpop0.obs` as
this is the derived SFS. It also expects that file to be in the same
folder as your `tpl` file - so you need to move it there.

You can then launch an inference as follows:

``` bash
./fastsimcoal2 -M -t bottleneck.tpl -e bottleneck.est -d -n100000  -L 25 -c 4 -q
```

Here, `-M` specifies that we wish to conduct parameter inference, `-t`
and `-e` are used to specify the `tpl` and `est` file used for the
inference, `-d` specifies that we provide a derived SFS (the name of
which is hard-coded), `-n100000` that 100,000 genealogies shall be
simulated to approximate the expected SFS (larger numbers = more
accurate = slower), `-M` that the maximum likelihood shall be calculated
and printed for each iteration, `-L 25` that 25 optimization rounds
shall be run, `-c4` that 4 cores are to be used and `-q` that minimal
output should be written. (run `./fastsimcoal2` to see all possible
options).

If you run the above command line, `fastsimcoal2` will print the current
parameter estimate and the likelihood for each iteration. Is the final
estimate close to the true parameter value you used for the simulation?

The run will also produce three files inside the `bottleneck` folder:

1.  The file `bottleneck.bestlhoods` that contains the final estimate.
2.  The file `bottleneck.brent_lhoods` that contains the estimates in
    each iteration.
3.  The file `bottleneck_maxL.par` that contains a version of the `tpl`
    file with the best parameter values inserted instead of the tags.

**Q3)** Modify the files `bottleneck.tpl` and `bottleneck.est` such that
also the population size change gets estimated (i.e. use the tag
`SIZE_CHANGE` instead of the factor 20). Can these parameters both be
inferred?

<details>
<summary>
Solutions
</summary>

Your `bottleneck.tpl` file should look like this:

    //Parameters for the coalescence simulation program : simcoal.exe
    1 samples to simulate
    //Deme sizes (haploid number of genes)
    500
    //Sample sizes
    20
    //Growth rates
    0
    //Number of migration matrices : If 0 : No migration between demes
    0
    //Historical event: time, source, sink, migrants, new deme size, new growth rate, new migration matrix
    1 historical events
    T_BOTTLENECK 0 0 1 SIZE_CHANGE 0 0
    //Number of independent chromosome
    1 0
    //Number of contiguous linkage blocks
    1
    //Per Block: Data type, No. of loci, Recombination rate to the right-side locus, plus optional parameters
    FREQ 1 0 1.0e-07

Your `bottlneck.est` file should look like this:

    // Parameter file for the bottleneck model
    // ***************************************
    [PARAMETERS]
    //#isInt? #name #dist. #min #max
    1 T_BOTTLENECK unif 1 10000 output
    0 SIZE_CHANGE unif 0 50 output

With these files, estimates should run very well!
</details>

<br/>

### Challenging the optimizer of `fastsimcoal2`

**Q4)**\* Let us next illustrate that optimizing demographic parameters
is difficult by also inferring the current size. Modify the files
`bottleneck.tpl` and `bottleneck.est` accordingly (i.e. use the tag
`N_CURRENT` instead of the fixed size 500) and rerun the estimation.
What do you observe?

<details>
<summary>
Solutions
</summary>

Your `bottleneck.tpl` file should look like this:

    //Parameters for the coalescence simulation program : simcoal.exe
    1 samples to simulate
    //Deme sizes (haploid number of genes)
    N_CURRENT
    //Sample sizes
    20
    //Growth rates
    0
    //Number of migration matrices : If 0 : No migration between demes
    0
    //Historical event: time, source, sink, migrants, new deme size, new growth rate, new migration matrix
    1 historical events
    T_BOTTLENECK 0 0 1 20 0 0
    //Number of independent chromosome
    1 0
    //Number of contiguous linkage blocks
    1
    //Per Block: Data type, No. of loci, Recombination rate to the right-side locus, plus optional parameters
    FREQ 1 0 1.0e-07

Your `bottlneck.est` file should look like this:

    // Parameter file for the bottleneck model
    // ***************************************
    [PARAMETERS]
    //#isInt? #name #dist. #min #max
    1 T_BOTTLENECK unif 1 10000 output
    0 SIZE_CHANGE unif 0 10 output
    1 N_CURRENT unif 1 10000 output

</details>

<br/>

In this last estimation you probably noticed that the estimates are much
worse than before. This has two reasons:

1.  There are actually many parameter combinations that are compatible
    with the data. You can check that by running the optimization
    multiple times and comparing the final estimates and obtained
    likelihoods. While the estimates will vary, the likelihoods are not
    so dissimilar. That is in part because `fastsimcoal2` (by default)
    only fits the polymorphic part of the SFS. It may does make sense to
    fix one reference parameter (e.g. the current size, as we did
    before), which ensures that only one combination is valid. Or to
    help the optimizer by setting a more narrow initial search range.

2.  In addition, the Brent optimization algorithm implemented in
    `fastsimcoal2` is not the most efficient. You should see this, as
    the values often stop improving after a few iterations while a new
    run may result in a better likelihood. It is thus generally
    recommended to perform many optimization runs and pick the best
    likelihood among them.

## Inferring an SFS from sequence data

You find sequence data (as `bam` files) for 10 diploid individuals (20
alleles) in the repository under `Daniel_Wegmann/day5/data`. This data
is from a population that went through a recent bottleneck. It has a
current effective size of 500 individuals and the goal of this series of
exercises is to make use of the sequence data to learn about the the
timing of the bottleneck and the relative ancestral size.

There are several bioinformatic steps to get from raw sequence data to
an estimated SFS:

1.  We will use `atlas` to estimate sequencing errors and PMD, i.e. to
    recalibrate quality scores.
2.  We will use `atlas` to compute genotype likelihoods using the
    learned errors.
3.  We will use `atlas` to calculate site allele frequency likelihoods
    for the full set of samples.
4.  We will use `winsfs` to estimate the SFS from those samples.

### Installing `atlas`

The easiest way to install `atlas` is via conda. To install packages on
your VM via conda, you first need to create a new environment (I use the
name `demoinf`, but you are free to use any name):

``` bash
conda create -n demoinf python=3.11
```

This will guide you through the process of creating a new environment
(just say yes to everything it wants to install).

Then, activate that new environment and install `atlas`:

``` bash
conda activate demoinf
conda install atlas
```

Check that `atlas`is properly installed by running it without arguments,
in which case you should see a list of all available tasks:

``` bash
atlas
```

**Step 1) estimate errors**: Let us now use `atlas` to learn about PMD
and recalibrate sequencing errors. This is done with task
`estimateErrors` and should be done on each `bam` file individually. To
speed up this tutorial, however, we just do it on the very first `bam`
file and then used the learned error parameters for all the others (this
is OK here since I simulated all `bam` files with the same error model).
Let us work directly in the `Daniel_Wegmann/day5/data` folder (this will
take a while, so read on):

``` bash
atlas estimateErrors --minDeltaLL 0.1 --fasta reference.fasta --bam mysterious_population_ind1.bam
```

Apart from the `bam` file, we also provide the reference with `--fasta`
and use `--minDeltaLL 0.1` to tell `atlas` to use a somewhat less
stringent criterion for optimization to speed up things for this
tutorial. Once finished, this will result in a file
`mysterious_population_ind1_RGInfo.json` that contains the error and PMD
parameters for this sample.

**Step 2) compute genotype likelihoods**: The bam-file format is a
per-read file format, however, for most downstream analysis, and
especially when comparing different samples, we want to do this per site
(position on chromosome). For this, the Genotype Likelihood File Format
(`glf`) is convenient, which stores for each site the likelihoods of the
ten genotypes: AA, AC, AG, AT, CC, CG, CT, GG, GT, TT.

Let us create the `glf` files for each individual using the task `GLF`
and by providing learned error parameters:

``` bash
for bam in *.bam; do
    atlas GLF --bam $bam --RGInfo mysterious_population_ind1_RGInfo.json
done
```

This will result in `glf` files with ending `*.glf.gz` and corresponding
index-files with ending `*.glf.idx`. You can have a look at one of them
using the `atlas` task `printGLF`:

``` bash
atlas printGLF --glf mysterious_population_ind1_GLF.glf.gz | less
```

**Step 3) calculate site allele frequency likelihoods**: We next make
use of the individual genotype likelihoods to calculate the likelihoods
of all possible allele frequencies for all sites. We do so using the
task `saf` of `atlas`. This task takes as input the list of all `glf`
files to use, so let’s first compile that list and then pass it on to
`atlas`:

``` bash
samples=$(ls -1 *.glf.gz | paste -s -d ',' -)
atlas saf --glf $samples --fasta reference.fasta
```

This will result in a file called

**Step 4) infer SFS with `winsfs`**: The last step is to infer the SFS
itself. For this we will make use of `winsfs`
(<https://github.com/malthesr/winsfs>). To install it, you need to first
install a Rust toolchain:

``` bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
```

You will be asked to install several things along - just agree to all of
them ;-)

You can then install `winsfs`:

``` bash
cargo install winsfs-cli
```

Once installed, you can use it to infer the SFS:

``` bash
winsfs saf.saf.idx > obs.sfs
```

## Inferring bottleneck parameters from the inferred SFS

Let us now use the inferred SFS to infer the time of the bottleneck
(T_BOTTLENECK) and the population size change (SIZE_CHANGE) while fixing
the current population size to 500 (copy the relevant files from the
exercises). But we first need to do three changes to the SFS we
inferred:

1.  If you look at the inferred SFS found in the file `obs.sfs` you will
    notice that `winsfs` uses a slightly different format than expected
    by `fastsimcoal2`. Let us thus modify that format.
2.  One issue when inferring an SFS is that we often do not have
    information about which allele is ancestral and which one is
    derived. Since we did not provide that information when inferring
    the `saf` file, it is safest to fold the SFS.
3.  `fastsicoal2` expects the SFS to be in a file with a particular
    name.

So let’s change all these things:

``` bash
(
echo "1 observations"
echo "d0_0  d0_1    d0_2    d0_3    d0_4    d0_5    d0_6    d0_7    d0_8    d0_9    d0_10   d0_11   d0_12   d0_13   d0_14   d0_15   d0_16   d0_17   d0_18   d0_19   d0_20"
tail -n+2 obs.sfs |  awk '{for (i = 1; i <= int((NF+1)/2); i++){ if(i==NF - i + 1){ print $i } else { print $i + $(NF - i + 1) } }}' | tr '\n' '\t'
) > bottleneck_MAFpop0.obs
```

You can now infer the demographic parameter from this folded SFS:

``` bash
./fastsimcoal2 -M -t bottleneck.tpl -e bottleneck.est -m -n100000  -L 25 -c 4 -q
```

Note that we now use the tag `-m` rather than `-d` to tell
`fastsimcoal2` that we use a folded SFS and hence the minor allele
rather than the derived allele.

**Q5)** What are the estimates you get? And are they robust to different
runs?

<details>
<summary>
Solutions
</summary>

The data was actually simulated with the same parameters we used above:
a bottleneck 100 generations ago that reduced the population by a factor
20. However, your probably will estimate an older bottleneck that was a
bit less extreme, although there might be some variation across runs due
to a somewhat flat likelihood.

Why are the estimate smore challangeing? There are two reasons: 1. I
only simulated 2Mb of data, which results in a nosier SFS. 2.
Recalibration may not be perfect (also because there is little data),
further introducing noise.

However, you shoudl still reliably identify the presence of a
bottleneck.
</details>

<br/>

**Q6)**: Also obtain demographic parameters from an SFS estimated
without error recalibration (without providing the `json` file when
calculating `glf` files). Does recalibration matter?
