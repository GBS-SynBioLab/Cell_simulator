# Cell_simulator

Purpose
--
Integrated development environment for the design of synthetic biological systems, from DNA sequence to systems level


Introduction
--
Mathematical modeling module includes a software tool for host aware design of synthetic gene circuits. Host context is implemented based on the coarse-grained mechanistic whole cell model from Weisse, 2015.


Documentation
--

```function diluted_species = degradation(species)```
* variable ```species``` denotes the diluted molecular species (M) diluted

```function p = translation(ribocomplex,codons)```
* variable ```ribocomplex``` denotes the ribosome/mRNA complex species (M)
* variable ```codons``` denotes size of mRNA in bases triplets

```function mRNA = transcription(tx_rate)```
* variable ```tx_rate``` denotes the promoter transcription rate (M/sec)

```function r_or_mRNA = ribosome_unbind(complex)```
* variable ```complex```  denotes the  dissociating ribocomplex (M) 

```function ribocomplex = ribosome_bind(species)```
* variable ```species``` denotes the  molecular species (M) bound

```function dil_species = dilution(species)```
* variable ```species``` denotes the diluted molecular species (M/sec) diluted
