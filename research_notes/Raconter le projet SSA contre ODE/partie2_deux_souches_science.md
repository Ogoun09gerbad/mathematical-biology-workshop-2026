# Part 2 scientific background: two-strain SIRS with waning immunity, ODE vs Gillespie SSA

Method note: sources were reached by web search and by fetching abstracts or full text for a few pages (Minayev & Ferguson 2009, Fowler 2021, Gog & Grenfell 2002). Several classic papers (Bremermann & Thieme 1989, Castillo-Chavez et al. 1996, Kamo & Sasaki 2002, Andreasen et al. 1997, Mollison 1991) were NOT read in full; the claims below about them come from bibliographic records or later citing papers, as flagged.

Key framing point for the report writer: our model is strain-specific immunity with NO cross-protection (R_i is fully susceptible to strain j; there appears to be no coinfection and no cross-protective term). Most classical "competitive exclusion" theorems are for TOTAL cross-immunity or no recovery (SIS), which is a different structural assumption. The literature below is split accordingly.

## Q1. Competitive exclusion in multi-strain models: strain with largest R0 wins

### Takeaway
The classical result (Bremermann & Thieme 1989) is for models where hosts infected with, or immune to, one strain are protected against the others (total cross-immunity / shared susceptible pool); there the strain with the largest R0 persists and the others die out. Our "no cross infection" variant (strain 2 excludes strain 1) is an instance of this structure, not of the strain-specific-immunity case.

### Cited Findings
- Bremermann & Thieme, "A competitive exclusion principle for pathogen virulence", J. Math. Biol. 27, 179-190 (1989) (record: [Springer](https://link.springer.com/doi/10.1007/BF00276102)). Described in secondary sources as a multi-strain SIR (modified Anderson-May) ODE model in which the strain with the largest reproduction number survives ([search summary of citing works, e.g. Bichara et al. / arXiv 2308.10348](https://arxiv.org/pdf/2308.10348)). Full text not accessed.
- A Purdue-hosted tribute describes the result as holding "without multiple infections and with complete cross-immunity": only the parasite strain with the maximal basic replacement ratio survives — [Purdue MBE special issue](https://www.math.purdue.edu/~fengz/documents/publications/mbe_special_issue.pdf) (reported through a search summary; not read directly).
- Castillo-Chavez, Huang & Li, "Competitive exclusion in gonorrhea models and other sexually transmitted diseases", SIAM J. Appl. Math. 56(2), 494-508 (1996): two-strain two-sex homogeneous SIS model; coexistence not possible except under very special circumstances; includes global stability of the stationary states — [SIAM](https://epubs.siam.org/doi/10.1137/S003613999325419X) (abstract-level, via search summary). Earlier Cornell tech report: [BU-1217-M](https://ecommons.cornell.edu/bitstream/1813/31801/1/BU-1217-M.pdf).
- Ackleh & Allen, "Competitive exclusion in SIS and SIR epidemic models with total cross immunity and density-dependent host mortality", Discrete Contin. Dyn. Syst. B 5 (2005) 175-188, doi 10.3934/dcdsb.2005.5.175 (citation as given in a reference list in [Bichara, Iggidr & Sallet](https://ouci.dntb.gov.ua/works/40EeWY3l); theorems not read). Their earlier paper: "Competitive exclusion and coexistence for pathogens in an epidemic model with variable population size", J. Math. Biol. 47(2):153-168 (2003) (same source). Note: variable-population/density-dependent mortality can alter the outcome, relevant since our model has demography (but with constant birth/death mu; the 2003/2005 results concern density-dependent mortality, so not directly applicable).
- Martcheva's textbook "An Introduction to Mathematical Epidemiology" (Springer 2015) has a chapter "Multi-strain Disease Dynamics" — [table of contents listing](https://tenlong.com.tw/products/9781489978325). Chapter content not read.
- Exclusion principle reported to be extended to vector-borne, patchy, etc. settings, and to fail in some cases (periodic environments, host age structure, superinfection, coinfection, cross-immunity) — [search summary of Martcheva-related papers](https://people.clas.ufl.edu/maia/files/Main-V7.pdf) (not read in detail).
- Review on this line of work cited by that summary: Li, Liu & Martcheva, Math. Biosci. Eng. 2010 (listed at [AIMS](https://www.aimspress.com/aimspress-data/mbe/2010/1/PDF/1551-0018_2010_1_123.pdf); content not verified here).

### Inferences
- In our no-cross-infection variant (each strain only infects S, with I1/I2 competing for one susceptible pool), strain 2 (R0=3.84) excluding strain 1 (R0=2.31) matches the Bremermann-Thieme prediction. Whether the sheet's equations exactly satisfy the theorem's hypotheses (SIRS rather than SIR, waning) should be checked by the student; I did not find a source covering exactly this SIRS case.
- Castillo-Chavez et al. concerns SIS; do not cite it as covering SIRS.

### Gaps
- Original Bremermann-Thieme and Castillo-Chavez full texts were not read; no direct source for the SIRS-with-waning exclusion theorem.
- Andreasen / Lipsitch / Levin exclusion papers named in the brief were not located individually (see Q2 for Andreasen, Lin & Levin 1997 which is about partial cross-immunity).

## Q2. Coexistence with strain-specific (no cross-protective) immunity; role of waning immunity

### Takeaway
Literature on strain-specific immunity mostly treats cross-immunity as a parameter sigma between 0 and 1; coexistence is favoured when cross-protection is weak and when infections are long relative to host lifetime (Gog & Grenfell). Cross-protection and waning can also produce oscillations. I did not find a paper treating exactly our model (two strains, R_i susceptible to the other strain, SIRS with waning, demography); the coexistence we find is consistent with the general idea but the citation is by analogy.

### Cited Findings
- Gog & Grenfell, "Dynamics and selection of many-strain pathogens", PNAS 99(26):17209-17214 (2002), doi 10.1073/pnas.252512799 — [PMC139294](https://pmc.ncbi.nlm.nih.gov/articles/PMC139294). Model tracks each host's immune status with one variable per strain, cross-immunity matrix sigma_ij ("polarized" immunity: partial cross-immunity either fully protects or does nothing); strains cluster; long infectious period relative to host lifespan lets clusters coexist, short infectious period leaves a single dominant cluster. Source also says long-lived immunity tends to limit circulation to one cluster at a time. (As reported from a fetch of the abstract/page, paraphrased.) Note: their infection is SIR-type with immunity, so sigma=0 corresponds to our independent-strain case, but this was not stated explicitly in what I read.
- Andreasen, Lin & Levin, "The dynamics of cocirculating influenza strains conferring partial cross-immunity", J. Math. Biol. 35, 825-842 (1997), doi 10.1007/s002850050079 (citation confirmed via [reference list at MMNP](https://www.mmnp-journal.org/articles/mmnp/ref/2012/03/mmnp201273p49/mmnp201273p49.html)). Findings not read.  Follow-up: Lin, Andreasen & Levin, "Dynamics of influenza A drift: the linear three-strain model", Math. Biosci. 162, 33-51 (1999) (same reference list).
- Kamo & Sasaki, "The effect of cross-immunity and seasonal forcing in a multi-strain epidemic model", Physica D 165, 228-241 (2002): two-strain SIR with seasonal forcing and a cross-immunity parameter; period-doubling bifurcations with increasing forcing — bibliographic details and abstract-level summary via a library catalogue (search summary; full text not accessed). Not a waning-immunity result.
- Gupta, Maiden, Feavers, Nee, May & Anderson, "The maintenance of strain structure in populations of recombining infectious agents", Nature Medicine 2(4):437-442 (1996): strain structure with non-overlapping antigenic repertoires maintained via immune selection (citation confirmed via PubMed-style listing in a search result; description from the later [Buckee et al. PNAS 2011 summary](https://pubmed.ncbi.nlm.nih.gov/21876129/)). Concerns strain structure of recombining pathogens, not directly our model.
- Nuno, Feng, Martcheva & Castillo-Chavez, "Dynamics of two-strain influenza with isolation and partial cross-immunity", SIAM J. Appl. Math. 65(3), 964-982 (2005): cross-immunity and isolation can produce sustained oscillations — [record](https://redi.cedia.edu.ec/document/159283). (Brief mentions the authors as "Nuno et al."; the actual author list is Nuno, Feng, Martcheva, Castillo-Chavez, per the record.)
- Bianco, Shaw & Schwartz, "Epidemics with multistrain interactions: the interplay between cross immunity and antibody-dependent enhancement", Chaos (2009): weak cross immunity stabilizing, strong cross immunity numerically produces oscillations and chaos — [arXiv:0812.3840](https://arxiv.org/pdf/0812.3840) (via search summary).
- A 2021 arXiv paper summarizes Gog & Grenfell as showing that coexistence of multiple cross-immune strains requires sufficiently long immunity duration — [arXiv:2111.07762](https://arxiv.org/pdf/2111.07762v1) (secondary).
- Waning immunity and oscillations: Hethcote, Stech & van den Driessche, "Periodicity and stability in epidemic models: a survey" (1981), in Differential Equations and Applications in Ecology, Epidemics and Population Problems; the citation appears in a Springer reference list but I could not verify details — [search result](https://link.springer.com/article/10.1007/s10255-003-0093-3). Delayed-immunity SIRS models show Hopf bifurcation and periodic solutions in some cases — [van den Driessche survey (UVic)](https://dspace.library.uvic.ca/handle/1828/1767). Our system's damped oscillations around a stable endemic point are the generic ODE behaviour: "damped oscillations to a stable focus" for SIR with demography, and endemic equilibria approached via damped oscillations — [DIMACS slides](https://archive.dimacs.rutgers.edu/archive/Workshops/Analogies/protect/presentations/dimacs-alun.pdf) (secondary).

### Inferences
- With no cross-protection and no coinfection restriction, each strain sees a susceptible pool S + R_j; the strains interact only through the shared host pool/constraints, so coexistence is the expected default, in contrast with the total-cross-immunity case of Q1. I did not find a paper stating this exact claim; present it as an inference, not a cited result.
- The difference between "no cross infection" (equations as printed: exclusion) and "cross infection" (R_i can become infected by strain j: coexistence) is therefore the difference between total cross-protection and no cross-protection.

### Gaps
- No source found for a theorem on SIRS two-strain with zero cross-protection. Kamo & Sasaki, Andreasen et al. and Castillo-Chavez-style papers' specific sigma=0 results not read.
- Whether the Nuno et al. "partial" cross-immunity model reduces to ours at the zero limit was not checked.

## Q3. Invasion analysis / mutual invasibility as the coexistence criterion

### Takeaway
Mutual invasibility (each strain can invade the other's endemic equilibrium) is the standard route to a coexistence equilibrium in two-strain models; invasion number > 1 for both means coexistence in the models where it applies.

### Cited Findings
- In a two-strain spatial SIS model, invasion numbers measure whether each strain can invade the single-strain equilibrium of the other; if both exceed one a coexistence equilibrium exists; an invasion number above one makes the other strain's equilibrium unstable — summarised from a search of Martcheva-linked papers, [Paper-SIS_rev2 (UF)](https://people.clas.ufl.edu/maia/files/Paper-SIS_rev2.pdf) and [Main-V7 (UF)](https://people.clas.ufl.edu/maia/files/Main-V7.pdf). Not read in detail; this is for SIS/spatial, not our SIRS.
- Chapter on multi-strain disease dynamics in Martcheva's textbook (see Q1) is the likely textbook source for invasion reproduction numbers; I did not verify content.
- Diekmann & Heesterbeek, "Mathematical Epidemiology of Infectious Diseases" (2000) and Keeling & Rohani (2008) are plausible textbook sources for invasion R0, but I did not retrieve them; no verified page-level citation.

### Inferences
- Our invasion numbers R0_j (1 - I_i*) = 3.83 and 2.27 > 1 are consistent with the formula R_inv_j = R0_j * (fraction susceptible to j at strain i's equilibrium). Note: in our model the fraction susceptible to j when i is at equilibrium is S + R_i (approx 1 - I_i* - I_j contributions... ), so 1 - I_i* is a consistent approximation only if the other classes are as described; the student should state the exact expression used.

### Gaps
- No verified peer-reviewed source for invasion analysis specific to SIRS with waning and zero cross-protection; textbook page references not obtained.

## Q4. Stochastic extinction of strains in small populations

### Takeaway
Demographic stochasticity generally reduces strain persistence and diversity relative to deterministic models, since deterministic strains can never go fully extinct; extinction risk increases as population size falls. This matches our SSA results (loss at N=100-1000, none at N=10000). The "two strains lost together" result is our own observation; I found no literature on it.

### Cited Findings
- Minayev & Ferguson, "Incorporating demographic stochasticity into multi-strain epidemic models: application to influenza A", J. R. Soc. Interface 6(40):989-996 (2009), doi 10.1098/rsif.2008.0467 — [PMC2827439](https://pmc.ncbi.nlm.nih.gov/articles/PMC2827439/). In deterministic models "a strain, once it becomes established, never goes completely extinct"; adding an extinction threshold and stochastic mutation lowers strain diversity, prevalence and strain lifetimes; lifetimes shorten (slowly) as population size decreases, roughly logarithmically. This is a semi-stochastic model with cross-immunity, not full SSA, and concerns influenza A.
- Lopes & Luczak, "Extinction time for the weaker of two competing SIS epidemics", Ann. Appl. Probab. 30(6), 2880-2922 (2020) — [arXiv:1802.04037](https://arxiv.org/pdf/1802.04037). Two strains in a closed population of size N with perfect cross-immunity; asymptotics for the time at which the weaker strain disappears (stronger strain supercritical). Summary taken from the search result abstract; this is the closest "stochastic competitive exclusion" reference found, but it is SIS with total cross-immunity, i.e. the exclusion regime, not coexistence.
- Stochastic multistrain dengue work: Bianco, Shaw & Schwartz (2009) as above is deterministic; an APS 2009 abstract describes an agent-based ADE model and "the role that fluctuations due to system size have on disease extinction paths" — [APS abstract](https://meetings-archive.aps.org/mar/2009/t9/5/) (abstract only; authors and results not confirmed). The search did not give a stochastic extinction paper by "Shaw & Schwartz" that I could cite beyond this.
- General fade-out/critical community size: smaller or isolated populations show larger extinction probabilities and more violent oscillations — [arXiv:2104.03254](https://arxiv.org/pdf/2104.03254); Bartlett-style fade-out after a major epidemic when few susceptibles remain (same). Keeling & Grenfell 1997 Science on critical community size via [Princeton portal](https://collaborate.princeton.edu/en/publications/disease-extinction-and-community-size-modeling-the-persistence-of/) (the entry titled "Disease extinction and community size: modeling the persistence of measles" appears in search results; I did not verify authors/year from it, beyond that Keeling & Grenfell 1997 was named in the summary).
- Another multi-strain stochastic SIS result: all strains are eventually eliminated since the disease-free state is absorbing; coexistence predicted by deterministic models is unlikely in stochastic models even for short times — search summary of a stochastic SIS multi-strain paper; the identifying source for this was [Deterministic and stochastic epidemic models with multiple pathogens (TDL)](https://tdl-ir.tdl.org/items/f986d2ed-f41e-4e1e-93b6-4c2e559ed380) (summary only; not read, treat with caution).

### Inferences
- A small endemic equilibrium (I1*=0.103, I2*=0.211 means 1 and 2 infected individuals at N=10000?) — note the student should state whether I* are fractions or counts; if prevalence is ~1e-2 to 1e-5 of N the expected infected numbers are tiny and stochastic fade-out is expected. This reasoning is mine, not cited.
- "Two strains lost together" cannot be traced to literature; likely an effect of the mechanism that the two strains are coupled through a shared susceptible pool and a common epidemic cycle (my speculation).

### Gaps
- No stochastic-SSA paper found on two-strain SIRS with zero cross-protection and demography.
- Gog & Grenfell and Bianco/Shaw/Schwartz stochastic extinction results were not found in a verifiable form.
- No reference found for the 85% loss figure comparison; that is a project result.

## Q5. ODE trajectory falls below one individual then rebounds: the "atto-fox" problem

### Takeaway
Yes, there is a standard name: the "atto-fox problem". The criticism of continuous models predicting absurdly small populations that then recover is attributed to Mollison (1991) on the rabies model, and the name itself is attributed to Lobry & Sari (2015). Sources disagree on who coined the phrase.

### Cited Findings
- Fowler, "Atto-Foxes and Other Minutiae", Bull. Math. Biol. 83(10):104 (2021), doi 10.1007/s11538-021-00936-x — [PMC8408093](https://pmc.ncbi.nlm.nih.gov/articles/PMC8408093). From the page: the term "atto-fox problem" is attributed to Lobry & Sari (2015); the issue traces to the rabies model of Anderson et al. (1981); Mollison (1991, p. 31) criticised the continuous model, said it cannot allow extinction and advised stochastic models; oscillation minima fall to about 1 atto-fox (1e-18) per km^2. Fowler's proposed fix: add a reservoir (refuge) for individuals immune from destruction.
- Mollison D., "Dependence of epidemic and population velocities on basic parameters", Math. Biosci. 107:255-287 (1991), doi 10.1016/0025-5564(91)90009-8 (reference as listed in Fowler). 
- Lobry & Sari, "Migrations in the Rosenzweig-MacArthur model and the 'atto-fox' problem", ARIMA 20 (2015) — [ARIMA](https://arima.episciences.org/1990/pdf) (search result; French title "Migrations dans le modele de Rosenzweig-MacArthur et le probleme de l'atto-fox"; I did not read the paper). Their approach: two-patch migration to prevent prey densities from becoming too small.
- Conflict: one search summary, relayed from a Strathclyde paper at [strathprints](https://strathprints.strath.ac.uk/41807/1/17025812.pdf), says Mollison (1991) himself "refers to this ... as the atto-fox problem". Fowler credits Lobry & Sari with the name. I could not resolve this; I recommend citing the concept as "the atto-fox problem (Mollison 1991; named in Lobry & Sari 2015 according to Fowler 2021)" with a hedge.
- A relevant quote paraphrase from an arXiv source summarised in search: continuous populations cannot reach zero, so the residual value can restart a second wave later — [arXiv:1111.6460](https://arxiv.org/pdf/1111.6460) (search summary only).

### Inferences
- Our ODE strain 2 dipping below one individual then rebounding in the no-cross-infection variant is a direct example of the atto-fox problem in an epidemic setting; the SSA's 85% loss is the stochastic counterpart.
- Related epidemiological concept: fade-out (Bartlett; critical community size), see Q4.

### Gaps
- Mollison 1991 text not accessed; whether it uses "atto-fox" is unresolved. The term "atto-fox" for epidemic (non-spatial) models specifically has no source found.
