Glossaries
###################

.. glossary::

    Tag
        Artificial sequence at the front of the pitcher. It is a signal for catcher that the target sequence is
        amplified. It was designed considering various factors such as mTm, self-dimer and GC contents.

    Pitcher
        Sequence consists of a tag and a probe. Tag is cleaved to send a signal to the catcher
        when the target sequence is amplified.

    Con
        Abbreviation of "C"leavage z"ON"e. 5' end region of a probe which is cleaved along with a tag.

    Catcher
        Sequence consists of catcher front, variability region and reverse-complemented tag. It contains a fluorescent
        dye and a quencher, so catcher could not emit light for itself. When a target is amplified, tag is cleaved from
        a pitcher, and catcher catches it. Then, DNA polymerase elongates the tag sequence with the catcher as a
        template. In result, catcher forms double-strand structure and it causes that the dye is separated with the
        quencher resulting in shine.

    Catcher front
        Artificial sequence at the 5' end of the catcher. It was designed considering various factors such as mTm,
        self-dimer and GC contents.

    Variability region
        Middle region of the catcher. It is designed not to bind with corresponding region of the pitcher.

    TOM
        TOM is an abbreviation of TOCE oligo mix (pitcher and catcher).

    Group
        In MultiTOM, group means a suite of probes whose Tm, dye, Con sequence,
        tag length and target type (PTOD, PMOD) are same. Probes in each group share
        same tag and catcher.

    mTm
        Modeling the melting temperature. Unlike a general Tm calculation method,
        mTm method use experimental values for Tm calculation, so more accurate.

    Position score
        The score to measure how much tag and catcher interact with upstreams and probes,
        respectively. Each position is weighted by its influence on interaction. Because the
        interaction is not desirable, The lower score, the better.

    Hairpin dG
        The simulated free energy of a hairpin structure. The larger value, the more unstable.
        Because a hairpin structure could hinder PCR, large value is preferred.

    Clique
        In the graph theory (Mathematics), a clique is a subset of vertices of an undirected
        graph such that every two distinct vertices in the clique are adjacent. In MultiTOM,
        It is adopted to pitcher and catcher combination, so Clique means that every pair of
        oligonucleotides do not form dimer with each other.
