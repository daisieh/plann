# plann
Plann: a command-line application for annotating plastome sequences

Publication: Huang, Daisie I and Quentin C B Cronk (2015). Plann: A Command-Line Application for Annotating Plastome Sequences. Applications in Plant Sciences 3(8):1500026. doi:10.3732/apps.1500026

Code: Daisie Huang. (2015). plann: APPS publication version. Zenodo. doi:10.5281/zenodo.27369

Usage:
    plann.pl -reference gbfile.gb -fasta denovoplastome.fasta -out outfile
    [-organism "Genus species"] [-sample samplename]

Options:
      -reference:       a well-annotated plastome reference sequence, in genbank file format
      -fastafile:       the plastome sequence to be annotated, in fasta format
      -outfile:         the output name (default is "output")
      -organism:        [optional: scientific name for Genbank annotation]
      -sample:          the name of the plastome sample (default is the name in the fasta file)

# Updates
- Aug 28, 2015: I have recently discovered that the `/r` option for regular expressions was only introduced in Perl v5.14, so versions earlier than that will fail to compile Plann.pl properly. If this is the case, you will get this error:

```
Bareword found where operator expected at plann-master/lib/Subfunctions.pm line 200, near "tr/[AGCT]//dr"
Bareword found where operator expected at plann-master/lib/Subfunctions.pm line 222, near "tr/[AGCT]//dr"
syntax error at plann-master/lib/Subfunctions.pm line 200, near "tr/[AGCT]//dr"
syntax error at plann-master/lib/Subfunctions.pm line 222, near "tr/[AGCT]//dr"
```

You might have to install a new version of Perl; I recommend using Homebrew (`brew install perl`) if you have it. Make sure that, after installing a new Perl, that you're executing the correct version.
