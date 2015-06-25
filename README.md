# plann
Plann: a command-line application for annotating plastome sequences

Usage:
    plann.pl -reference gbfile.gb -fasta denovoplastome.fasta -out outfile
    [-organism "Genus species"] [-sample samplename]

Options:
      -reference:       a well-annotated plastome reference sequence, in genbank file format
      -fastafile:       the plastome sequence to be annotated, in fasta format
      -outfile:         the output name (default is "output")
      -organism:        [optional: scientific name for Genbank annotation]
      -sample:          the name of the plastome sample (default is the name in the fasta file)

