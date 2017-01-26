# plann
Plann: a command-line application for annotating plastome sequences

Publication: Huang, Daisie I and Quentin C B Cronk (2015). Plann: A Command-Line Application for Annotating Plastome Sequences. Applications in Plant Sciences 3(8):1500026. doi:10.3732/apps.1500026

Code: Daisie Huang. (2015). plann: APPS publication version. Zenodo. doi:10.5281/zenodo.27369

## Usage:
```
    plann.pl -reference gbfile.gb -fasta denovoplastome.fasta -out outfile
    [-organism "Genus species"] [-sample samplename]
```

## Options:
```
      -reference:       a well-annotated plastome reference sequence, in genbank file format
      -fastafile:       the plastome sequence to be annotated, in fasta format
      -outfile:         the output name (default is "output")
      -organism:        [optional: scientific name for Genbank annotation]
      -sample:          the name of the plastome sample (default is the name in the fasta file)
```
      
## Before you run:
If you don't already have it you will have to install a module from cpan (Comprehensive Perl Archive Network)
```Shell
$ cpan install XML::Simple
```	

## After you run:
Check the `.results.txt` file to see any missing genes that you might want to check prior to submission. If you fix these in your sequence, you can run plann again to make an updated annotation.
	
The files produced by plann can be assembled to make a complete .sqn archive to submit to GenBank, or to further edit with [Sequin](https://www.ncbi.nlm.nih.gov/Sequin/).

1. copy the `plann.sbt` template from the `test/` directory to the directory containing your output from plann
2. use a text editor to edit the `plann.sbt` file to replace the author's information with your information (or do the edits later in Sequin after you assemble the complete .sqn file, but don't forget to do this! Otherwise the credit for the sequence goes to Daisie Huang).
3. While in the directory with all the files use `tbl2asn` ([available from NCBI](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/)) to assemble the complete .sqn file:

	```Shell
	$ tbl2asn -Vbv -p . -t plann.sbt
	```
	
  This can be directly submitted, but it is better to check and edit with Sequin prior to submission.

# Updates
- Mar 25, 2016: 
	Added line breaks in `analysis.sh` to help certain versions of bash.
	Blast on Unix does not seem to like Windows line endings, which seems weird because the Genbank files downloaded from NCBI seem to come with them. Plann now rewrites temp versions of the genbank files with uniformly Unix line endings before running.
	
- Aug 28, 2015: 
	I have recently discovered that the `/r` option for regular expressions was only introduced in Perl v5.14, so versions earlier than that will fail to compile plann.pl properly. If this is the case, you will get this error:

```
Bareword found where operator expected at plann-master/lib/Subfunctions.pm line 200, near "tr/[AGCT]//dr"
Bareword found where operator expected at plann-master/lib/Subfunctions.pm line 222, near "tr/[AGCT]//dr"
syntax error at plann-master/lib/Subfunctions.pm line 200, near "tr/[AGCT]//dr"
syntax error at plann-master/lib/Subfunctions.pm line 222, near "tr/[AGCT]//dr"
```

You might have to install a new version of Perl; I recommend using Homebrew (`brew install perl`) if you have it. Make sure that, after installing a new Perl, that you're executing the correct version.
