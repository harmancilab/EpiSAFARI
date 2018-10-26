<html>
<font face="arial">
<!------------------- <title>EpiSAFARI</title> ---->
<div style="text-align:center"><img src="episafari_logo.png" alt="Could not load logo." width="1000" align="center"></div>
<h1>EpiSAFARI</h1>
EpiSAFARI is a command line tool for detection of signal features from the functional genomics signal profiles. <br><br>
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
- Signal Profile in bedGraph format, <br>
- Mapped reads in SAM/bowtie/eland/... format. <br>
</div><br>
and outputs:<br>
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
- The valleys and local minima and maxima in the signal profile. <br>
</div>
<br>

The output is formatted as an extended BED file with multiple columns. Please refer below for the specification of output file format.

<h2>Download and Installation</h2>
You can download EpiSAFARI C code <a href="https://github.com/gersteinlab/EpiSAFARI/archive/master.zip">here</a>. 

You need to have GSL libraries installed for building EpiSAFARI. If they are not installed, type:
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<i><font face="courier">
sudo yum -y install gsl-devel
sudo yum -y install gsl
<br>
</font></i>
</div><br>

This should install the necessary GSL libraries for building EpiSAFARI correctly.

<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<i><font face="courier">
unzip master.zip<br>
cd EpiSAFARI<br>
make clean<br>
make
</font></i>
</div><br>
to build EpiSAFARI. The executable is located under directory <font face="courier">bin/</font>. 

It may be useful to install <a href="http://samtools.sourceforge.net/">samtools</a> for processing BAM files.

To get help on which options are available, use:
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<i><font face="courier">
EpiSAFARI -help
</font></i>
</div>

<h2>Usage</h2>
EpiSAFARI run starts with setting up the input files. (Note that we use samtools for converting BAM file to SAM files.). EpiSAFARI can take bedGraph files and mapped reads
directly as input. It is necessary to divide the data into chromosomes.

<h2>Building input with bedGraph and bigWig files</h2>
We show an example from ENCODE project below:
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<i><font face="courier">
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneGm12878H3k04me3StdSigV2.bigWig <br>
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph <br>
chmod 755 bigWigToBedGraph <br>
./bigWigToBedGraph wgEncodeBroadHistoneGm12878H3k04me3StdSigV2.bigWig wgEncodeBroadHistoneGm12878H3k04me3StdSigV2.bigWig.bgr <br>
mkdir bedGraphs <br>
EpiSAFARI -separate_bedGraph_2_chromosomes wgEncodeBroadHistoneGm12878H3k04me3StdSigV2.bigWig.bgr bedGraphs <br>
</font></i>
</div><br>
If there are multiple replicates to be pooled, they can be done at once or separately. If done separately, EpiSAFARI pools the bedGraphs automatically.<br>

<h2>Building input with mapped read files</h2>
EpiSAFARI can also process mapped read files, for example in SAM format. We show again an example from the ENCODE Project: <br>
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<i><font face="courier">
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneGm12878H3k04me1StdAlnRep1V2.bam <br>
mkdir processed_reads <br>
samtools view wgEncodeBroadHistoneGm12878H3k04me1StdAlnRep1V2.bam | EpiSAFARI -preprocess_reads SAM stdin processed_reads <br>
</font></i>
</div><br>
If there are multiple replicates of reads to be pooled, they can be done at once or separately. If done separately, EpiSAFARI pools the reads automatically.<br><br>

We strongly suggest removing duplicates from the reads. This decreases feature identification time quite much:<br>
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<i><font face="courier">
mkdir processed_reads/sorted processed_reads/dedup <br>
EpiSAFARI -sort_reads processed_reads processed_reads/sorted <br>
EpiSAFARI -remove_duplicates processed_reads/sorted 2 processed_reads/dedup <br>
</font></i>
</div><br>

<h2>Spline Fitting to the Data</h2>
After input files are setup, we perform spline fitting to the read coverage signals: <br>
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<i><font face="courier">
n_spline_coeffs=20<br>
spline_order=20<br>
max_max_err=5<br>
max_avg_err=1<br>
l_win=1000<br>
<br>
# bedGraph files: <br>
EpiSAFARI -bspline_encode bedGraphs ${n_spline_coeffs} ${spline_order} ${max_max_err} ${max_avg_err} ${l_win}<br>
# mapped read files: <br>
EpiSAFARI -bspline_encode processed_reads/dedup ${n_spline_coeffs} ${spline_order} ${max_max_err} ${max_avg_err} ${l_win}<br>
</font></i>
</div><br>

We finally do feature identification. We first download the multi-mappability signal then identify the features:<br>

<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<i><font face="courier">
wget -c http://archive.gersteinlab.org/proj/MUSIC/multimap_profiles/hg19/hg19_36bp.tar.bz2<br>
tar -xvjf hg19_36bp.tar.bz2<br>
# bedGraph files: <br>
EpiSAFARI -get_significant_extrema bedGraphs 1000 5 1.2 125 hg19_36bp 1.2<br>
# mapped read files: <br>
EpiSAFARI -get_significant_extrema processed_reads/dedup 1000 5 1.2 125 hg19_36bp 1.2<br>
</font></i>
</div><br>

<h2>Output format</h2>
EpiSAFARI outputs the identified features to files named "significant_valleys_[chrosome id].bed" for each chromosome.<br><br>

This is an extended bed file with following columns:<br>
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<i><font face="courier">
<ol>
<li>[chromosome]: Chromosome ID </li>
<li>[Left maxima position]: Position the left maximum of the valley	</li>
<li>[Right maxima position]: Position the right maximum of the valley </li>
<li>[Minima position]: Position of the minimum of the valley </li>
<li>[Left maxima signal]: Signal at the left maxima	</li>
<li>[Right maxima signal]: Signal at the right maxima </li>
<li>[Minima signal]: Signal at the minima </li>
<li>[Average multi-mappability signal]: Average multi-mappability signal on the valley </li>
<li>[Maximum multi-mappability signal]: MAximum multi-mappability signal on the valley </li>
<li>[Left hill fraction]: Fraction of the nucleotides on the left hill</li>
<li>[Right hill fraction]: Fraction of the nucleotides on the right hill</li>
</font></i>
</div><br>

</html>