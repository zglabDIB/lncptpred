<h1>LncPTPred</h1>
LncPTPred is a target prediction tool to predict the interaction between lncRNA and Protein. Given the lncRNA sequence and the Protein name, this tool will extract the position specific binding probabilities across 
the given lncRNA loci using a shifting window of specific size. Web-server version corresponding to the LncPTPred tool can be accessed from  <a href='http://bicresources.jcbose.ac.in/zhumur/lncptpred/'>Here</a>. <br/>

<h2>Pre-requisite:</h2>
<ul>
<li>The query is provided in terms of lncRNA sequence specified in  <b>lncrna_seq.txt</b> file. The example of preformatted input sequence file is provided in <b>example</b> folder.</li>
<li>Sequence characters should be none other than <b>A, C, G, T</b> or <b>U</b>.</li>
<li>Input sequence should be provided in <b>single-line format</b>. Result from multi-line sequence may provide <b>erroneous result</b>.</li>
</ul>
<br/>
<b>Software:</b> For executing the codes, the required list of Softwares are given below.<br/>
<h4>Python programming Environment setup:</h4>
<ul>
<li>Download anaconda distribution from their official website (<a href='https://www.anaconda.com/download'>Here</a>)</li>
<li>Install the <b>git</b> application in order to clone the tool from github(<a href='https://git-scm.com/downloads'>git</a>)</li>
<li>Clone the whole application using following command
<pre>git clone https://github.com/zglabDIB/lncptpred.git</pre>
<li>Create a conda environment named viz. lncptpred using the following command in order to install both <b>python 3.9</b> and the <b>RNAfold</b> from <b>Viennarna</b> :</li>
<pre>conda create -n lncptpred -c bioconda -c conda-forge python=3.9 viennarna -y</pre>
<li>Activate the conda environment using</li>
<pre>conda activate lncptpred or source activate lncptpred</pre>
<li>Move to that downloaded lncptpred directory using</li>
<pre>cd "downloaded_path"/lncptpred/</pre>
<li>Download the necessary python packages from the requirement.txt file provided here using</li>
<pre>pip  install -r requirement.txt</pre>
<li>Unzip the model cloned in the <b>models</b> folder using following command:</li>
<pre>unzip models.zip</pre>
</ul>
<h2>Code Execution Procedure:</h2>
<ul>
<li>Copy-Paste the lncRNA Sequence in <b>lncrna_seq.txt</b> file, place it in <b>input</b> folder and execute it by the following code:</li>
<ul>
<li><b><i>python init.py</b></i></li>
<li>Select Window length of lncRNA Sequences with <b>Minimum Length=10, Maximum Length=40</b> and <b>Default Length=20</b>. Given <b>Length<10</b> and <b>>40</b> will be considered as 10 and 40 respectively.</li>
<li>Select Shifting size of moving window corresponding to lncRNA Sequence with <b>Minimum Size=1, Maximum Size=5</b> and <b>Default Size=1</b>. Given <b>Size<1</b> and <b>>5</b> will be considered as 1 and 5 respectively.</li>
<li>Select Strands either <b>'+'</b> or <b>'-'</b>(Provide exact strand information corresponding to the input lncRNA sequence)</li>
<li>Enter number to select protein from the list of 88 proteins can be found in <b>Final_Protein.txt</b> file located within <b>dataset</b> folder<br/>
</li>
</ul>

<li>After successful execution it will generate <b>prediction_output.txt</b> and <b>motif_logo.png</b> in <b>output</b> folder.</li>
<li><b>prediction_output.txt</b> gives the start & end location along with their sequences corresponindg to their binding sites where the <b>Final_Interacting_Score</b> score denotes probability of interaction associated with the binding segment.</li>
<li><b>motif_logo.png</b> generates the motif plot associated with the interacting lncRNA sequence segments. Those segments are considered for motif plots which are predicted as positively interacting by at least three models.</li>   
</ul>



