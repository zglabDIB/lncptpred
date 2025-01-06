<h1>LncPTPred</h1>
LncPTPred is a target prediction tool to predict the interaction between lncRNA and Protein. Given the lncRNA sequence and the Protein name, this tool will extract the position specific binding probabilities across 
the given lncRNA loci using a shifting window of specific size. Web-server version corresponding to the LncPTPred tool can be accessed from  <a href='http://bicresources.jcbose.ac.in/zhumur/lncptpred/'>Here</a> or <a href='http://dibresources.jcbose.ac.in/zhumur/lncptpred/'>Here</a>. <br/>

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
<li>Install the <b>git-lfs</b> in order to download the pre-trained models(<a href='https://git-lfs.com/'>git-lfs</a>)</li>
<li>Clone the whole application using following command
<ul><li><b><i>git clone https://github.com/zglabDIB/lncptpred.git</b></i></li></ul>
<li>Create a conda environment named viz. lncptpred using the following command:</li>
<ul><li><b><i>conda create -n lncptpred python==3.8.5</b></i></li></ul>
<li>Activate the conda environment using</li>
<ul><li><b><i>conda activate lncptpred or source activate lncptpred</b></i></li></ul>
<li>Move to that downloaded lncptpred directory using</li>
<ul><li><b><i>cd "downloaded_path"/lncptpred/</b></i></li></ul>
<li>Download the necessary python packages from the requirement.txt file provided here using</li>
<ul><li><b><i>pip  install -r requirement.txt</b></i></li></ul>
<li>Upon failing the installation of <b>git-lfs</b> user can download the pre-taring models using following command (linux):</li>
<ul><li><b><i>wget https://github.com/zglabDIB/lncptpred/blob/master/models/LR_Meta_lncRNA_Prot_Balanced.sav</b></i></li></ul>
<ul><li><b><i>wget https://github.com/zglabDIB/lncptpred/blob/master/models/LR_Meta_lncRNA_Prot_GLOBAL.sav</b></i></li></ul>
<ul><li><b><i>wget https://github.com/zglabDIB/lncptpred/blob/master/models/LR_Meta_lncRNA_Prot_Pos_Bias.sav</b></i></li></ul>
<ul><li><b><i>wget https://github.com/zglabDIB/lncptpred/blob/master/models/StandardScaler_lncRNA_Protein.sav</b></i></li></ul>
<li>In Windows go to the models folder and download individual pre-trained models directly and finally save them into the <b>models</b> folder</li>

</ul>
<h2>Code Execution Procedure:</h2>
<ul>
<li>Copy-Paste the lncRNA Sequence in <b>lncrna_seq.txt</b> file, place it in <b>input</b> folder and execute it by the following code:</li>
<ul>
<li><b><i>python init.py</b></i></li>
<li>Select Window length of lncRNA Sequences with <b>Minimum Length=10, Maximum Length=30</b> and <b>Default Length=20</b>. Given <b>Length<10</b> and <b>>30</b> will be considered as 10 and 30 respectively.</li>
<li>Select Shifting size of moving window corresponding to lncRNA Sequence with <b>Minimum Size=1, Maximum Size=5</b> and <b>Default Size=1</b>. Given <b>Size<1</b> and <b>>5</b> will be considered as 1 and 5 respectively.</li>
<li>Select Strands either <b>'+'</b> or <b>'-'</b>(Provide exact strand information corresponding to the input lncRNA sequence)</li>
<li>Enter number to select protein from following list:<br/>
<b>0=>YBX1<br/></b>
<b>1=>ALKBH5<br/></b>
<b>2=>G3BP1<br/></b>
<b>3=>G3BP2<br/></b>
<b>4=>SF3B1<br/></b>
<b>5=>C17orf85<br/></b>
<b>6=>LIN28B<br/></b>
<b>7=>ZCCHC4<br/></b>
<b>8=>HuR<br/></b>
<b>9=>ZFP36<br/></b>
<b>10=>SRSF2<br/></b>
<b>11=>SRSF1<br/></b>
<b>12=>FMR1<br/></b>
<b>13=>SNRPA1<br/></b>
<b>14=>TAF15<br/></b>
<b>15=>EWSR1<br/></b>
<b>16=>ESR1<br/></b>
<b>17=>CAPRIN1<br/></b>
<b>18=>C22orf28<br/></b>
<b>19=>MSI2<br/></b>
<b>20=>LIN28A<br/></b>
<b>21=>RC3H1<br/></b>
<b>22=>YTHDC1<br/></b>
<b>23=>QKI<br/></b>
<b>24=>ZC3H7B<br/></b>
<b>25=>YTHDF1<br/></b>
<b>26=>YTHDF2<br/></b>
<b>27=>FUS<br/></b>
<b>28=>AUF1<br/></b>
<b>29=>METTL14<br/></b>
<b>30=>PTBP2<br/></b>
<b>31=>METTL3<br/></b>
<b>32=>WTAP<br/></b>
</li>
</ul>

<li>After successful execution it will generate <b>prediction_output.txt</b> and <b>motif_logo.png</b> in <b>output</b> folder.</li>
<li><b>prediction_output.txt</b> gives the start & end location along with their sequences corresponindg to their binding sites where the <b>Final_Interacting_Score</b> score denotes probability of interaction associated with the binding segment.</li>
<li><b>motif_logo.png</b> generates the motif plot associated with the interacting lncRNA sequence segments. Those segments are considered for motif plots which are predicted as positively interacting by at least three models.</li>   
</ul>



