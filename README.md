# Plots repository
- This is a special Git repository which is specifically made to create plots for variable data analyses. These analyses are in a separate repo, can be found in *mvranovsky/star-upcDst*. To clone this repository, use the command below:

<pre><code> git clone git@github.com:mvranovsky/plots.git </pre></code>

- Switch to a different branch:

<pre><code> git checkout cleanCodes </pre></code>

- To compile the program, change directory to *build/* and use the commands below:

<pre><code> cd build </code></pre>
<pre><code> cmake ../ </code></pre>
<pre><code> make </code></pre>

- To create plots for specific data, use python script *RunPlot.py*. Before using it, one has to specify path *inputDirectory* to the directory with output from running analysis. 

<pre><code> ./RunPlot.py analysisTag </code></pre>

- Sometimes, the running takes quite a long time. That is why script Run.sh is there. By running:

<pre><code> nohup ./Run.sh analysisTag & </code></pre> 

- the running will not stop when the user leaves the node and the output and error are directed to *output.out* and *output.err*. To differ which analysis should run, the input tag should involve the tag of the analysis. 
- Once the running stops, the output is saved with the data from the analysis as *AnalysisOutput.root*.
- In directory *localDirectory/*, there are several scripts for fitting and downloading the output from RCF to a local computer. The paths have to be modified though.