# Nepre-Potential(Radius)
Scoring Function based on Neighbourhood Preference Statistics,include four kinds:
* a function consider Amino acid neighbor preference's position and type,
* a function just consider Amino acid neighbor preference's type combination(set 400 as a mark when calculate energy),
* Energy function for sigle structure(means when create energy matrix consider every amino acid in the pDB structure)
* Energy function for complex structure(means when create energy matrix just consider the neighbor between chain and chain)

Usage
----------
Some explanation about the file contained :
* AminoAcid.py(Class for establish amino acid).
* Energymatrixs(folder contain four kind Energy matrixs).
* EnergyForLoop.py is a program can used to calculate the loop's energy.
* EnergyForStandard.py is a program used to calculate the standard structure's energy.
* Nepre_AA400 return the sigle structure' energy include energy and energy400.
* Nepre_chain400.py return the complex structure'energy include energy and energt400.
* Procecode contain some program about how to create energy matrix and process date.
* A statistic of radius of different kind of amino acid have been done by us. We use gaussian distribute function
to fit them and use the gaussian mean data as the default vlaue. See details in mean_radius.txt.

For single protein potential energy calculate, go to linux shell and type:
<pre><code>
Nepre@liulab:~$ python Nepre.py ./example.pdb
</code></pre>
