# GeoSoftGPF_1
<pre>
Software Sharing for the 27 Jan - 18 May 2011 subset of the DISL FOCAL mooring dataset

The dataset can be found here: http://dx.doi.org/10.5281/zenodo.18943. 
Software is released under CC-BY-SA-NC license. 
Software (except stickplot.m) were written in 2011 by Mimi Tzeng (orcid.org/0000-0001-9396-3217), Dauphin Island Sea Lab, Kyeong Park Lab, with funding from the Fisheries Oceanography in Coastal Alabama (FOCAL) program.
The workflow that connects all of these scripts is illustrated in this diagram: http://dx.doi.org/10.5281/zenodo.20656. It will be described in an upcoming paper: Tzeng, M.W., K. Park, and B. Dzonwkowski. (in prep). A Subset of a Time-Series of Hydrographic and Current Data from a Permanent Moored Station Outside Mobile Bay, Alabama (27 Jan to 18 May 2011). AGU Earth and Space Science.

Software includes: 
1. mooring_all.pl - a perl script for preparing Seabird and YSI data files for import into Matlab (full metadata: http://ontosoft.org/portal/#browse/Software-uruml4oqtqp2). Step D in the workflow.
2. FindMaxDepth.m - contains scrap Matlab code to help locate the beginning and ending scan numbers for CTD vertical profiles. Step E1 in the workflow.
3. FindMoorEnds.m - contains scrap Matlab code to help locate the beginning and ending scan numbers for moored Seabird CTDs and thermistors, and YSI sondes. Step E1 in the workflow.
4. moorburst.m - a Matlab function that locates the beginning and ending scan numbers for calculating 20-minute averaged data from 1-minute measured data, in alignment with 20-minute measured YSI and ADCP data. Part of Step E2 in the workflow.
5. MOORprocess_all.m - a Matlab script for processing time-series data from moored CTDs, YSI sondes, and thermistors (full metadata: http://ontosoft.org/portal/#browse/Software-t4ct4jpj1hys). Part of Step E2 in the workflow.
6. stickplot.m - a Matlab function that generates stickplot figures. This function was not written by Mimi Tzeng at DISL. It was downloaded from the Internet in 2011 and its provenance is no longer known. Part of Step H1 and H2 in the workflow.
7. FindADCPendpoints.m - a Matlab script that generates a plot of pitch and roll over time, and stickplots for the twenty vertical bins closest to the surface, to help locate the beginning and ending scan numbers and surfacemost vertical bin for moored ADCPs. Part of Step H1 in the workflow.
8. MoorADCP.m - a Matlab script for processing time-series data from moored ADCPs (full metadata: http://ontosoft.org/portal/#browse/Software-1iv9uqkya98ol). Part of Step H2 in the workflow.

</pre>
