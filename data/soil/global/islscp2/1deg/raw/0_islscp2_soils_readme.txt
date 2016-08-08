AH Comment:

(i) WEB SIDE:
http://islscp2.sesda.com/ISLSCP2_1/html_pages/groups/hyd/islscp2_soils_1deg.html


############################################################################
#     Global Gridded Surfaces of Selected Soil Characteristics for the     #
#                  ISLSCP Initiative II Data Collection.                   #
############################################################################

The Global Soil Data Task was an international collaborative project with the 
objective of making accurate and appropriate data relating to soil properties 
accessible to the global change research community. The collaborators are holders 
of major international pedosphere data sets, such as the United States 
Department of Agriculture (USDA), the Food and Agriculture Organization (FAO) 
of the United Nations, and the International Soil Reference and Information 
Centre (ISRIC), as well as national soils institutes, individual soil scientists, 
and users of soil data. The task was coordinated by the now defunct Data and 
Information System framework activity of the International Geosphere-Biosphere 
Programme (IGBP-DIS).

The Global Soil Data Task assembled a reliable and accessible data set on 
pedosphere properties on a global scale. The data are sufficiently detailed to 
support rigorous analysis, and accessible to and understandable by soil 
scientists and non-soil scientists alike. The immediate goal was to supply 
pedosphere information to global change researchers, but the unification and 
distribution of the existing soils data bases is of great benefit to all 
researchers, not only those directly involved in global change studies. The 
task was enabled by the active participation of the principal international 
custodians of pedosphere data.

A CD-ROM disk containing extensive global pedon data, and the so-called 
"SoilData" software to create global gridded layers of selected soil 
parameters from the FAO Digital Soil Map of the World (DSMW) was one 
outcome of this IGBP-DIS activity which is currently available from the 
ORNL DAAC at http://daac.ornl.gov. The ISLSCP II staff has used the IGBP-DIS 
data and software to generate two-dimensional gridded maps of 18 selected 
soil parameters, including soil texture, at a 1 by 1 degree spatial 
resolution and for two soil depths. All data layers have been adjusted to 
match the ISLSCP II land/water mask.

Any problems encountered with the data set should be reported to Eric Brown de 
Colstoun at ericbdc@ltpmail.gsfc.nasa.gov.

###############################################################################
File Name Description:
----------------------
This data set provides gridded data for selected soil parameters derived from 
data and methods developed by the Global Soil Data Task, coordinated by the 
now defunct Data and Information System (DIS) of the International 
Geosphere-Biosphere Programme (IGBP), or IGBP-DIS. The ISLSCP II data sets 
have been produced by the ISLSCP II staff from an IGBP-DIS soil data CD-ROM 
distributed by the Oak Ridge National Laboratory Distributed Active Archive 
Center (ORNL DAAC, http://daac.ornl.gov/). For the ISLSCP II data collection, 
gridded global maps of 18 selected soil parameters, including soil texture, 
are provided on a 1 by 1 degree Earth grid and for two soil depths (0-30cm 
and 0-150cm). Note that additional maps at finer spatial resolutions are 
being produced and will be available in a future revision of this data set. 
The soils data files for this data collection are named using the following 
naming convention:

      soil_parameterdepth_1d.asc

where:
parameter  is the particular soil parameter (see Section 3.2 for a list of 
           parameters).
depth 		   is the soil depth in cm for the estimated soil parameter (i.e. 
           0-30 or 0-150).
1d 	       identifies the spatial resolution of the data as 1 degree in 
           both latitude and longitude.
.asc 	    	identifies the format of the data as ASCII, or text format.

NOTE: The files for soil thermal capacity are named "soil_therm_capWC_depth_1d.asc", 
where WC is the percent water content of the soil: 0, 10, 50 or 100 percent.

##############################################################################
ASCII File Format:
------------------
All of the files in the ISLSCP Initiative II data collection are in the ASCII, 
or text format. The file format consists of numerical fields of varying length, 
which are delimited by a single space and arranged in columns and rows. The 
files each contain 360 columns by 180 rows. All values in these files are 
written as real numbers. Cells over water or permanent ice are assigned the 
value -99 or -77, respectively, on all data layers except soil texture, where 
they are assigned the values 0 and 13.

All files are gridded to a common equal-angle lat./long. grid, where the 
coordinates of the upper left corner of the files are located at 180oW, 90oN and 
the lower right corner coordinates are located at 180oE, 90oS. Data in the files 
are ordered from North to South and from West to East beginning at 180o West and 
90o North. 

The data files are PKZip compressed. On UNIX, they can be decompressed using the 
"unzip -a" command. On Windows, they can be decompressed using WinZip or other 
PKZip software. On the Macintosh, they can be decompressed using Stuffit 
Expander.

##############################################################################
