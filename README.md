Read and plot CISMIP formated v2 OR COSMOS formated V2c seismic ground motion records such as posted on CESMD website (Center of Strong Motion Data https://www.strongmotioncenter.org/). The code can also be used to view free field instrument records on the CESMD site and those also available at the HCAI website (https://hcai.ca.gov/construction-finance/facility-detail/ - navigate to a hospital that has instrumented buildings and look under the Instrumented Buildings Tab). Python code reads a zip file containing .v2 or .V2c files that contains one or three channels (Free-Field instruments have 3 channels)

This is a Streamlit version of the ReadV2 application posted as separate application.  Does not have all the functions of the ReadV2 application yet.

3/17/2025 Absolute trigger for autoranging was failing for earthquake records that have very low accelerations - revised for a dyanmic trigger based on 1/10 maximum absolute value.

3/29/2025 Added new file name scheme used by USGS in records from the Burma Earthquake.

3/30/2025 Changed the way v2C files read, added error catching if file names dont match standards.

4/4/2025 Changed ADRS spectra to be a vertical arrangement rather than horizontal to avoid titles overlapping.  Also implemented a logspace interval for spectra.

5/26/2025 Added ROTD50 spectra and associated ASI computation (with polar gragh plotted against azimuth angles).

Try out at https://appreadv2-8tcju9gckv5rnfcrja59nj.streamlit.app/
