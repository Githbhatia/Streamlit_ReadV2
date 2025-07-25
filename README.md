Read and plot CISMIP formated v2 OR COSMOS formated V2c seismic ground motion records such as posted on CESMD website (Center of Strong Motion Data https://www.strongmotioncenter.org/). The code can also be used to view free field instrument records on the CESMD site and those also available at the HCAI website (https://hcai.ca.gov/construction-finance/facility-detail/ - navigate to a hospital that has instrumented buildings and look under the Instrumented Buildings Tab). Python code reads a zip file containing .v2 or .V2c files that contains one or three channels (Free-Field instruments have 3 channels)

This is a Streamlit version of the ReadV2 application posted as separate application.  Does not have all the functions of the ReadV2 application yet.

3/17/2025 Absolute trigger for autoranging was failing for earthquake records that have very low accelerations - revised for a dyanmic trigger based on 1/10 maximum absolute value.

3/29/2025 Added new file name scheme used by USGS in records from the Burma Earthquake.

3/30/2025 Changed the way v2C files read, added error catching if file names dont match standards.

4/4/2025 Changed ADRS spectra to be a vertical arrangement rather than horizontal to avoid titles overlapping.  Also implemented a logspace interval for spectra.

5/26/2025 Added ROTD50 spectra and associated ASI computation (with polar gragh plotted against azimuth angles).

5/30/2025 Added Arais Intensity Computation.   Added SI and DSI computation.

6/1/2025 Animation of displacements in 3D added.

6/6/2025 Somewhat optimized Animation - put at the end of application to prevent rerun of animation with other options (still looking for a better option).

6/13/2025 Improved map - included marker for epicenter and location of station.  Added sliders for elevation and azimuth angles for animation.

7/1/2025 Revised as v2c files dont seem to have consistent names.

7/9/2025 Improved compatibility with older COSMOS formats - can now read CESMD files from the Northridge and Loma Prieta Earthquakes.

7/13/2025 Removed requirement that the zip archive contain only one free-field record.  Can now download multiple records at a time from CESMD website and pick one to display.

7/18/2025 Can now read PEER formated files form https://ngawest2.berkeley.edu/ - had to revise many parts to accomodate this new format.

7/19/2025 Added option to write RotD spectra to a file.

Try out at https://appreadv2-8tcju9gckv5rnfcrja59nj.streamlit.app/
