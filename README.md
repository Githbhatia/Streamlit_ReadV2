Read and plot CISMIP formated v2 OR COSMOS formated V2c seismic ground motion records such as posted on CESMD website (Center of Strong Motion Data https://www.strongmotioncenter.org/). The code can also be used to view free field instrument records on the CESMD site and those also available at the HCAI website (https://hcai.ca.gov/construction-finance/facility-detail/ - navigate to a hospital that has instrumented buildings and look under the Instrumented Buildings Tab). Python code reads a zip file containing .v2 or .V2c files that contains one or three channels (Free-Field instruments have 3 channels)

This is a Streamlit version of the ReadV2 application posted as separate application.  Does not have all the functions of the ReadV2 application yet.

3/17/2025 Absolute trigger for autoranging was failing for earthquake records that have very low accelerations - revised for a dyanmic trigger based on 1/10 maximum absolute value.

Try out at https://appreadv2-8tcju9gckv5rnfcrja59nj.streamlit.app/
