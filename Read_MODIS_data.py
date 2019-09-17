#
# **************************************************************************
# Read MODIS data by Jose's modis_data tool
# **************************************************************************
#
import modis_data
import gdal, osr
import numpy as np
import struct
import pdb

# Hainich
#tile = "h18v03"
#loc = '/space/ucfamc3/hainich/MODIS/MOD09GA/'
#map_x = 730191.05
#map_y = 5679754.30
#lon = 10.4530
#lat = 51.0792
#site_name = 'hainich'
#vrt_dir = "vrt_hainich/"

# Wytham
#site_name = 'wytham'
#tile = "h17v03"
#map_x = -92362.62
#map_y = 5757123.82
# Geo. coordinates of origin (?)
#lon = -1.339
#lat = 51.775
#vrt_dir = "vrt_wytham/"

# Somalia
#site_name = 'somalia'
#tile = "h22v08"
#map_x = 5252830.91 #5203067.25
#map_y = 667170.31
# Geo. coordinates of origin (?)
#lon = 47.05
#lat = 6.0
#vrt_dir = "vrt_somalia/"

# Viterbo
site_name = 'viterbo'
tile = "h18v04"
map_x = 987840.2
map_y = 4712492.02
# Geo. coordinates of origin (?)
lon = 12.02656111
lat = 42.38041111
vrt_dir = "vrt_viterbo/"

save_dir = "data/"
year_range = np.arange(2001, 2016)


#US-Ne1
#tile = "h10v04"
#loc = '/home/ucfamc3/DATA/mead/MODIS/MOD09GQ/'
#map_x = -8076086.29
#map_y = 4577355.43

# Determine 250 or 500 meters product
resolution = 500



def init_proj():
    """
    Initialize WGS84 and MODIS sinusoidal projections
    """

    # Create an instance of the SpaialReference class
    # where definitions of coordinate systems
    wgs84 = osr.SpatialReference()

    # Initialize it with WFS84 by EPSG code
    wgs84.ImportFromEPSG(4326)
    #print wgs84

    # from prof. Lewis/Jose
    modis_sinu = osr.SpatialReference() # define the SpatialReference object
    # In this case, we get the projection from a Proj4 string
    modis_sinu.ImportFromProj4( \
                        "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext")

    return wgs84, modis_sinu

#for year in year_range:
#    print year
#    m = modis_data.MODISFinder( tile=tile, year=year, head_loc=loc)
#    m.create_virtual_data()

for year in year_range:
        g = gdal.Open(vrt_dir + "statekm_%d_%s.vrt" % (year, tile))
        g2 = g.GetGeoTransform()

        wgs84, modis_sinu = init_proj()

        # Transform from wgs84 to MODIS
        tx = osr.CoordinateTransformation(wgs84, modis_sinu)
        (x, y, z) = tx.TransformPoint(lon, lat)

        #x = map_x; y = map_y
        pixel_y2 = int(np.round((y - g2[3])/(g2[5])))
        pixel_x2 = int(np.round((x - g2[0])/g2[1]))

        print 'pixel_x, pixel_y = ', pixel_x2, pixel_y2

	if resolution == 500:
		datasets = dict(zip(
		        ['sza', 'saa', 'vza', 'vaa', 'qa', 'b01', 'b02', 'b03', 'b04', 'b05', 'b06', 'b07'],
		        [vrt_dir + "SolarZenith_%d_%s.vrt"%(year, tile),\
                  vrt_dir + "SolarAzimuth_%d_%s.vrt"%(year, tile),\
		         vrt_dir + "SensorZenith_%d_%s.vrt"%(year, tile),\
                  vrt_dir + "SensorAzimuth_%d_%s.vrt"%(year, tile),\
		         vrt_dir + "statekm_%d_%s.vrt"%(year, tile),\
                  vrt_dir + "brdf_%d_%s_b01.vrt"%(year, tile),\
                  vrt_dir + "brdf_%d_%s_b02.vrt"%(year, tile),\
		         vrt_dir + "brdf_%d_%s_b03.vrt"%(year, tile),\
                  vrt_dir + "brdf_%d_%s_b04.vrt"%(year, tile),\
                  vrt_dir + "brdf_%d_%s_b05.vrt"%(year, tile),\
		         vrt_dir + "brdf_%d_%s_b06.vrt"%(year, tile),\
                  vrt_dir + "brdf_%d_%s_b07.vrt"%(year, tile)]))
	else:
		datasets = dict(zip(
		        ['sza', 'saa', 'vza', 'vaa', 'qa', 'b01', 'b02'],
		        [vrt_dir + "SolarZenith_%d_%s.vrt"%(year, tile),\
                  vrt_dir + "SolarAzimuth_%d_%s.vrt"%(year, tile),\
		         vrt_dir + "SensorZenith_%d_%s.vrt"%(year, tile),\
                  vrt_dir + "SensorAzimuth_%d_%s.vrt"%(year, tile),\
		         vrt_dir + "statekm_%d_%s.vrt"%(year, tile),\
                  vrt_dir + "brdf_%d_%s_b01.vrt"%(year, tile),\
                  vrt_dir + "brdf_%d_%s_b02.vrt"%(year, tile)]))


        output = {}
        for ds, fich in datasets.iteritems():
            g = gdal.Open( fich )
            s = g.ReadRaster( xoff=pixel_x2, yoff=pixel_y2, xsize=1, ysize=1,
                         band_list=[b+1 for b in xrange(g.RasterCount)] )
            #pdb.set_trace()
            if ds.find("qa") >= 0:
                x = np.array(struct.unpack(g.RasterCount*"H" , s))
            else:
                x = np.array(struct.unpack(g.RasterCount*"h" , s))

            output[ds] = x

        # QA=200 should be excluded. I think.
        QA_OK=np.array( [ 8, 72, 136, 200, 1288, 2056, 2120, 2184, 2248 ] )
        qa_passer = np.logical_or.reduce( [ output['qa'] == x for x in QA_OK] )

	# reflectance in any band shouldn't be equal or less than zero
        # Except band six which has defect values
	if resolution == 500:
		bb = [1,2,3,4,5,7]
	else:
		bb = [1,2]
	for b in bb:
		ii = np.where(output["b0%d"%b] <= 0.)
		qa_passer[ii] = 0.

        output['qa_passer'] = qa_passer
        [bin(b) for b in QA_OK]

        # Save it for later. I think this can be important even if not too sophisticated.
        # for i in range(len(output['b01'])):
        #         if output['b01'][i] > 1000:
        #                 qa_passer[i] = 0

        doys = np.array([int(g.GetRasterBand(b+1).GetMetadata()['DoY']) for b in xrange(g.RasterCount)])
        
    
	if resolution == 500:
		# header = "DoY QA QA_PASS SZA SAA VZA VAA RAA B01 B02 B03 B04 B05 B06 B07"
		header = "DoY QA QA_PASSER SZA SAA VZA VAA RAA RHO_B1 RHO_B2 RHO_B3 RHO_B4 RHO_B5 RHO_B6 RHO_B7"
		np.savetxt(save_dir + "MODIS_%s_%d.txt" % (site_name, year),\
                   np.c_[ doys, output['qa'], qa_passer, output["sza"]/100., output["saa"]/100.,
		          output["vza"]/100., output["vaa"]/100., (output["saa"]-output["vaa"])/100.,
		          output['b01']/10000., output['b02']/10000.,output['b03']/10000.,output['b04']/10000.,
		          output['b05']/10000., output['b06']/10000., output['b07']/10000.], fmt="%10.4f", header=header )
	else:
		header = "DoY QA QA_PASSER SZA SAA VZA VAA RAA RHO_B1 RHO_B2"
		np.savetxt(save_dir + "MODIS_%s_%d.txt" % (site_name, year),\
                   np.c_[ doys, output['qa'], qa_passer, output["sza"]/100., output["saa"]/100.,
		          output["vza"]/100., output["vaa"]/100., (output["saa"]-output["vaa"])/100.,
		          output['b01']/10000., output['b02']/10000.], fmt="%10.4f", header=header )

print 'Read MODIS data is done!!!'
