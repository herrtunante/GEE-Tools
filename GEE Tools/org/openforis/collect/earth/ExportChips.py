#!/usr/bin/env python
'''
Created on Jan 23, 2017

@author: Alfonso Sanchez-Paus
'''

import ee
import ee.mapclient
from ee.ee_date import Date
from ee.batch import Export


ee.Initialize()

plots = ee.FeatureCollection('ft:1QMa9xE6ePmEe1gaSShxBCQqVvzKheX30mr8qT64D')

def addBufferChip(feature):
    """A function to add a buffer of 35 meters around the plot location."""
    return feature.buffer(1500)

def addBufferIndices(feature):
    """A function to add a buffer of 35 meters around the plot location."""
    return feature.buffer(35)

plots_for_chips = plots.map(lambda f: f.buffer(2000))

plots_for_indices = plots.map(lambda f: f.buffer(35))

#Where is the work carried out ?
country_name = 'CubaImprovedV1';

#Do you want to export the roi OR the convex hull of the ROI ? 
export_roi = "yes";
export_hul = "no";

#Which years to look at ? 
year_start = 2015;
year_end   = 2016;

# Which days of the year to use to search the LANDSAT collection? 
julianDayStart = 1;
julianDayEnd   = 365;

#Which cloud cover is accepted?
cloudcovthres = 5;

#ee.mapclient.addToMap(small_counties, {'color': '900000'})

#//////////////////// Create functions to harmonize band names
my_band_names = ['B2','B3','B4','B5','B6','B7', 'B8','B10', 'B11'];

def harmonize_l57(image):
    return image.select(['B1','B2','B3','B4','B5','B7','B8','B6_VCID_1','B6_VCID_2'], my_band_names)

def harmonize_l8(image):
    return image.select(['B2','B3','B4','B5','B6','B7','B8','B10','B11'], my_band_names)

sentinelVisBands = ['B8_mean','B4_mean','B3_mean'];
l7VisBands = ['B5_mean', 'B6_mean', 'B4_mean'];

def cloudMask_jrc(im):
    """Implement Sentinel cloud mask
    Opaque and cirrus cloud masks cause bits 10 and 11 in QA60 to be set,
    so values less than 1024 are cloud-free
    """
    mask = ee.Image(0).where(im.select('QA60').gte(1024), 1).Not();
    return im.updateMask(mask);
#    cloud1 = ee.Image(0).where(im.select('QA60').gte(1024), 1)
#     cloud2 = ee.Image(0).where(im.select('B2').gte(1200), 1)
#     cloud2b = cloud2.focal_min(200,'circle','meters')
#     #compute the n of cloudfree observations
#     clouds=((cloud1.eq(1)).Or(cloud2b.eq(1)))
#     im = im.addBands(clouds.eq(0).select([0],["cloudfree"])); # and add a band with this information
#     return im.updateMask(clouds.eq(0));

def generateChips(chip):
    # Which years ?
    for year in range(year_start,year_end):
        print(year);
        
        tmin = str(year) + '-11-01';
        tmax = str(year+1) + '-04-30';
    
        #     /////////////////////////////////
        #     // Combine FILTERS           ////
        #     /////////////////////////////////
        #Create a location filter
        locFilter = ee.Filter.geometry(chip.geometry() )
    
        # Create a range of date to filter imagery to use
        dateFilter = ee.Filter.date( Date(tmin), Date(tmax) )
    
        dateLocFilter = ee.Filter.And(locFilter,dateFilter)
    
        #Which days of the year to use ? Which cloud cover is accepted?
        doyFilter = ee.Filter().calendarRange(julianDayStart,julianDayEnd);
        cloudCoverFilter = ee.Filter().lt('CLOUD_COVER',cloudcovthres);
        
        #Combine the filters 
        combinedyFilter = ee.Filter.And(dateLocFilter,doyFilter, cloudCoverFilter);
        
        
        #     /////////////////////////////////
        #     // Select Imagery | Filters  ////
        #     /////////////////////////////////
        #     // Look into the archives and select the imagery following the filters
        inputl5 = ee.ImageCollection('LANDSAT/LT5_L1T_TOA').filter(combinedyFilter);
        inputl7 = ee.ImageCollection('LANDSAT/LE7_L1T_TOA').filter(combinedyFilter);
        inputl8 = ee.ImageCollection('LANDSAT/LC8_L1T_TOA').filter(combinedyFilter);
    
        #     /////////////////////////////////////////
        #     //// Check availability of archives /////
        #     /////////////////////////////////////////
        count_l5 = inputl5.aggregate_count('system:index')
        count_l7 = inputl7.aggregate_count('system:index')
        count_l8 = inputl8.aggregate_count('system:index')
        
        
        #     /////////////////////////////////////////
        #     //// Normalize band names           /////
        #     /////////////////////////////////////////
        l5 = inputl5.map(harmonize_l57)
        l7 = inputl7.map(harmonize_l57)
        l8 = inputl8.map(harmonize_l8)
    
        #     /////////////////////////////////////////
        #     //// Merge 3 collections into one   /////
        #     /////////////////////////////////////////
        merged = ee.ImageCollection(l5.merge(l7).merge(l8))
       
        #     /////////////////////////////////////////
        #     //// Create median values           /////
        #     /////////////////////////////////////////
        #median_al = merged.reduce( ee.Reducer.intervalMean(25,75) )
                  
        idChip = chip.get('id').getInfo() 
        
        #     ////////////////////////////////////////////////////////////////////////////////////////////////
        #     ///////////////////// EXPORT RESULTS USING THE DOWNLOAD GRID YOU SPECIFIED         /////////////
        #     ////////////////////////////////////////////////////////////////////////////////////////////////
        if export_roi == 'yes' and year >= 2015:
            # Call the sentinel 2 collection
            inputs2 = ee.ImageCollection('COPERNICUS/S2').filter(dateLocFilter).filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', 30) )
        
            # Mask out the clouds
            s2clouds=inputs2.map(cloudMask_jrc);
        
            # Take only the median observations and select 10m bands
            median_s2 = s2clouds.select(['B2','B3','B4','B8']).reduce(ee.Reducer.intervalMean(25,75));
        
            s2task = ee.batch.Export.image.toDrive(
                                 median_s2.select(sentinelVisBands).clip(chip.geometry() ), 
                                 description = 'median_roi_clip_s2_bb_'+country_name ,#+ "_" + str( chip.get('id') ), 
                                 scale = 10, 
                                 folder = country_name+ '_Sentinel_Chips_' + str(year),
                                 fileNamePrefix = idChip
                                 )
            s2task.start()
            
            print ( s2task.status()['id'] )

    
        if export_roi == 'yes' :
            landsatTask = ee.batch.Export.image.toDrive( 
                                 panSharpen(merged).select(l7VisBands).clip(chip.geometry() ), 
                                 description= '_median_roi_clip_'+ str(year) +'_bb_'+ country_name , # "_" + str( chip.get('id') ), 
                                 scale= 30, 
                                 folder = country_name+'_Landsat_Chips_' + str(year) , 
                                 fileNamePrefix = idChip 
                                 )
            landsatTask.start();
            
            print ( landsatTask.status()['id'] )
            

def panSharpen( landsatHarmonizedcollection ):

    # Mask clouds by mapping the cloudMask function over the collection.
    # This will add a cloud score band called 'cloud' to every image.
    collection = landsatHarmonizedcollection.map(lambda f: ee.Algorithms.Landsat.simpleCloudScore(f) );

    # Convert the collection to an array.
    array = collection.toArray();

    # Label of the axes.
    imageAxis = 0;
    bandAxis = 1;

    # Get the cloud slice and the bands of interest.
    bands = array.arraySlice(bandAxis, 0, len(my_band_names))
    clouds = array.arraySlice(bandAxis, len(my_band_names))

    # Sort by cloudiness.
    sortedArray = bands.arraySort(clouds)

    # Get the least cloudy images, 20% of the total.
    numImages = sortedArray.arrayLength(imageAxis).multiply(0.2).int()
    leastCloudy = sortedArray.arraySlice(imageAxis, 0, numImages)

    # Get the mean of the least cloudy images by reducing along the image axis.
    mean = leastCloudy.arrayReduce(reducer= ee.Reducer.mean(), axes = [imageAxis])
    
    bandNames = ee.List(my_band_names);
    
    # Turn the reduced array image into a multi-band image for display.
    meanImage = mean.arrayProject([bandAxis]).arrayFlatten([bandNames])
    # Convert the RGB bands to the HSV color space.
    hsv = meanImage.select( l7VisBands ).rgbToHsv()

    # Swap in the panchromatic band and convert back to RGB.
    sharpened = ee.Image.cat([hsv.select('hue'), hsv.select('saturation'), meanImage.select('B8')]).hsvToRgb()

    return sharpened


#plots_for_chips.map( generateChips );

chipList = plots_for_chips.toList(100)
size = chipList.size().getInfo()

for i in range(0, size):
    generateChips( ee.Feature( chipList.get(i) ) )

