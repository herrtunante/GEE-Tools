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
country_name = 'CubaxxxFinal';

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
my_band_names = ['B3','B4','B5','B6','B7'];

def harmonize_l57(image):
    return image.select(['B2','B3', 'B4', 'B5', 'B7'], my_band_names)

def harmonize_l8(image):
    return image.select(['B3','B4', 'B5', 'B6', 'B7'], my_band_names)

sentinelVisBands = ['B4_mean','B3_mean','B2_mean'];

def cloudMask_jrc(im):
    """Implement Sentinel cloud mask
    Opaque and cirrus cloud masks cause bits 10 and 11 in QA60 to be set,
    so values less than 1024 are cloud-free
    """
    cloud1 = ee.Image(0).where(im.select('QA60').gte(1024), 1)
    cloud2 = ee.Image(0).where(im.select('B2').gte(1200), 1)
    cloud2b = cloud2.focal_min(200,'circle','meters')
    #compute the n of cloudfree observations
    clouds=((cloud1.eq(1)).Or(cloud2b.eq(1)))
    im = im.addBands(clouds.eq(0).select([0],["cloudfree"])); # and add a band with this information
    return im.updateMask(clouds.eq(0));

l7VisBands = ['B5_mean', 'B6_mean', 'B4_mean'];


def generateChips(chip):
    # Which years ?
    for year in range(year_start,year_end):
        print(year);
        
        tmin = str(year) + '-01-01';
        tmax = str(year) + '-12-31';
    
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
        inputl5 = ee.ImageCollection('LT5_L1T_TOA').filter(combinedyFilter);
        inputl7 = ee.ImageCollection('LE7_L1T_TOA').filter(combinedyFilter);
        inputl8 = ee.ImageCollection('LC8_L1T_TOA').filter(combinedyFilter);
    
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
        median_al = merged.reduce( ee.Reducer.intervalMean(25,75) )
        
    #     //Changed addToMap() it became a deprecated method  since September 2016, instead use Map.addLayer() 
        
    #     Map.addLayer(median_al.clipToCollection(roi),{bands: l7VisBands, max: 0.3}, country_name + '_median_roi_'+year,false);
    #     ////////////////////////////////////////////////////////////////////////////////////////////////
    #     ///////////////////// EXPORT RESULTS USING THE COUNTRY SHAPE YOU SPECIFIED        /////////////
    #     ////////////////////////////////////////////////////////////////////////////////////////////////
    
        # Call the sentinel 2 collection
        inputs2 = ee.ImageCollection('COPERNICUS/S2').filter(dateLocFilter);
        
        # Mask out the clouds
        s2clouds=inputs2.map(cloudMask_jrc);
        
        # Take only the median observations and select 10m bands
        median_s2 = s2clouds.select(['B2','B3','B4','B8']).reduce(ee.Reducer.intervalMean(25,75));
        
        
        idChip = chip.get('id').getInfo() 
        print (idChip)
        
    #     ////////////////////////////////////////////////////////////////////////////////////////////////
    #     ///////////////////// EXPORT RESULTS USING THE DOWNLOAD GRID YOU SPECIFIED         /////////////
    #     ////////////////////////////////////////////////////////////////////////////////////////////////
        if export_roi == 'yes' and year >= 2015:
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
                                 median_al.select(l7VisBands).clip(chip.geometry() ), 
                                 description= '_median_roi_clip_'+ str(year) +'_bb_'+ country_name , # "_" + str( chip.get('id') ), 
                                 scale= 30, 
                                 folder = country_name+'_Landsat_Chips_' + str(year) , 
                                 fileNamePrefix = idChip 
                                 )
            landsatTask.start();
            
            print ( landsatTask.status()['id'] )
            


#plots_for_chips.map( generateChips );

chipList = plots_for_chips.toList(100)
size = chipList.size().getInfo()


for i in range(0, size):
    generateChips( ee.Feature( chipList.get(i) ) )

