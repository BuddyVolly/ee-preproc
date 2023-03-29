import ee

def apply(images, options):

    def mapMeanSpace(image):

        reducer = ee.Reducer.mean()
        kernel = ee.Kernel.square(options['radius'], options['units'])
        mean = image.reduceNeighborhood(reducer, kernel).rename(meanBand)
        ratio = image.divide(mean).rename(ratioBand)
        return image.addBands(mean).addBands(ratio).copyProperties(image)


    bands = images.first().bandNames()
    meanBand = bands.map(lambda b: ee.String(b).cat('_mean'))
    ratioBand = bands.map(lambda b: ee.String(b).cat('_ratio'))

    meanSpace = images.map(mapMeanSpace)

    def mtDespeckSingle(image):

        meanSpace2 = ee.ImageCollection(meanSpace).select(ratioBand)
        b = image.select(meanBand)
        filtered = b.multiply(meanSpace2.sum()).divide(meanSpace2.count()).rename(bands).select(['HH', 'HV'])

        return image.addBands(filtered, None, True).select(bands)

    return meanSpace.map(mtDespeckSingle);
