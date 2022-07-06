# -*- coding:UTF-8 -*-

# The script allows to  prolongate supplementary contour lines. Contour lines with z value are used as an input data.

# activate libraries 
import arcpy
from arcpy.sa import *
import arcpy.cartography as CA



def contour_prolongation(input_contours, 
                         height_field, 
                         intermediate_raster_cellsize, 
                         output_contours):
    # Description of the contour lines features
    desc = arcpy.Describe(input_contours)
    xmin = desc.extent.XMin
    xmax = desc.extent.XMax
    ymin = desc.extent.YMin
    ymax = desc.extent.YMax
    # extent of the file with contour lines
    arcpy.env.extent = desc.extent

    # X, Y dimensions, full extent
    Xext = xmax - xmin
    Yext = ymax - ymin
    full_ext = round((Xext + Yext), 0)

    # a cellsize for intermediate rasters 
    cellSize = 3

    arcpy.AddMessage("Creation of the rasters with constant value")
    # Creation of 2 rasters with constant value. 
    # The value is an extremely big number of euclidean distance (X + Y extents)
    minVal_1 = CreateConstantRaster(full_ext, "FLOAT", cellSize) 
    minVal_1.save("minVal_1")

    # conversion from rasters to arrays 
    # Here will be stored the first minimum value after loop (where comparison of euclidean distances will be performed)
    minVal_1_numpy_array = arcpy.RasterToNumPyArray(minVal_1) 

    # the same for the second raster 
    minVal_2 = CreateConstantRaster(full_ext, "FLOAT", cellSize) 
    minVal_2.save("minVal_2")

    # Here will be stored the second minimum value after loop
    minVal_2_numpy_array = arcpy.RasterToNumPyArray(minVal_2) 

    arcpy.AddMessage("Calculation of the contour interval")
    # the name of the field where z values are stored
    Hvalues = [row[0] for row in arcpy.da.SearchCursor(input_contours, height_field)]
    # find the z-values of contour lines
    uniqueValues = set(Hvalues)
    # unique z-values
    Hvalues = list(uniqueValues) 

    heights = Hvalues
    heights_unique = list(set(heights))
    # ascending sorting of unique z-values
    heights_unique.sort()
    arcpy.AddMessage("Unique z-values: {0}".format(heights_unique))
    print(heights_unique)

    # the list where contour interval between neighboring contour lines will be stored
    steps = []
    # find difference between current and following z-value. 
    # The number is the contour interval 
    for i in range(len(heights_unique) - 1):
        step = heights_unique[i + 1] - heights_unique[i]
        # add to the list
        steps.append(step)
    # print(steps)
    steps_unique = list(set(steps))
    # ascending sorting of contour intervals
    steps_unique.sort()
    # the last value corresponds to the contour interval between intermediate lines
    contour_step = steps_unique[-1]
    arcpy.AddMessage("The contour interval is {0} m".format(contour_step))

    arcpy.AddMessage("Creation of the new field to divide contour lines into intermediate and supplementary")
    arcpy.AddField_management(input_contours, "line_type", "TEXT", field_alias = "line_type")
    fields = [height_field, "line_type"]
    with arcpy.da.UpdateCursor(input_contours, fields) as cursor:
        for row in cursor:
            # if remainder of the division of z-value into contour interval is equal to 0,
            # then it is the intermediate contour line
            if(row[0] % contour_step == 0):
                row[1] = "intermediate"
            else:
                # if remainder of the division of z-value into contour interval is not equal to 0, 
                # then it is the supplementary contour line
                row[1] = "supplementary"
            cursor.updateRow(row)

    # Make a layer from input layer
    arcpy.MakeFeatureLayer_management(input_contours, "lyr") 

    # Selection of only intermediate lines. It allows to find middle lines (future supplementary lines) 
    arcpy.SelectLayerByAttribute_management("lyr", "NEW_SELECTION", "\"line_type\" = 'intermediate'")


    arcpy.AddMessage("Euclidean distance calculation")
    # Cells will be marked 'no data' if their value is bigger than maximum distance (treshhold)
    maxDistance = full_ext
    # Euclidean distance are used for creation of middle lines between 2 intermediate contours
    outEucDistance = EucDistance("lyr", maxDistance, cellSize) 
    # save output 
    outEucDistance.save("EuDist")

    subset_layer = "lyr"
    # rasters extent
    nrows = outEucDistance.height
    ncols = outEucDistance.width

    # List of unique z-values for the subset layer
    subHvalues = [row[0] for row in arcpy.da.SearchCursor(subset_layer, height_field)]
    uniqueSubHValues = set(subHvalues)
    subHvalues = list(uniqueSubHValues) 
    arcpy.AddMessage("Unique intermediate z-values: {0}".format(subHvalues))
    print(subHvalues)

    arcpy.AddMessage("Comparison new euclidean distance with primary one")
    # In the loop comparison between new and primary euclidean distances is provided
    for i in subHvalues:
        # Select objects within subset layer
        arcpy.SelectLayerByAttribute_management("lyr", "NEW_SELECTION", "%s = %s" % (height_field, i)) 
        # %s - row format, then substitution
        # new euclidean distance
        outEucDistance = EucDistance("lyr", maxDistance, cellSize)  
        EucDist_numpy = arcpy.RasterToNumPyArray(outEucDistance)
        # Comparison new euclidean distance with primary one
        for row in range(nrows):
            for col in range(ncols):
                # if euclidean distance (ED) for current point is less than ED in constant value raster (1),
                # then rewrite new value to (1). Then write it to constant value raster (2) 
                if(EucDist_numpy[row, col] <= minVal_1_numpy_array[row, col]): 
                    minVal_2_numpy_array[row, col] = minVal_1_numpy_array[row, col]
                    minVal_1_numpy_array[row, col] = EucDist_numpy[row, col] 
                # if ED for current point is more than ED in constant value raster (1)
                # and less than ED in constant value raster (2), then new value  
                # is written to the constant value raster
                elif(EucDist_numpy[row, col] > minVal_1_numpy_array[row, col] and 
                EucDist_numpy[row, col] < minVal_2_numpy_array[row, col]):
                    minVal_2_numpy_array[row, col] = EucDist_numpy[row, col]

    arcpy.AddMessage("Middle lines creation")                
    # Rewrite constant value rasters with new euclidean distance
    # convert arrays to rasters
    minVal_1_new = arcpy.NumPyArrayToRaster(minVal_1_numpy_array, arcpy.Point(xmin, ymin), cellSize, cellSize)
    minVal_1_new.save("minVal_1_new")

    minVal_2_new = arcpy.NumPyArrayToRaster(minVal_2_numpy_array, arcpy.Point(xmin, ymin), cellSize, cellSize)
    minVal_2_new.save("minVal_2_new")
            
    # difference between 1 and 2 constant value rasters: cells with values close to 0
    # mark middle lines between 2 intermediate lines
    diff_array = minVal_1_numpy_array - minVal_2_numpy_array
    diff = arcpy.NumPyArrayToRaster(diff_array, arcpy.Point(xmin, ymin), cellSize, cellSize)
    # Convert to integer
    diff_int = Int(diff)  

    # Set null to separate middle lines from background
    # the value was choseen empirically
    set0 = SetNull(diff_int, diff_int, '"VALUE" < -5')
    set0_int = Int(set0)

    arcpy.AddMessage("Middle line processing") 
    # Raster thinning
    # background value are marked as nodata; rounding angles
    ThinRaster = Thin(set0_int, background_value = "NODATA", filter = "NO_FILTER", corners = "ROUND") 
    # Convert raster to lines
    midlines = "midlines"
    # background values = 0; simplify geometry
    arcpy.RasterToPolyline_conversion(ThinRaster, midlines, background_value = "ZERO", 
    simplify = "SIMPLIFY", raster_field = "VALUE") 

    # Lines form optimization 
    # Values were choseen empirically
    midlines_smpl = "midlines_smpl"
    CA.SimplifyLine(midlines, midlines_smpl, "BEND_SIMPLIFY", tolerance = 25) 
    midlines_smooth = "midlines_smooth"
    CA.SmoothLine(midlines_smpl, midlines_smooth, "PAEK", tolerance = 5)

    arcpy.AddMessage("Selection of middle lines") 
    arcpy.MakeFeatureLayer_management(input_contours, "suppl_lines")  
    # Choose only supplementary contour lines 
    suppl_lines = arcpy.SelectLayerByAttribute_management("suppl_lines", "NEW_SELECTION",
    "\"line_type\" = 'supplementary'")

    arcpy.MakeFeatureLayer_management(midlines_smooth, "midlines") 
    # Choose middle lines which intersect primary supplementary lines
    midl_intersect = arcpy.SelectLayerByLocation_management("midlines", "INTERSECT", suppl_lines)

    arcpy.AddMessage("Intersect middle lines with supplementary lines") 
    # Intersect middle lines with supplementary lines, points will be as output format
    for_inters = [suppl_lines, midl_intersect]
    inters_point = "intersection_point"
    arcpy.Intersect_analysis(for_inters, inters_point, "ALL", output_type = "point")

    # Intersection points converts from multipart object to singlepart
    inters_point_single = "intersection_point_single"
    arcpy.MultipartToSinglepart_management(inters_point, inters_point_single)

    # Split middle lines at intersection points
    split_lns = "splitted_lines"
    arcpy.SplitLineAtPoint_management(midl_intersect, inters_point_single, split_lns, 
                                      search_radius= "1 Meters")

    arcpy.AddMessage("Spatial join of intersection points and middle lines")
    # Spatial join of intersection points and selected middle lines to get field with z-values
    sp_join_midl_intersect = "sp_join_midl_intersect"
    arcpy.SpatialJoin_analysis(split_lns, inters_point_single, match_option= "INTERSECT", 
                               out_feature_class = sp_join_midl_intersect, join_operation = "JOIN_ONE_TO_ONE", 
                               join_type = "KEEP_ALL", search_radius= "1 Meters")

    arcpy.AddMessage("Processing of every segment of middle lines")
    # Selection of necessary layers and fields
    arcpy.MakeFeatureLayer_management(sp_join_midl_intersect, "split_lns")
    arcpy.MakeFeatureLayer_management(inters_point_single, "intersect_pts")
    arcpy.da.SearchCursor(inters_point_single, "FID_new_lines_full")
    midlines_split = arcpy.da.SearchCursor(sp_join_midl_intersect, "OBJECTID")

    arcpy.AddMessage("Intersection point selection")
    for i in midlines_split:
        # Selection of every fragment of middle lines
        selection = "OBJECTID = {}".format(i[0])       
        lns = arcpy.SelectLayerByAttribute_management("split_lns", "NEW_SELECTION", selection)
        # Selection of points which intersect with fragments
        sel = arcpy.SelectLayerByLocation_management("intersect_pts", "INTERSECT", lns, 
        search_distance = "1 Meters", selection_type = "NEW_SELECTION")
        # Count selected points
        counts = int(arcpy.GetCount_management(sel)[0])
        
        # If there is 1 point then save the point to a new layer. If there are 2 points, then find ID  
        # of primary supplementary lines. If IDs of 2 points are different, then save the line
        if(counts == 1):
            # Check if the layer exists. If doesn't exist create a layer, if it exists add the line 
            if arcpy.Exists("selected_lns"):
                arcpy.Append_management(lns, "selected_lns", "NO_TEST")
            else:
                arcpy.FeatureClassToFeatureClass_conversion(lns, arcpy.env.workspace, "selected_lns")
        if(counts == 2):
            # Select middle lines, if intersection points (middle lines with supplementary lines) have different IDs., 
            # i.e. different supplementary lines which are joined
            ids = [row[0] for row in arcpy.da.SearchCursor(sel, "OBJECTID")]
            if(ids[0] != ids[1]):
                if arcpy.Exists("selected_lns"):
                    arcpy.Append_management(lns, "selected_lns", "NO_TEST")
                else:
                    arcpy.FeatureClassToFeatureClass_conversion(lns, arcpy.env.workspace, "selected_lns")
                    


    arcpy.AddMessage("Creation of Voronoy diagram")
    # Create Voronoy diagram for all input contours. It will be used as a mask for clipping middle lines.
    # Firstly, it is necessary to densify vertices in contours, convert lines to points. 
    # Copy layer because a new layer is not created during densify operation
    full_lns_copy = "full_lns_copy"
    arcpy.CopyFeatures_management(input_contours, full_lns_copy)
    arcpy.Densify_edit(full_lns_copy, "DISTANCE")

    pts_from_all_lines = "pts_from_all_lines"
    arcpy.FeatureVerticesToPoints_management(full_lns_copy, pts_from_all_lines, "ALL")
    
    # Creation of Voronoy diagram
    thiess_poly = "thiess_poly"
    arcpy.CreateThiessenPolygons_analysis(pts_from_all_lines, thiess_poly, fields_to_copy = "ALL")
    arcpy.MakeFeatureLayer_management(thiess_poly, "thiess_poly")
    # Choose polygons which correspond to supplementary contours
    suppl_poly = arcpy.SelectLayerByAttribute_management("thiess_poly", "NEW_SELECTION",
    "\"line_type\" = 'supplementary'")
    # Clip middle lines along selected polygons 
    clipped_lns_by_polyg = "clipped_lns_by_polyg"
    arcpy.Erase_analysis("selected_lns", suppl_poly, clipped_lns_by_polyg)
    arcpy.Delete_management(full_lns_copy)
    arcpy.Delete_management(pts_from_all_lines)
    arcpy.Delete_management(thiess_poly)
    arcpy.Delete_management("selected_lns")

    arcpy.AddMessage("Merge middli lines with supplemenatry lines")
    # Join middle lines with supplementary contour lines
    full_lines_copy = 'full_lines_copy'
    arcpy.CopyFeatures_management(input_contours, full_lines_copy)
    fieldMappings = arcpy.FieldMappings()
    # Add all attributes from 2 layers
    fieldMappings.addTable(full_lines_copy)
    fieldMappings.addTable(clipped_lns_by_polyg)
    fldMap = arcpy.FieldMap()
    fldMap.addInputField(full_lines_copy, height_field)
    fieldMappings.addFieldMap(fldMap)

    # Delete all attributes except z-value and ID
    for field in fieldMappings.fields:
        if field.name not in ["OBJECTID", height_field]:
            fieldMappings.removeFieldMap(fieldMappings.findFieldMapIndex(field.name))

    joined_contours = "joined_contours"
    arcpy.Merge_management([full_lines_copy, clipped_lns_by_polyg], joined_contours, fieldMappings)
    # Delete duplicated field with z-value
    arcpy.DeleteField_management("joined_contours", "Habs_1")

    # Create end points of contours for joining
    end_pts_from_lns = "end_pts_from_lns"
    arcpy.FeatureVerticesToPoints_management(joined_contours, end_pts_from_lns, "BOTH_ENDS")
    # find the nearest countours' ends whithin the layer
    # input и near - end points
    arcpy.Near_analysis(end_pts_from_lns, end_pts_from_lns, search_radius = "500 Meters")

    # if in the point layer OBJECTID_1[0] = NEAR_FID[1] и OBJECTID_1[1] = NEAR_FID[0], 
    # join these points
    arcpy.AddMessage("Find point pairs:")
    fields = ["OBJECTID", "NEAR_FID", height_field, "NEAR_DIST"]
    fids_list = []
    # list with points which are NOT appropriate 
    skip_list = []
    for k in arcpy.da.SearchCursor("end_pts_from_lns", fields):
        # current row, check if the distance is not equal to 0
        if(k[0] in skip_list or k[3] < 1e-7):
            continue
        # Choose an appropriate point
        # Condition: id of the nearest point (NEAR_FID) = id of the current point (OBJECTID_1). 
        # If it is true on the contrary, then these 2 points are nearest to each other. They should be joined
        where_clause = "OBJECTID = %s" % (k[1]) 
        for j in arcpy.da.SearchCursor("end_pts_from_lns", fields, where_clause):
            if(j[1] == k[0] and k[2] == j[2]):  # check the equality of z-values
                fids_list.append([k[0], j[0], k[2]])
                # to avoid duplicating of point pairs skip the second one
                skip_list.append(j[0])
    # print(fids_list)

    arcpy.AddMessage("Convert points to lines and write z-values")
    pts_to_lns = "pts_to_lns"
    arcpy.CopyFeatures_management(joined_contours, "joined_cont_copy")
    arcpy.MakeFeatureLayer_management(end_pts_from_lns, "end_pts_from_lns")
    for i in fids_list:
        # Round brackets are necessary for selection 
        pts_ids = "OBJECTID IN " + str([i[0], i[1]]).replace('[','(').replace(']',')')
        selected_pairs = arcpy.SelectLayerByAttribute_management("end_pts_from_lns", "NEW_SELECTION", pts_ids)
        arcpy.PointsToLine_management(selected_pairs, "pts_to_lns")
        # Create field to write z-value
        arcpy.AddField_management(pts_to_lns, height_field, "FLOAT")
        addH = arcpy.da.UpdateCursor("pts_to_lns", ["OBJECTID", "Shape", "Shape_Length", height_field])
        # write z-value to the field Habs
        for h in addH:
            h[3] = i[2]
            addH.updateRow(h)
        del h, addH

        # Add every line to the main layer
        arcpy.Append_management("pts_to_lns", "joined_cont_copy", "NO_TEST")

    arcpy.AddMessage("Join input supplementary lines with segments of middle lines")
    # Merge primary contours with line segments by z-value 
    dissolved_final_lns = "dissolved_final_lns"
    arcpy.Dissolve_management("joined_cont_copy", dissolved_final_lns, height_field)

    # Convert contour lines from multipart to singlepart objects
    final_lines = output_contours
    arcpy.MultipartToSinglepart_management(dissolved_final_lns, final_lines)

    # Delete unnecessary attribute fields 
    arcpy.DeleteField_management(output_contours, "ORIG_FID")

    # Create shape-file
    # arcpy.FeatureClassToFeatureClass_conversion(final_lines, , "finally_processed_lns.shp") 


if __name__ == '__main__':
    # division of the calculation execution to several processes for acceleration
    arcpy.env.parallelProcessingFactor = "0"
    # activate the file overwriting
    arcpy.env.overwriteOutput = True
    args = tuple(arcpy.GetParameterAsText(i) for i in range(arcpy.GetArgumentCount()))
    contour_prolongation(*args)
